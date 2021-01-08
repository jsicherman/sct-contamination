#' CPM
#' Compute counts per million for a sparse matrix
#' 
#' @param y l `matrix` of counts
#' @param log logical, if `TRUE` then `log2` values are returned
#'
#' @seealso @source{edgeR::cpm}
cpm <- function(y, log = F) {
  require(Matrix)
  
  sf <- Matrix::colSums(y) / 1e6
  
  if(class(y) == 'dgCMatrix') {
    sep <- y@p
    sep <- sep[-1] - sep[-length(sep)]
    j <- S4Vectors::Rle(1:length(sep), sep)
    y@x <- y@x / sf[as.integer(j)]
    
    if(log)
      y@x <- log2(y@x + 1)
  } else
    stop(paste('CPM for', class(y)[1], 'not supported.'))
  
  y
}

#' Make a Seurat data object ready for these analyses
#'
#' @param counts A dgCMatrix of counts
#' @param metadata A data frame of metadata. Must have a "cluster" column
#'
#' @return The experiment in Seurat format with cpm in the data slot and identities set to the cluster column from the meta
as_Seurat <- function(counts, metadata) {
  require(Seurat)
  
  seurat <- CreateSeuratObject(counts)
  seurat@assays$RNA@data <- cpm(seurat@assays$RNA@counts, T)
  seurat@meta.data <- cbind(seurat@meta.data, metadata)
  Idents(seurat) <- seurat$cluster
  seurat
}

#' Find cell type markers and aggregate reference data by cell type
#'
#' @param ref.data A Seurat object from @method{as_Seurat}
#' @param min.pct The minimum fraction of cells that must express a gene for it to be tested
#' @param min.diff.pct The minimum difference in fraction of detection for gene to be detected
#' @param logfc.threshold The minimum logfc to call a marker gene
#' @param return.thresh The maximum p-value of markers
#'
#' @return A list with marker genes in "markers" and cluster-wise gene expression medians in "collapsed"
simplify_reference <- function(ref, min.pct = 0.1, min.diff.pct = 0.4, return.thresh = 0.01, logfc.threshold = 1,
                               cores = 4) {
  require(Seurat)
  require(MAST)
  require(dplyr)
  
  cores.orig = getOption('mc.cores')
  options(mc.cores = cores)
  markers <- FindAllMarkers(ref, min.pct = min.pct, min.diff.pct = min.diff.pct,
                            logfc.threshold = logfc.threshold, return.thresh = return.thresh, test.use = 'MAST')
  options(mc.cores = cores.orig)
  
  list(markers = markers,
       collapsed = do.call(cbind, lapply(levels(Idents(ref)), function(cluster) {
         subdata <- subset(ref, idents = cluster)@assays$RNA@data
         trim <- floor(0.25 * ncol(subdata))
         
         apply(subdata, 1, function(x) {
           mean(sort(x)[-c(1:trim, (length(x)-trim):length(x))])
         })
       })) %>% `colnames<-`(levels(Idents(ref))))
}

#' Map samples to a collapsed reference
#'
#' @param test A dgCMatrix of CPMs
#' @param ref.collapsed A collapsed reference dataset from @method{simplify_reference}
#' @param markers.perc The percentage of markers to use in each iteration, relative to the smallest cluster
#' @param n.iter The number of bootstrap iterations
#' @param n.cores The number of CPU cores
#' @param collapse Whether or not to collapse into a cumulative correlation
map_samples <- function(test, ref.collapsed, clusters = colnames(ref.collapsed$collapsed),
                        markers.perc = 0.6, n.iter = 1e3, n.cores = 4, collapse = T, cor.method = 'pearson') {
  require(dplyr)
  require(parallel)
  
  ref <- ref.collapsed$collapsed[, clusters]
  markers <- ref.collapsed$markers[ref.collapsed$markers$cluster %in% clusters, ]
  
  markers.use <- lapply(colnames(ref), function(x) {
    intersect(rownames(ref),
              intersect(rownames(test),
                        rownames(markers)[markers$cluster == x]))
  }) %>% setNames(colnames(ref))
  
  n.markers <- min(sapply(markers.use, function(cluster) floor(markers.perc * length(cluster))))
  message(paste('Using', n.markers, 'markers per cluster.'))
  
  cores.original <- getOption('mc.cores')
  options(mc.cores = n.cores)
  
  remapped <- mclapply(1:n.iter, function(x) {
    markers.iter <- unlist(lapply(markers.use, sample, n.markers))
    
    do.call(rbind, lapply(1:ncol(test), function(sample) {
      cor(test[markers.iter, sample], ref[markers.iter, ], method = cor.method)
    })) %>% `colnames<-`(colnames(ref)) %>% `rownames<-`(colnames(test)) -> correlated
    
    list(correlation = correlated,
         mapping = setNames(colnames(correlated)[apply(correlated, 1, which.max) %>% lapply(max, 1) %>% unlist],
                            rownames(correlated)))
  })
  
  options(mc.cores = cores.original)
  
  if(!collapse)
    return(remapped)
  
  mapping <- matrix(0, ncol(test), ncol(ref)) %>%
    `rownames<-`(colnames(test)) %>% `colnames<-`(colnames(ref))
  
  for(i in 1:length(remapped)) {
    for(j in 1:nrow(mapping)) {
      mapping[j, remapped[[i]]$mapping[j]] <- mapping[j, remapped[[i]]$mapping[j]] + 1
    }
  }
  
  stdev.matrices <- function(matrices) {
    apply(array(unlist(matrices), c(dim(matrices[[1]]), length(matrices))), 1:2, sd) %>%
      `rownames<-`(rownames(matrices[[1]])) %>% `colnames<-`(colnames(matrices[[1]]))
  }
  
  list(cum.correlation = Reduce('+', lapply(remapped, '[[', 1)),
       max.correlation = Reduce(pmax, lapply(remapped, '[[', 1)),
       min.correlation = Reduce(pmin, lapply(remapped, '[[', 1)),
       stdev.correlation = stdev.matrices(lapply(remapped, '[[', 1)),
       mapping = mapping,
       iterations = n.iter)
}

#' Normalize covariates based on an empirical prior
#'
#' @param cum.cor The scaled output of @method{map_samples(..., collapse = T)$cum.correlation}
#' @param clusters The cluster definitions (corresponding to the rows of @param{cum.cor})
#' @param prior The prior from @method{(generate_prior)}
#' @param remappings A named character of remappings to ensure that cluster names in @param{clusters} and @param{prior} are normalized
#' @param shift Whether or not to push the mean of the on-target cell type to 0 (and everything with it)
#'
#' @return Covariates from @param{cum.cor} downscaled by the appropriate @param{prior}
normalize_covariates <- function(cum.cor, clusters, prior,
                                 remappings = c(PVALB = 'Pvalb', `PV Neurons` = 'Pvalb',
                                                Pyr = 'L[1-6]', `PYC Neurons` = 'L[1-6]',
                                                SST = 'Sst', `SST Neurons` = 'Sst',
                                                VIP = 'Vip', `VIP Neurons` = 'Vip'),
                                 shift = F) {
  do.call(rbind, lapply(as.character(unique(clusters)), function(cluster) {
    p <- remappings[cluster]
    if(is.null(p) || is.na(p))
      p <- cluster
    
    scaled <- t(t(cum.cor[clusters == cluster, ]) / prior[p, colnames(cum.cor)])
    
    scaled + mean(1 - scaled[, p]) * shift
  }))
}

#' Generate an empirical prior from Tasic 2016 data
#'
#' @param ref.collapsed A collapsed reference dataset from @method{simplify_reference} 
#' @param clusters A character vector of regex's that will be grepl'd to the "primary_type_label" of Tasic 2016
#' @param markers.perc The percentage of markers to use in each iteration, relative to the smallest cluster
#' @param n.iter The number of bootstrap iterations
#' @param n.cores The number of CPU cores
#'
#' @return The Tasic 2016 data for the specified clusters mapped onto the provided reference
generate_prior <- function(ref.collapsed,
                           clusters = colnames(ref.collapsed$collapsed), markers.perc = 0.6,
                           n.iter = 1e3, n.cores = 4) {
  require(Matrix)
  require(dplyr)
  require(tasic2016data)
  require(matrixStats)
  
  lapply(clusters, function(type) {
    which(grepl(type, tasic_2016_anno$primary_type_label))
  }) %>% setNames(clusters) -> t2016_samples
  
  tasic <- Matrix(tasic_2016_counts[, tasic_2016_anno$sample_name][, unlist(t2016_samples)], sparse = T)
  map_samples(cpm(tasic, T), ref.collapsed, markers.perc, n.iter, n.cores) -> mapping_ref
  
  mapping.ref <- t((t(mapping_ref$cum.correlation) / mapping_ref$iterations))
  mapping.ref.meta <- data.frame(cluster = unlist(lapply(names(t2016_samples), function(type) rep(type, length(t2016_samples[[type]])))))
  
  lapply(unique(mapping.ref.meta$cluster), function(cluster) {
    colMedians(mapping.ref[mapping.ref.meta$cluster == cluster, ])
  }) %>% setNames(unique(mapping.ref.meta$cluster))
}

#' Generate an empirical prior from Tasic 2016 data. This version of generate_prior uses a within-dataset correlation
#' instead of correlation to Tasic 2016.
#'
#' @param ref The reference data
#' @param ref.collapsed The collapsed reference dataset from @method{simplify_reference} 
#' @param clusters A character vector of cell types to consider
#' @param markers.perc The percentage of markers to use in each iteration, relative to the smallest cluster
#' @param max.cells.per.ident The maximum number of cells to try to map per cluster
#' @param n.iter The number of bootstrap iterations
#' @param n.cores The number of CPU cores
#'
#' @return The Tasic 2016 data for the specified clusters mapped onto the provided reference
generate_prior <- function(ref, ref.collapsed, clusters = levels(Idents(ref)), markers.perc = 0.6,
                           max.cells.per.ident = 200, random.seed = 1,
                           n.iter = 1e3, n.cores = 4) {
  require(dplyr)
  require(matrixStats)
  
  # Subset a smaller number of cells to use
  ref <- SubsetData(ref, ident.use = clusters, max.cells.per.ident = max.cells.per.ident, random.seed = random.seed)
  
  map_samples(ref@assays$RNA@data, ref.collapsed, clusters, markers.perc, n.iter, n.cores) -> mapping_ref
  
  mapping.ref <- mapping_ref$cum.correlation / mapping_ref$iterations
  
  mapping.collapsed <- do.call(rbind, lapply(clusters, function(cluster) {
    colMedians(mapping.ref[ref$cluster == cluster, ], na.rm = T)
  })) %>% `rownames<-`(clusters) %>% `colnames<-`(clusters)
  
  list(full = mapping.ref, collapsed = mapping.collapsed)
}

#' Generate contamination covariates. This method will run the entire pipeline and return results in a named list.
#'
#' @param test.data A dgCMatrix of counts
#' @param test.meta A data frame of metadata. Must have a "cluster" column
#' @param ref.collapsed A collapsed reference dataset from @method{simplify_reference}
#' @param prior The prior from @method{(generate_prior)}, or NULL to calculate it automatically
#' @param markers.perc The percentage of markers to use in each iteration, relative to the smallest cluster
#' @param n.iter The number of bootstrap iterations
#' @param n.cores The number of CPU cores
#' @param remappings A named character of remappings to ensure that cluster names in @param{clusters} and @param{prior} are normalized
#' @param seed A random seed, or NULL for a random seed
#'
#' @return A list with the (un-normalized) mapping, prior, and mapping normalized to the prior
generate_covariates <- function(test.data, test.meta, ref.collapsed, prior = NULL,
                                markers.perc = 0.6, n.iter = 1e4, n.cores = 4,
                                remappings = c(PVALB = 'Pvalb', `PV Neurons` = 'Pvalb',
                                               Pyr = 'L[1-6]', `PYC Neurons` = 'L[1-6]',
                                               SST = 'Sst', `SST Neurons` = 'Sst',
                                               VIP = 'Vip', `VIP Neurons` = 'Vip'),
                                seed = NULL) {
  require(Seurat)
  require(tasic2016data)
  
  if(!is.null(seed))
    set.seed(seed)
  
  test <- as_Seurat(test.data, test.meta)
  
  if(is.null(prior))
    prior <- generate_prior(ref.collapsed, markers.perc, n.iter, n.cores)
  
  mapping <- map_samples(test@assays$RNA@data, ref.collapsed, markers.perc, n.iter, n.cores)
  
  list(mapping = mapping,
       prior = prior,
       normalized = normalize_covariates(t(t(mapping$cum.correlation) / mapping$iterations),
                                         clusters, prior))
}

#' Analyze synthetic data
#'
#' @param celltypes A list of cell types to generate data for (defaults to the column names of @param{ref.collapsed}$collapsed)
#' @param ref.collapsed A collapsed reference dataset from @method{simplify_reference}
#' @param prior The prior from @method{(generate_prior)}
#' @param markers.perc The percentage of markers to use in each iteration, relative to the smallest cluster
#' @param n.iter The number of bootstrap iterations
#' @param n.cores The number of CPU cores
#' @param remappings A named character of remappings to ensure that cluster names in @param{celltypes} and @param{prior} are normalized
#' @param n.genes How many genes to simulate
#' @param n.samples.per.cond How many samples to simulate per condition
#' @param n.DEGs How many DEGs to simulate
#' @param effect.size The mean effect size for DEGs
#' @param contamination A square matrix for contamination with the same length as @code{length(celltypes)}. Each row represents the amount of contamination by that row onto each column.
#' @param identity.min
#' @param identity.max Parameters to control how similar the simulated clusters are to their real counterparts
#' @param within.sample.min
#' @param within.sample.max Parameters to control how much the contamination differs per expressed gene
#' @param between.sample.min
#' @param between.sample.max Parameters to control how much the contamination differs between samples
#' 
#' @return The synthetic data, both contaminated and not, with normalized and unnormalized contamination coefficients
run_synthetic <- function(synthetic = NULL, celltypes = colnames(ref.collapsed$collapsed),
                          ref.collapsed, prior = NULL,
                          markers.perc = 0.6, n.iter = 1e4, n.cores = 4,
                          remappings = c(PVALB = 'Pvalb', `PV Neurons` = 'Pvalb',
                                         Pyr = 'L[1-6]', `PYC Neurons` = 'L[1-6]',
                                         SST = 'Sst', `SST Neurons` = 'Sst',
                                         VIP = 'Vip', `VIP Neurons` = 'Vip'),
                          n.genes = 3e4, n.samples.per.cond = 10, n.DEGs = n.genes / 3e2, effect.size = 4,
                          contamination = matrix(0, length(celltypes), length(celltypes)), confound = 0.5,
                          identity.min = 0.9, identity.max = 1.1,
                          within.sample.min = 0.4, within.sample.max = 1.6,
                          between.sample.min = 0.7, between.sample.max = 1.3) {
  require(compcodeR)
  require(Seurat)
  require(dplyr)
  require(Matrix)
  require(dqrng)
  
  # Redefine this because the implementation as provided was too slow to use
  generateSyntheticData <- function (dataset, n.vars, samples.per.cond, n.diffexp, repl.id = 1, 
            seqdepth = 1e+07, minfact = 0.7, maxfact = 1.4, relmeans = "auto", 
            dispersions = "auto", fraction.upregulated = 1, between.group.diffdisp = FALSE, 
            filter.threshold.total = 1, filter.threshold.mediancpm = 0, 
            fraction.non.overdispersed = 0, random.outlier.high.prob = 0, 
            random.outlier.low.prob = 0, single.outlier.high.prob = 0, 
            single.outlier.low.prob = 0, effect.size = 1.5, output.file = NULL) 
  {
    if (!is.null(output.file)) {
      if (!(substr(output.file, nchar(output.file) - 3, nchar(output.file)) == 
            ".rds")) {
        stop("output.file must be an .rds file.")
      }
    }
    uID <- paste(sample(c(0:9, letters, LETTERS), 10, replace = TRUE), 
                 collapse = "")
    condition <- rep(c(1, 2), each = samples.per.cond)
    S1 <- which(condition == 1)
    S2 <- which(condition == 2)
    if (length(effect.size) == 1) {
      n.upregulated <- floor(fraction.upregulated * n.diffexp)
      if (fraction.upregulated != 0 & n.diffexp != 0) {
        genes.upreg <- 1:n.upregulated
      }
      else {
        genes.upreg <- NULL
      }
      if (fraction.upregulated != 1 & n.diffexp != 0) {
        genes.downreg <- (n.upregulated + 1):n.diffexp
      }
      else {
        genes.downreg <- NULL
      }
      genes.nonreg <- setdiff(1:n.vars, union(genes.upreg, 
                                              genes.downreg))
    }
    else {
      if (length(effect.size) != n.vars) {
        stop("The length of the effect.size vector must be the same as the number of simulated genes.")
      }
      else {
        genes.upreg <- which(effect.size > 1)
        genes.downreg <- which(effect.size < 1)
        genes.nonreg <- which(effect.size == 1)
        n.upregulated <- length(genes.upreg)
        n.diffexp <- length(genes.upreg) + length(genes.downreg)
        fraction.upregulated <- n.upregulated/n.diffexp
      }
    }
    differential.expression <- rep(0, n.vars)
    differential.expression[genes.upreg] <- 1
    differential.expression[genes.downreg] <- 1
    upregulation <- rep(0, n.vars)
    upregulation[genes.upreg] <- 1
    downregulation <- rep(0, n.vars)
    downregulation[genes.downreg] <- 1
    if (is.character(relmeans) | is.character(dispersions)) {
      mu.phi.estimates <- system.file("extdata", "Pickrell.Cheung.Mu.Phi.Estimates.rds", 
                                      package = "compcodeR")
      mu.phi.estimates <- readRDS(mu.phi.estimates)
      mu.estimates <- mu.phi.estimates$pickrell.cheung.mu
      phi.estimates <- mu.phi.estimates$pickrell.cheung.phi
      to.include <- sample(1:length(mu.estimates), n.vars, 
                           replace = ifelse(n.vars > length(mu.estimates), TRUE, 
                                            FALSE))
      truedispersions.S1 <- phi.estimates[to.include]
      truemeans.S1 <- mu.estimates[to.include]
    }
    if (!is.character(relmeans)) {
      if (length(relmeans) != n.vars) 
        stop("The length of the relmeans vector must be the same as the number of simulated genes.")
      truemeans.S1 <- c(relmeans)
    }
    if (!is.character(dispersions)) {
      if (nrow(cbind(dispersions)) != n.vars) 
        stop("The number of provided dispersions must be the same as the number of simulated genes.")
      truedispersions.S1 <- cbind(dispersions)[, 1]
      if (ncol(cbind(dispersions)) > 1) {
        truedispersions.S2 <- cbind(dispersions)[, 2]
      }
      else {
        truedispersions.S2 <- truedispersions.S1
      }
    }
    nfacts <- runif(2 * samples.per.cond, min = minfact, max = maxfact)
    seq.depths <- nfacts * seqdepth
    overdispersed <- rep(1, n.vars)
    if (fraction.non.overdispersed > 0) {
      overdispersed[genes.upreg[1:round(fraction.non.overdispersed * 
                                          length(genes.upreg))]] <- 0
      overdispersed[genes.downreg[1:round(fraction.non.overdispersed * 
                                            length(genes.downreg))]] <- 0
      overdispersed[genes.nonreg[1:round(fraction.non.overdispersed * 
                                           length(genes.nonreg))]] <- 0
    }
    prob.S1 <- truemeans.S1
    prob.S2 <- rep(0, length(prob.S1))
    if (length(effect.size) == 1) {
      for (i in 1:n.vars) {
        if (i %in% genes.upreg) {
          prob.S2[i] <- (effect.size + rexp(1, rate = 1)) * 
            prob.S1[i]
        }
        else {
          if (i %in% genes.downreg) {
            prob.S2[i] <- 1/(effect.size + rexp(1, rate = 1)) * 
              prob.S1[i]
          }
          else {
            prob.S2[i] <- prob.S1[i]
          }
        }
      }
    }
    else {
      prob.S2 <- c(effect.size) * prob.S1
    }
    true.log2foldchange <- log2(prob.S2/prob.S1)
    sum.S1 <- sum(prob.S1)
    sum.S2 <- sum(prob.S2)
    if (is.character(dispersions)) {
      truedispersions.S2 <- truedispersions.S1
      if (between.group.diffdisp == TRUE) {
        for (i in 1:length(truedispersions.S2)) {
          sample.base <- phi.estimates[abs(log10(mu.estimates) - 
                                             log10(prob.S2[i])) < 0.05]
          if (length(sample.base) < 50) {
            sample.base <- phi.estimates[order(abs(log10(mu.estimates) - 
                                                     log10(prob.S2[i])))][1:500]
          }
          truedispersions.S2[i] <- sample(sample.base, 
                                          1)
        }
      }
    }
    truedispersions.S1 <- truedispersions.S1 * overdispersed
    truedispersions.S2 <- truedispersions.S2 * overdispersed
    do.call(cbind, lapply(1:(length(S1) + length(S2)), function(j) {
      if(j %in% S1) {
        if (overdispersed[i] == 1) {
          rnbinom(n.vars, mu = prob.S1/sum.S1 * 
                               seq.depths[j], size = 1/truedispersions.S1)
        }
        else {
          rpois(n.vars, lambda = prob.S1/sum.S1 * 
                             seq.depths[j])
        }
      } else {
        if (overdispersed[i] == 1) {
          rnbinom(n.vars, mu = prob.S2/sum.S2 * 
                               seq.depths[j], size = 1/truedispersions.S2)
        }
        else {
          rpois(n.vars, lambda = prob.S2/sum.S2 * 
                             seq.depths[j])
        }
      }
    })) -> Z
    
    random.outliers <- matrix(0, nrow(Z), ncol(Z))
    random.outliers.factor <- matrix(1, nrow(Z), ncol(Z))
    if (random.outlier.high.prob != 0 | random.outlier.low.prob != 
        0) {
      for (i in 1:nrow(Z)) {
        for (j in 1:ncol(Z)) {
          tmp <- runif(1)
          if (tmp < random.outlier.high.prob) {
            random.outliers[i, j] <- 1
            random.outliers.factor[i, j] <- runif(1, min = 5, 
                                                  max = 10)
          }
          else if (tmp < random.outlier.low.prob + random.outlier.high.prob) {
            random.outliers[i, j] <- (-1)
            random.outliers.factor[i, j] <- 1/runif(1, 
                                                    min = 5, max = 10)
          }
        }
      }
      Z <- round(random.outliers.factor * Z)
    }
    has.single.outlier <- rep(0, n.vars)
    single.outliers <- matrix(0, nrow(Z), ncol(Z))
    single.outliers.factor <- matrix(1, nrow(Z), ncol(Z))
    if (single.outlier.high.prob != 0 | single.outlier.low.prob != 
        0) {
      has.single.outlier[genes.upreg[1:floor((single.outlier.high.prob + 
                                                single.outlier.low.prob) * length(genes.upreg))]] <- 1
      has.single.outlier[genes.downreg[1:floor((single.outlier.high.prob + 
                                                  single.outlier.low.prob) * length(genes.downreg))]] <- 1
      has.single.outlier[genes.nonreg[1:floor((single.outlier.high.prob + 
                                                 single.outlier.low.prob) * length(genes.nonreg))]] <- 1
      for (i in 1:nrow(Z)) {
        if (has.single.outlier[i] == 1) {
          the.sample <- sample(1:(ncol(Z)), 1)
          if (runif(1) < (single.outlier.high.prob/(single.outlier.high.prob + 
                                                    single.outlier.low.prob))) {
            single.outliers[i, the.sample] <- 1
            single.outliers.factor[i, the.sample] <- runif(1, 
                                                           min = 5, max = 10)
          }
          else {
            single.outliers[i, the.sample] <- (-1)
            single.outliers.factor[i, the.sample] <- 1/runif(1, 
                                                             min = 5, max = 10)
          }
        }
      }
      Z <- round(single.outliers.factor * Z)
    }
    rownames(Z) <- 1:n.vars
    n.random.outliers.up.S1 <- apply(random.outliers[, S1] > 
                                       0, 1, sum)
    n.random.outliers.up.S2 <- apply(random.outliers[, S2] > 
                                       0, 1, sum)
    n.random.outliers.down.S1 <- apply(random.outliers[, S1] < 
                                         0, 1, sum)
    n.random.outliers.down.S2 <- apply(random.outliers[, S2] < 
                                         0, 1, sum)
    n.single.outliers.up.S1 <- apply(single.outliers[, S1] > 
                                       0, 1, sum)
    n.single.outliers.up.S2 <- apply(single.outliers[, S2] > 
                                       0, 1, sum)
    n.single.outliers.down.S1 <- apply(single.outliers[, S1] < 
                                         0, 1, sum)
    n.single.outliers.down.S2 <- apply(single.outliers[, S2] < 
                                         0, 1, sum)
    nf <- calcNormFactors(Z)
    norm.factors <- nf * colSums(Z)
    common.libsize <- exp(mean(log(colSums(Z))))
    pseudocounts <- sweep(Z + 0.5, 2, norm.factors, "/") * common.libsize
    log2.pseudocounts <- log2(pseudocounts)
    M.value <- apply(log2.pseudocounts[, S2], 1, mean) - apply(log2.pseudocounts[, 
                                                                                 S1], 1, mean)
    A.value <- 0.5 * (apply(log2.pseudocounts[, S2], 1, mean) + 
                        apply(log2.pseudocounts[, S1], 1, mean))
    variable.annotations <- data.frame(truedispersions.S1 = truedispersions.S1, 
                                       truedispersions.S2 = truedispersions.S2, truemeans.S1 = prob.S1, 
                                       truemeans.S2 = prob.S2, n.random.outliers.up.S1 = n.random.outliers.up.S1, 
                                       n.random.outliers.up.S2 = n.random.outliers.up.S2, n.random.outliers.down.S1 = n.random.outliers.down.S1, 
                                       n.random.outliers.down.S2 = n.random.outliers.down.S2, 
                                       n.single.outliers.up.S1 = n.single.outliers.up.S1, n.single.outliers.up.S2 = n.single.outliers.up.S2, 
                                       n.single.outliers.down.S1 = n.single.outliers.down.S1, 
                                       n.single.outliers.down.S2 = n.single.outliers.down.S2, 
                                       M.value = M.value, A.value = A.value, truelog2foldchanges = true.log2foldchange, 
                                       upregulation = upregulation, downregulation = downregulation, 
                                       differential.expression = differential.expression)
    rownames(variable.annotations) <- rownames(Z)
    sample.annotations <- data.frame(condition = condition, depth.factor = nfacts)
    info.parameters <- list(n.diffexp = n.diffexp, fraction.upregulated = fraction.upregulated, 
                            between.group.diffdisp = between.group.diffdisp, filter.threshold.total = filter.threshold.total, 
                            filter.threshold.mediancpm = filter.threshold.mediancpm, 
                            fraction.non.overdispersed = fraction.non.overdispersed, 
                            random.outlier.high.prob = random.outlier.high.prob, 
                            random.outlier.low.prob = random.outlier.low.prob, single.outlier.high.prob = single.outlier.high.prob, 
                            single.outlier.low.prob = single.outlier.low.prob, effect.size = effect.size, 
                            samples.per.cond = samples.per.cond, repl.id = repl.id, 
                            dataset = dataset, uID = uID, seqdepth = seqdepth, minfact = minfact, 
                            maxfact = maxfact)
    s <- apply(Z, 1, sum)
    keep.T <- which(s >= filter.threshold.total)
    Z.T <- Z[keep.T, ]
    variable.annotations.T <- variable.annotations[keep.T, ]
    filtering <- paste("total count >=", filter.threshold.total)
    cpm <- sweep(Z.T, 2, apply(Z.T, 2, sum), "/") * 1e+06
    m <- apply(cpm, 1, median)
    keep.C <- which(m >= filter.threshold.mediancpm)
    Z.TC <- Z.T[keep.C, ]
    variable.annotations.TC <- variable.annotations.T[keep.C, 
    ]
    filtering <- paste(filtering, "; ", paste("median cpm >=", 
                                              filter.threshold.mediancpm))
    rownames(Z.TC) <- paste("g", 1:nrow(Z.TC), sep = "")
    colnames(Z.TC) <- paste("sample", 1:ncol(Z.TC), sep = "")
    rownames(sample.annotations) <- colnames(Z.TC)
    rownames(variable.annotations.TC) <- rownames(Z.TC)
    data.object <- compData(count.matrix = Z.TC, variable.annotations = variable.annotations.TC, 
                            sample.annotations = sample.annotations, filtering = filtering, 
                            info.parameters = info.parameters)
    if (!is.null(output.file)) {
      saveRDS(data.object, file = output.file)
    }
    return(invisible(data.object))
  }
  
  if(is.null(synthetic)) {
    message('Generating synthetic data...')
    lapply(celltypes, function(celltype) {
      message(paste('Cell type:', celltype))
      synth <- generateSyntheticData('synthetic', n.genes, n.samples.per.cond, n.DEGs,
                                     fraction.upregulated = 0.5,
                                     filter.threshold.total = 0,
                                     effect.size = effect.size)
      subsample <- dqsample(1:n.genes, n.genes)
      
      synth@count.matrix <- synth@count.matrix[subsample, ] %>% `rownames<-`(paste0('g', 1:n.genes))
      synth@variable.annotations <- synth@variable.annotations[subsample, ] %>% `rownames<-`(paste0('g', 1:n.genes))
      
      colnames(synth@count.matrix) <- paste0(celltype, paste0('-', 1:(n.samples.per.cond * 2)))
      synth@sample.annotations$condition <- as.logical(synth@sample.annotations$condition - 1)
      rownames(synth@sample.annotations) <- colnames(synth@count.matrix)
      
      list(count = synth@count.matrix, meta = synth@sample.annotations, vars = synth@variable.annotations)
    }) %>% setNames(celltypes) -> synthetics
    
    synthetic <- CreateSeuratObject(Matrix(do.call(cbind, lapply(synthetics, '[[', 1)), sparse = T))
    synthetic@assays$RNA@data <- cpm(synthetic@assays$RNA@counts, log = T)
    synthetic@meta.data <- cbind(synthetic@meta.data, do.call(rbind, lapply(synthetics, '[[', 2)))
    synthetic@meta.data$cluster <- rep(names(synthetics), each = n.samples.per.cond * 2)
    synthetic@meta.data$sample <- rownames(synthetic@meta.data)
    Idents(synthetic) <- synthetic@meta.data$cluster
    
    synthetic@assays$RNA@data <- rbind(synthetic@assays$RNA@data,
                                       matrix(nrow = nrow(ref.collapsed$markers),
                                              ncol = ncol(synthetic@assays$RNA@data)) %>%
                                         `rownames<-`(ref.collapsed$markers$gene))
    
    message('Adding cellular identities...')
    for(celltype in names(synthetics)) {
      synthetic@assays$RNA@data[ref.collapsed$markers$gene, synthetic@meta.data$cluster == celltype] <-
        ref.collapsed$collapsed[ref.collapsed$markers$gene, celltype] *
        matrix(runif(length(ref.collapsed$markers$gene) * sum(synthetic@meta.data$cluster == celltype),
                     identity.min, identity.max),
               nrow = length(ref.collapsed$markers$gene), byrow = T)
    }
  }
  
  contaminate <- function(data, contamination, confound,
                          within.sample.min, within.sample.max,
                          between.sample.min, between.sample.max) {
    contam <- function(data, pure, contamination, i, celltypes, exprA, exprB, confounding) {
      for(j in 1:ncol(contamination)) {
        if(i != j && contamination[i, j] > 0) {
          message(paste('Contaminating', celltypes[[j]]))
          
          pureExprA <- pure[, data$cluster == celltypes[j] & data$condition]
          pureExprA[is.na(pureExprA)] <- 0
          pureExprB <- pure[, data$cluster == celltypes[j] & !data$condition]
          pureExprB[is.na(pureExprB)] <- 0

          # The new expression level is:
          # Expr_old + contamination * expr_diff(noised within)(columnwise) * noise_between(rowwise) * confound
          data@assays$RNA@data[, data$cluster == celltypes[j] & data$condition] <-
            (data@assays$RNA@data[, data$cluster == celltypes[j] & data$condition] +
               # Old expression
               contamination[i, j] *
               t(t((exprA - pureExprA) * dqrunif(nrow(pure), within.sample.min, within.sample.max)) * # Randomness within columns
                   dqrunif(sum(data$cluster == celltypes[j] & data$condition), between.sample.min, between.sample.max) * # Randomness between columns
                   confounding[data$cluster == celltypes[j] & data$condition])) %>% # Confounding
            pmax(0)
          
          data@assays$RNA@data[, data$cluster == celltypes[j] & !data$condition] <-
            (data@assays$RNA@data[, data$cluster == celltypes[j] & !data$condition] +
               # Old expression
               contamination[i, j] *
               t(t((exprB - pureExprB) * dqrunif(nrow(pure), within.sample.min, within.sample.max)) * # Randomness within columns
                   dqrunif(sum(data$cluster == celltypes[j] & !data$condition), between.sample.min, between.sample.max) * # Randomness between columns
                   confounding[data$cluster == celltypes[j] & !data$condition])) %>% # Confounding
            pmax(0)
        }
      }
      
      data
    }
    
    tmp <- data@assays$RNA@data
    
    if(confound == 0.5)
      confounding <- rep(T, ncol(tmp))
    else
      confounding <- (data$condition == sample(c(T, F), ncol(tmp), T, c(confound, 1 - confound))) + runif(ncol(tmp), -0.1, 0.1)
    
    for(i in 1:nrow(contamination)) {
      message(paste('Contaminating with', celltypes[[i]]))
      data <- contam(data, tmp, contamination, i, celltypes,
                     Matrix::rowMeans(data@assays$RNA@data[, data$cluster == celltypes[i] & data$condition]),
                     Matrix::rowMeans(data@assays$RNA@data[, data$cluster == celltypes[i] & !data$condition]),
                     confounding)
    }
    
    data
  }
  
  synthetic.contaminated <- contaminate(synthetic, contamination, confound,
                                        within.sample.min, within.sample.max, between.sample.min, between.sample.max)
  
  synthetic.contaminated@assays$RNA@counts <- t(t(2^(synthetic.contaminated@assays$RNA@data) - 1) *
                                                  Matrix::colSums(synthetic.contaminated@assays$RNA@counts) / 1e6) %>%
    round %>% as('dgCMatrix') %>% pmax(0)
  
  if(is.null(prior))
    prior <- generate_prior(ref.collapsed, markers.perc, n.iter, n.cores)
  
  mapping.prior <- map_samples(synthetic@assays$RNA@data, ref.collapsed, colnames(ref.collapsed$collapsed), markers.perc, n.iter, n.cores)
  mapping <- map_samples(synthetic.contaminated@assays$RNA@data, ref.collapsed, colnames(ref.collapsed$collapsed), markers.perc, n.iter, n.cores)
  
  ret <- list(pure = synthetic,
              contamination.parameters = list(confound = confound,
                                              within.sample.min = within.sample.min,
                                              within.sample.max = within.sample.max,
                                              between.sample.min = between.sample.max,
                                              contamination = contamination),
              contaminated = synthetic.contaminated,
              mapping.contaminated = mapping,
              mapping.pure = mapping.prior,
              normalized = normalize_covariates(mapping$cum.correlation / n.iter,
                                                synthetic.contaminated$cluster, prior, remappings))
  
  if(exists('synthetics'))
    ret$meta <- lapply(synthetics, '[[', 3) %>% setNames(names(synthetics))
  
  ret
}