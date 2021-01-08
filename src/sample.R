loadAging <- function(exons = T, introns = T) {
  require(Matrix)
  
  newton.aging <- 0
  if(exons)
    newton.aging <- newton.aging + readRDS(paste0('data/Aging/counts_exon.rds'))
  if(introns)
    newton.aging <- newton.aging + readRDS(paste0('data/Aging/counts_intron.rds'))
  
  colnames(newton.aging) <- gsub('(_|S[1-9].*)', '', colnames(newton.aging))
  
  # Clean up some problematic samples
  newton.aging <- newton.aging[, -2]
  newton.aging <- newton.aging[, -c(36, 40)]
  
  newton.aging.meta <- readRDS('data/Aging/GSE119183.rds')
  
  qc <- read.csv('data/Aging/qc.tsv', sep = '\t')
  newton.aging.meta$mitochondrial <- qc$mitochondrial[-c(2, 37, 41)] / Matrix::colSums(newton.aging)
  
  rownames(newton.aging.meta) <- newton.aging.meta$sample
  list(data = newton.aging, meta = newton.aging.meta)
}

loadStress <- function(exons = T, introns = T) {
  require(Matrix)
  
  newton.stress <- 0
  if(exons)
    newton.stress <- newton.stress + readRDS('data/UCMS/counts_exon.rds')
  if(introns)
    newton.stress <- newton.stress + readRDS('data/UCMS/counts_intron.rds')
  
  colnames(newton.stress) <- gsub('(_|S[1-9].*)', '', colnames(newton.stress))
  
  newton.stress.meta <- data.frame(sample = colnames(newton.stress),
                                   cluster = gsub('-.*', '', colnames(newton.stress)),
                                   stress = grepl('-U6[0-9]{2}', colnames(newton.stress)))
  
  qc <- read.csv('data/UCMS/qc.tsv', sep = '\t')
  newton.stress.meta$mitochondrial <- qc$mitochondrial / Matrix::colSums(newton.stress)
  
  rownames(newton.stress.meta) <- newton.stress.meta$sample
  list(data = newton.stress, meta = newton.stress.meta)
}

loadReference <- function(exons = T, introns = T,
                          clusters = c('L[1-6]', 'Astro', 'Endo', 'Pvalb', 'Vip', 'Sst', 'Oligo', 'DG', 'Macrophage'),
                          regions = '.*') {
  require(scrattch.io)
  require(dplyr)
  
  tasic <- 0
  if(exons)
    tasic <- tasic + read_tome_dgCMatrix('data/AllenInstitute/Mouse/transcrip.tome', 'data/t_exon')
  if(introns)
    tasic <- tasic + read_tome_dgCMatrix('data/AllenInstitute/Mouse/transcrip.tome', 'data/t_intron')
  
  tasic.meta <- read.csv('data/AllenInstitute/Mouse/sample_annotations.csv') %>%
    .[(grepl(paste0('(', paste0(clusters, collapse = '|'), ')'), .[, 'subclass_label']) |
         (.[, 'class_label'] == 'Glutamatergic' &
            grepl(paste0('(', paste0(clusters, collapse = '|'), ')'), 
                  sapply(as.character(.[, 'cluster_label']), function(label) strsplit(label, ' ', T)[[1]][2])))) &
        (grepl(paste0('(', paste0(regions, collapse = '|'), ')'), .[, 'region_label']) |
           grepl('(Astro|Endo|Oligo)', .[, 'subclass_label'])), ]
  
  rownames(tasic) <- read_tome_gene_names('data/AllenInstitute/Mouse/transcrip.tome')
  colnames(tasic) <- read_tome_sample_names('data/AllenInstitute/Mouse/transcrip.tome')
  
  tasic <- tasic[, tasic.meta$sample_name]
  
  tasic.meta$cluster <- ''
  for(cluster in clusters) {
    tasic.meta$cluster[grepl(cluster, tasic.meta$subclass_label)] <- cluster
  }
  
  tasic.meta$cluster <- factor(tasic.meta$cluster)
  list(data = tasic, meta = tasic.meta)
}

deg <- function(data, meta, model) {
  register(MulticoreParam(8))
  DESeq(DESeqDataSetFromMatrix(countData = data, colData = meta, design = model) %>%
          .[rowSums(counts(.[]) > 0) >= 10, ], parallel = T)
}

do.DEG <- function(data, model, mapping = NULL) {
  if(!is.null(mapping)) {
    data@meta.data <- merge(data@meta.data,
                            data.frame(apply(mapping, 2, scale),
                                       sample = rownames(mapping)), by = 'sample') %>%
      `rownames<-`(.[, 'sample'])
  }
  
  lapply(Idents(data) %>% unique, function(cluster) {
    dataSubset <- SubsetData(data, ident.use = cluster)
    deg(dataSubset@assays$RNA@counts,
        dataSubset@meta.data,
        model)
  }) -> DEGs
  
  names(DEGs) <- Idents(data) %>% unique
  
  DEGs
}

trace(glmer, edit = T)

TASIC <- loadReference(clusters = c('L[1-6]', 'Astro', 'Endo', 'Pvalb', 'Vip', 'Sst', 'Oligo', 'Macrophage'), regions = c('ACA', 'CA'))
TASIC <- as_Seurat(TASIC$data, TASIC$meta)
REFERENCE <- simplify_reference(TASIC)
PRIOR <- generate_prior(TASIC, REFERENCE, clusters = c('Astro', 'L[1-6]', 'Endo', 'Pvalb', 'Vip', 'Sst', 'Oligo'),
                        max.cells.per.ident = 900)

STRESS <- loadStress()
AGING <- loadAging()

STRESS <- as_Seurat(STRESS$data, STRESS$meta)
AGING <- as_Seurat(AGING$data, AGING$meta)

STRESS.MAPPING <- map_samples(STRESS@assays$RNA@data, REFERENCE)
AGING.MAPPING <- map_samples(AGING@assays$RNA@data, REFERENCE)

t2016_anno <- tasic_2016_anno
t2016_anno$sample <- t2016_anno$sample_name
t2016_anno$cluster <- 'Other'
for(cluster in c('L[1-6]', 'Astro', 'Endo', 'Micro', 'Oligo', 'Pvalb', 'Sst', 'Vip')) {
  t2016_anno$cluster[grepl(cluster, t2016_anno$primary_type_label)] <- cluster
}

t2016 <- as_Seurat(Matrix(tasic_2016_counts[, t2016_anno$sample]), t2016_anno)
t2016.MAPPING <- map_samples(t2016@assays$RNA@data, REFERENCE)

rm(t2016_anno, cluster)

####

TASIC@assays$RNA@data[c('Slc17a7', 'Sst', 'Pvalb', 'Vip'), ] %>% as.matrix %>% melt %>%
  merge(TASIC@meta.data[, c('sample_name', 'cluster')], by.x = 'Var2', by.y = 'sample_name', all.x = T) %>%
  .[.[, 'cluster'] %in% c('L[1-6]', 'Pvalb', 'Sst', 'Vip'), ] %>% cbind(experiment = 'Single Cell') %>%
  mutate(cluster = factor(cluster, labels = c('PYC Neurons',
                                              'PV Neurons',
                                              'SST Neurons',
                                              'VIP Neurons'))) %>%
  mutate(cluster = factor(cluster, levels = c('PYC Neurons', 'SST Neurons', 'PV Neurons', 'VIP Neurons'))) %>%
  filter(Var1 %in% c('Slc17a7', 'Sst')) %>%
  ggplot(aes(cluster, value, fill = Var1)) + geom_point() + geom_boxplot()+ xlab(element_blank()) +
  scale_y_continuous(name = 'Expression Level', expand = c(0, 0)) +
  facet_rep_grid(vars(Var1), vars(experiment), scales = 'free_y', switch ='y') +
  theme_set(theme_classic(base_size = 20, base_line_size = 0.25, base_rect_size = 0.25)) +
  theme(strip.background = element_blank(), strip.text.x = element_text(hjust = 0), strip.text.y = element_text(face = 'italic', hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.text = element_text(face = 'italic'), legend.position = 'none',
        strip.placement = 'outside') +
  scale_fill_manual(values = setNames(c('#ee2210',
                                        '#fa6b0a',
                                        '#14983d',
                                        '#149bed'), c('Pvalb', 'Slc17a7', 'Sst', 'Vip')), name = 'Marker Gene') +
  labs(tag = '')

(cpm(AGING@assays$RNA@counts[c('Slc17a7', 'Sst', 'Pvalb', 'Vip'), ])/1e5) %>% as.matrix %>% melt %>%
  merge(AGING@meta.data[, c('sample', 'cluster')], by.x = 'Var2', by.y = 'sample', all.x = T) %>%
  cbind(experiment = 'Aging') %>%
  rbind(
    (cpm(STRESS@assays$RNA@counts[c('Slc17a7', 'Sst', 'Pvalb', 'Vip'), ])/1e5) %>% as.matrix %>% melt %>%
      merge(STRESS@meta.data[, c('sample', 'cluster')], by.x = 'Var2', by.y = 'sample', all.x = T) %>%
      cbind(experiment = 'UCMS') %>% mutate(cluster = factor(cluster, labels = c('PV Neurons',
                                                                                 'PYC Neurons',
                                                                                 'SST Neurons',
                                                                                 'VIP Neurons')))
  ) %>%
  mutate(cluster = factor(cluster, levels = c('PYC Neurons', 'SST Neurons', 'PV Neurons', 'VIP Neurons'))) %>%
  filter(Var1 %in% c('Slc17a7', 'Sst')) %>%
  ggplot(aes(cluster, value, fill = Var1)) + geom_boxplot() + xlab(element_blank()) +
  scale_y_continuous(name = '', expand = c(0, 0)) +
  facet_rep_wrap(~Var1+experiment, ncol = 2, scales = 'free_y', labeller = function(x) { list(c('Aging', 'Stress', '' ,'')) }) +
  theme_set(theme_classic(base_size = 20, base_line_size = 0.25, base_rect_size = 0.25)) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.text = element_text(face = 'italic'), legend.position = 'none',
        strip.placement = 'outside') +
  scale_fill_manual(values = setNames(c('#ee2210',
                                        '#fa6b0a',
                                        '#14983d',
                                        '#149bed'), c('Pvalb', 'Slc17a7', 'Sst', 'Vip')), name = 'Marker Gene')

normalize_covariates(t2016.MAPPING$cum.correlation/t2016.MAPPING$iterations, t2016$cluster, PRIOR$collapsed) %>% melt %>%
  merge(t2016@meta.data[, c('sample_name', 'cluster')], by.x = 'Var1', by.y = 'sample_name', all.x = T) %>%
  mutate(Var2 = factor(Var2, labels = c('Astrocytes', 'Endothelial Cells', 'PYC Neurons', 'Oligodendrocytes',
                                        'PV Neurons', 'SST Neurons', 'VIP Neurons'))) %>%
  mutate(Var2 = factor(Var2, levels = c('Astrocytes', 'Endothelial Cells', 'Oligodendrocytes',
                                        'PYC Neurons', 'SST Neurons', 'PV Neurons', 'VIP Neurons'))) %>%
  .[.[, 'cluster'] %in% c('L[1-6]', 'Pvalb', 'Sst', 'Vip'), ] %>%
  mutate(cluster = factor(cluster, labels = c('PYC Neurons',
                                              'PV Neurons',
                                              'SST Neurons',
                                              'VIP Neurons'))) %>%
  mutate(cluster = factor(cluster, levels = c('PYC Neurons', 'PV Neurons', 'SST Neurons', 'VIP Neurons'))) %>%
  cbind(experiment = 'Single Cell') %>%
  ggplot(aes(cluster, value, fill = Var2)) + geom_boxplot(outlier.alpha = 0.2) +
  scale_x_discrete(expand = c(0.15, 0)) +
  scale_y_continuous(name = 'Scaled Coefficient', expand = c(0, 0), limits = c(0.8, 1.6)) + xlab(element_blank()) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_fill_manual(values = setNames(c('#fec10c',
                                        '#9a703e',
                                        '#0c5bb0',
                                        '#fa6b0a',
                                        '#ee2210',
                                        '#14983d',
                                        '#149bed'), c('Astrocytes', 'Endothelial Cells',
                                                      'Oligodendrocytes',
                                                      'PYC Neurons',
                                                      'PV Neurons',
                                                      'SST Neurons',
                                                      'VIP Neurons')), name = 'Mapped Cell Type') +
  facet_wrap(~experiment, nrow = 1, scales = 'free') +
  theme_set(theme_classic(base_size = 20, base_line_size = 0.25, base_rect_size = 0.25)) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0),
        legend.position = 'none', axis.text.x = element_blank())

normalize_covariates(AGING.MAPPING$cum.correlation/AGING.MAPPING$iterations, AGING$cluster, PRIOR$collapsed) %>% melt %>%
  merge(AGING@meta.data[, c('sample', 'cluster')], by.x = 'Var1', by.y = 'sample', all.x = T) %>%
  mutate(Var2 = factor(Var2, labels = c('Astrocytes', 'Endothelial Cells', 'PYC Neurons', 'Oligodendrocytes',
                                        'PV Neurons', 'SST Neurons', 'VIP Neurons'))) %>%
  mutate(Var2 = factor(Var2, levels = c('Astrocytes', 'Endothelial Cells', 'Oligodendrocytes',
                                        'PYC Neurons', 'SST Neurons', 'PV Neurons', 'VIP Neurons'))) %>%
  mutate(cluster = factor(cluster, levels = c('PYC Neurons', 'SST Neurons', 'PV Neurons', 'VIP Neurons'))) %>%
  cbind(experiment = 'Aging') %>%
  rbind(
    normalize_covariates(STRESS.MAPPING$cum.correlation/STRESS.MAPPING$iterations, STRESS$cluster, PRIOR$collapsed) %>% melt %>%
      merge(STRESS@meta.data[, c('sample', 'cluster')], by.x = 'Var1', by.y = 'sample', all.x = T) %>%
      mutate(Var2 = factor(Var2, labels = c('Astrocytes', 'Endothelial Cells', 'PYC Neurons', 'Oligodendrocytes',
                                            'PV Neurons', 'SST Neurons', 'VIP Neurons'))) %>%
      mutate(Var2 = factor(Var2, levels = c('Astrocytes', 'Endothelial Cells', 'Oligodendrocytes',
                                            'PYC Neurons', 'SST Neurons', 'PV Neurons', 'VIP Neurons'))) %>%
      mutate(cluster = factor(cluster, labels = c('PV Neurons', 'PYC Neurons', 'SST Neurons', 'VIP Neurons'))) %>%
      mutate(cluster = factor(cluster, levels = c('PYC Neurons', 'SST Neurons', 'PV Neurons', 'VIP Neurons'))) %>%
      cbind(experiment = 'Stress')
  ) %>%
  ggplot(aes(cluster, value, fill = Var2)) + geom_boxplot(outlier.alpha = 0) +
  scale_x_discrete(expand = c(0.15, 0)) +
  scale_y_continuous(name = 'Scaled Coefficient', expand = c(0, 0), limits = c(0.8, 2.4)) + xlab(element_blank()) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_fill_manual(values = setNames(c('#fec10c',
                                        '#9a703e',
                                        '#0c5bb0',
                                        '#fa6b0a',
                                        '#ee2210',
                                        '#14983d',
                                        '#149bed'), c('Astrocytes', 'Endothelial Cells',
                                                      'Oligodendrocytes',
                                                      'PYC Neurons',
                                                      'PV Neurons',
                                                      'SST Neurons',
                                                      'VIP Neurons')), name = 'Mapped Cell Type') +
  facet_wrap(~experiment, nrow = 1, scales = 'free') +
  theme_set(theme_classic(base_size = 20, base_line_size = 0.25, base_rect_size = 0.25)) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0),
        legend.position = 'none', axis.text.x = element_blank())

normalize_covariates(STRESS.MAPPING$cum.correlation/STRESS.MAPPING$iterations, STRESS$cluster, PRIOR$collapsed) %>% melt %>%
  merge(STRESS@meta.data[, c('sample', 'cluster')], by.x = 'Var1', by.y = 'sample') %>%
  group_by(Var2, cluster) %>% summarize(A = mean(value)) %>% ungroup() %>%
  mutate(Var2 = factor(Var2, labels = c('Astrocytes', 'Endothelial Cells', 'PYC Neurons', 'Oligodendrocytes', 'PV Neurons', 'SST Neurons', 'VIP Neurons'))) %>%
  mutate(Var2 = factor(Var2, levels = c('Astrocytes', 'Endothelial Cells', 'Oligodendrocytes', 'PYC Neurons', 'SST Neurons', 'PV Neurons', 'VIP Neurons'))) %>%
  cbind(B = 
          normalize_covariates(AGING.MAPPING$cum.correlation/AGING.MAPPING$iterations, AGING$cluster, PRIOR$collapsed) %>% melt %>%
          merge(AGING@meta.data[, c('sample', 'cluster')], by.x = 'Var1', by.y = 'sample') %>%
          group_by(Var2, cluster) %>% summarize(B = mean(value)) %>% ungroup() %>%
          mutate(Var2 = factor(Var2, labels = c('Astrocytes', 'Endothelial Cells', 'PYC Neurons', 'Oligodendrocytes', 'PV Neurons', 'SST Neurons', 'VIP Neurons'))) %>%
          mutate(Var2 = factor(Var2, levels = c('Astrocytes', 'Endothelial Cells', 'Oligodendrocytes', 'PYC Neurons', 'SST Neurons', 'PV Neurons', 'VIP Neurons'))) %>%
          .[, 'B']) %>% mutate(cluster = factor(cluster, labels = c('PV Neurons', 'PYC Neurons', 'SST Neurons', 'VIP Neurons'))) %>%
  mutate(cluster = factor(cluster, levels = c('PYC Neurons', 'SST Neurons', 'PV Neurons', 'VIP Neurons'))) %>%
  group_by(cluster) %>% mutate(cor = cor(A, B))

normalize_covariates(STRESS.MAPPING$cum.correlation/STRESS.MAPPING$iterations, STRESS$cluster, PRIOR$collapsed) %>% melt %>%
  merge(STRESS@meta.data[, c('sample', 'cluster')], by.x = 'Var1', by.y = 'sample') %>%
  merge(STRESS.MAPPING$stdev.correlation %>% melt(value.name = 'stdev'), by = c('Var1', 'Var2')) %>%
  group_by(Var2, cluster) %>% summarize(A = mean(value), A_err = sqrt(sum(stdev^2)/n()^2)) %>% ungroup() %>%
  mutate(Var2 = factor(Var2, labels = c('Astrocytes', 'Endothelial Cells', 'PYC Neurons', 'Oligodendrocytes', 'PV Neurons', 'SST Neurons', 'VIP Neurons'))) %>%
  mutate(Var2 = factor(Var2, levels = c('Astrocytes', 'Endothelial Cells', 'Oligodendrocytes', 'PYC Neurons', 'SST Neurons', 'PV Neurons', 'VIP Neurons'))) %>%
  cbind(B = 
          normalize_covariates(AGING.MAPPING$cum.correlation/AGING.MAPPING$iterations, AGING$cluster, PRIOR$collapsed) %>% melt %>%
          merge(AGING@meta.data[, c('sample', 'cluster')], by.x = 'Var1', by.y = 'sample') %>%
          merge(AGING.MAPPING$stdev.correlation %>% melt(value.name = 'stdev'), by = c('Var1', 'Var2')) %>%
          group_by(Var2, cluster) %>% summarize(B = mean(value), B_err = sqrt(sum(stdev^2)/n()^2)) %>% ungroup() %>%
          mutate(Var2 = factor(Var2, labels = c('Astrocytes', 'Endothelial Cells', 'PYC Neurons', 'Oligodendrocytes', 'PV Neurons', 'SST Neurons', 'VIP Neurons'))) %>%
          mutate(Var2 = factor(Var2, levels = c('Astrocytes', 'Endothelial Cells', 'Oligodendrocytes', 'PYC Neurons', 'SST Neurons', 'PV Neurons', 'VIP Neurons'))) %>%
          .[, c('B', 'B_err')]) %>% mutate(cluster = factor(cluster, labels = c('PV Neurons', 'PYC Neurons', 'SST Neurons', 'VIP Neurons'))) %>%
  mutate(cluster = factor(cluster, levels = c('PYC Neurons', 'SST Neurons', 'PV Neurons', 'VIP Neurons'))) %>%
  filter(cluster %in% c('SST Neurons')) %>%
  ggplot(aes(A, `B.B`, color = Var2)) + geom_point(size = 3) +
  xlim(0.8, 2.5) + ylim(0.8, 2.5) +
  geom_abline(slope = 1, linetype = 'dashed', color = 'red') +
  scale_color_manual(values = setNames(c('#fec10c', '#9a703e',
                                         '#0c5bb0', '#fa6b0a',
                                         '#ee2210', '#14983d',
                                         '#149bed'),
                                       c('Astrocytes', 'Endothelial Cells', 'Oligodendrocytes',
                                         'PYC Neurons', 'PV Neurons', 'SST Neurons',
                                         'VIP Neurons')), name = 'Mapped Cell Type') +
  theme(legend.position = 'none', strip.background = element_blank(), aspect.ratio = 1) +
  scale_shape(name = 'Actual Cell Type') + xlab('Stress') + ylab(element_blank())

GSE141464 <- getGemma('data/15971_GSE141464_expmat.data.txt')
colnames(GSE141464) <- c(paste0(rep('Saline', 3), paste0('-', 1:3)),
                         paste0(rep('Yoke.Control', 3), paste0('-', 1:3)),
                         paste0(rep('Cocaine', 3), paste0('-', 1:3)),
                         'Saline-4')

GSE141464 <- as_Seurat(Matrix(as.matrix(GSE141464), sparse = T),
                       data.frame(cluster = rep('L[1-6]', 10),
                                  sample = colnames(GSE141464)))

GSE141464.MAPPING <- map_samples(GSE141464@assays$RNA@counts, REFERENCE)

GSE141337 <- getGemma('data/15746_GSE141337_expmat.data.txt')
colnames(GSE141337) <- c('Cocaine-1', 'Saline-1',
                         paste0(rep('Yoke.Control', 3), paste0('-', 1:3)),
                         paste0(rep('Cocaine', 3), paste0('-', 2:4)),
                         paste0(rep('Saline', 3), paste0('-', 2:4)))

GSE141337 <- as_Seurat(Matrix(as.matrix(GSE141337), sparse = T),
                       data.frame(cluster = rep('L[1-6]', 11),
                                  sample = colnames(GSE141337)))

GSE141337.MAPPING <- map_samples(GSE141337@assays$RNA@counts, REFERENCE)

GSE74985 <- getGemma('data/9689_GSE74985_expmat.data.txt')
colnames(GSE74985) <- c(paste0(rep('CA1D', 3), paste0('-', 1:3)),
                        paste0(rep('CA1V', 3), paste0('-', 1:3)),
                        paste0(rep('CA4', 3), paste0('-', 1:3)),
                        paste0(rep('DGD', 3), paste0('-', 1:3)),
                        paste0(rep('DGV', 3), paste0('-', 1:3)),
                        paste0(rep('CA2', 3), paste0('-', 1:3)),
                        paste0(rep('CA3D', 3), paste0('-', 1:3)),
                        paste0(rep('CA3V', 3), paste0('-', 1:3)))

GSE74985 <- as_Seurat(Matrix(as.matrix(GSE74985), sparse = T),
                      data.frame(cluster = c(rep('CA1sp', 6),
                                             rep('DG', 9),
                                             rep('CA2sp', 3),
                                             rep('CA3sp', 6)),
                                 sample = colnames(GSE74985)))

TASIC.HIPPO <- loadReference(regions = 'HIP', clusters = c('DG', 'CA1sp', 'CA2sp', 'CA3sp', 'Astro', 'Vip', 'Pvalb', 'Macrophage', 'Sst', 'Oligo'))
TASIC.HIPPO <- as_Seurat(TASIC.HIPPO$data, TASIC.HIPPO$meta)
REFERENCE.HIPPO <- simplify_reference(TASIC.HIPPO)
PRIOR.HIPPO <- generate_prior(TASIC.HIPPO, REFERENCE.HIPPO, clusters = c('DG', 'CA1sp', 'CA2sp', 'CA3sp', 'Astro', 'Vip', 'Pvalb', 'Macrophage', 'Sst', 'Oligo'))

GSE74985.MAPPING <- map_samples(GSE74985@assays$RNA@counts, REFERENCE.HIPPO, clusters = c('DG', 'CA1sp', 'CA2sp', 'CA3sp', 'Astro', 'Vip', 'Pvalb', 'Macrophage', 'Sst', 'Oligo'))

exp_data <- normalize_covariates(AGING.MAPPING$cum.correlation/AGING.MAPPING$iterations, AGING$cluster, PRIOR$collapsed) %>%
  melt %>% merge(AGING@meta.data[, c('sample', 'cluster')], by.x = 'Var1', by.y = 'sample') %>% cbind(exp = 'LCM-seq') %>%
  rbind(
    normalize_covariates(STRESS.MAPPING$cum.correlation/STRESS.MAPPING$iterations, STRESS$cluster, PRIOR$collapsed) %>%
      melt %>% merge(STRESS@meta.data[, c('sample', 'cluster')], by.x = 'Var1', by.y = 'sample') %>% cbind(exp = 'LCM-seq')
  ) %>%
  rbind(
    normalize_covariates(GSE141337.MAPPING$cum.correlation/GSE141337.MAPPING$iterations, GSE141337$cluster, PRIOR$collapsed) %>%
      melt %>% merge(GSE141337@meta.data[, c('sample', 'cluster')], by.x = 'Var1', by.y = 'sample') %>% cbind(exp = 'TRAP-seq')
  ) %>%
  rbind(
    normalize_covariates(GSE141464.MAPPING$cum.correlation[1:9,]/GSE141464.MAPPING$iterations, GSE141464$cluster[1:9], PRIOR$collapsed) %>%
      melt %>% merge(GSE141464@meta.data[, c('sample', 'cluster')], by.x = 'Var1', by.y = 'sample') %>% cbind(exp = 'TRAP-seq')
  ) %>%
  rbind(
    normalize_covariates(GSE74985.MAPPING$cum.correlation/GSE74985.MAPPING$iterations, GSE74985$cluster, PRIOR.HIPPO$collapsed) %>%
      melt %>% merge(GSE74985@meta.data[, c('sample', 'cluster')], by.x = 'Var1', by.y = 'sample') %>% cbind(exp = 'Manual Sorting')
  ) %>%
  filter(Var2 %in% c('Astro', 'Oligo', 'Sst'), cluster %in% c('L[1-6]', 'Pyr', 'PYC Neurons', 'CA1sp', 'CA2sp', 'CA3sp')) %>%
  mutate(Var2 = factor(Var2, labels = c('Astrocytes', 'Oligodendrocytes', 'SST Neurons'))) %>%
  mutate(exp = factor(exp, levels = c('scRNA-seq', 'Manual Sorting', 'LCM-seq', 'TRAP-seq')))

t.test(exp_data %>% filter(exp == 'Manual Sorting' & Var2 == 'Astrocytes') %>% .[, 'value'],
       exp_data %>% filter(exp == 'TRAP-seq' & Var2 == 'Astrocytes') %>% .[, 'value'])
t.test(exp_data %>% filter(exp == 'Manual Sorting' & Var2 == 'Astrocytes') %>% .[, 'value'],
       exp_data %>% filter(exp == 'LCM-seq' & Var2 == 'Astrocytes') %>% .[, 'value'])
t.test(exp_data %>% filter(exp == 'LCM-seq' & Var2 == 'Astrocytes') %>% .[, 'value'],
       exp_data %>% filter(exp == 'TRAP-seq' & Var2 == 'Astrocytes') %>% .[, 'value'])

wilcox.test(exp_data %>% filter(exp == 'Manual Sorting' & Var2 == 'Oligodendrocytes') %>% .[, 'value'],
       exp_data %>% filter(exp == 'TRAP-seq' & Var2 == 'Oligodendrocytes') %>% .[, 'value'])
wilcox.test(exp_data %>% filter(exp == 'Manual Sorting' & Var2 == 'Oligodendrocytes') %>% .[, 'value'],
       exp_data %>% filter(exp == 'LCM-seq' & Var2 == 'Oligodendrocytes') %>% .[, 'value'])
wilcox.test(exp_data %>% filter(exp == 'LCM-seq' & Var2 == 'Oligodendrocytes') %>% .[, 'value'],
       exp_data %>% filter(exp == 'TRAP-seq' & Var2 == 'Oligodendrocytes') %>% .[, 'value'])

t.test(exp_data %>% filter(exp == 'Manual Sorting' & Var2 == 'SST Neurons') %>% .[, 'value'],
       exp_data %>% filter(exp == 'TRAP-seq' & Var2 == 'SST Neurons') %>% .[, 'value'])
t.test(exp_data %>% filter(exp == 'Manual Sorting' & Var2 == 'SST Neurons') %>% .[, 'value'],
       exp_data %>% filter(exp == 'LCM-seq' & Var2 == 'SST Neurons') %>% .[, 'value'])
t.test(exp_data %>% filter(exp == 'LCM-seq' & Var2 == 'SST Neurons') %>% .[, 'value'],
       exp_data %>% filter(exp == 'TRAP-seq' & Var2 == 'SST Neurons') %>% .[, 'value'])

exp_data %>%
  ggplot(aes(exp, value, color = Var2)) + 
  geom_violin() + geom_point(pch = 21, position = position_jitterdodge()) +
  facet_wrap(~Var2, scales = 'free_y') +
  theme(legend.justification = c(1, 1),#, legend.position = c(1, 1),
        strip.background = element_blank(), strip.text = element_text(hjust = 0),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_color_manual(name = 'Off-target Cell Type', values = setNames(c('#fec10c', '#9a703e',
                                         '#0c5bb0', '#fa6b0a',
                                         '#ee2210', '#14983d',
                                         '#149bed'),
                                       c('Astrocytes', 'Endo', 'Oligodendrocytes',
                                         'PYC Neurons', 'Pvalb', 'SST Neurons',
                                         'VIP Neurons'))) +
  xlab(element_blank()) + ylab('Scaled Contamination Coefficient') + scale_y_continuous(expand = c(0, 0))

get_DEG.variants <- function(data, contam.formula = ~`L.1.6.` + Sst + Vip + Pvalb + condition) {
  list(noqc = get_DEGs(data$contaminated, ~condition),
       raw  = get_DEGs(data$contaminated, contam.formula, data$mapping.contaminated$cum.correlation / data$mapping.contaminated$iterations))
}

lapply(c(0, 0.05, 0.1, 0.2), function(contam) {
  lapply(c(0.5, 0.6, 0.7, 1), function(confound) {
    if(contam == 0 && confound != 0.5) return(NULL)
    
    run_synthetic(celltypes = c('L[1-6]', 'Sst'), ref.collapsed = REFERENCE,
                  prior = PRIOR$collapsed,
                  contamination = matrix(rep(contam, 4), byrow = T, nrow = 2), confound = confound, n.genes = 20000, effect.size = 3, n.DEGs = 3000)
  })
}) -> SYNTHETIC

lapply(1:4, function(A) {
  lapply(1:4, function(B) {
    if(A == 1 && B != 1) return(NULL)
    
    get_DEG.variants(SYNTHETIC[[A]][[B]], ~`L.1.6.` + Sst + condition)
  })
}) -> DEGs

thresh <- 0.1
n <- c(4,4)

DEGs[[n[1]]][[n[2]]]$noqc$Sst %>% results %>% as.data.table(keep.rownames = T) %>%
  merge(SYNTHETIC[[n[1]]][[n[2]]]$meta$Sst %>% as.data.table(keep.rownames = T) %>% .[, .(rn, differential.expression)],
        by = 'rn', all = T) %>%
  .[grepl('g[0-9]', rn) & !is.na(differential.expression), .(rn, V2 = padj, differential.expression)] %>%
  mutate(V2 = ifelse(is.na(V2), 1, V2)) %>%
  summarize(TP = sum(V2 < thresh & differential.expression == 1)/sum(V2 < thresh),
            FPF = sum(V2 < thresh & differential.expression == 0)/sum(V2 < thresh),
            TN = sum(V2 >= thresh & differential.expression == 0)/sum(V2 >= thresh),
            FNF = sum(V2 >= thresh & differential.expression == 1)/sum(V2 >= thresh))
DEGs[[n[1]]][[n[2]]]$raw$Sst %>% results %>% as.data.table(keep.rownames = T) %>%
  merge(SYNTHETIC[[n[1]]][[n[2]]]$meta$Sst %>% as.data.table(keep.rownames = T) %>% .[, .(rn, differential.expression)],
        by = 'rn', all = T) %>%
  .[grepl('g[0-9]', rn) & !is.na(differential.expression), .(rn, V2 = padj, differential.expression)] %>%
  mutate(V2 = ifelse(is.na(V2), 1, V2)) %>%
  summarize(TP = sum(V2 < thresh & differential.expression == 1)/sum(V2 < thresh),
            FPF = sum(V2 < thresh & differential.expression == 0)/sum(V2 < thresh),
            TN = sum(V2 >= thresh & differential.expression == 0)/sum(V2 >= thresh),
            FNF = sum(V2 >= thresh & differential.expression == 1)/sum(V2 >= thresh))

lapply(1:length(DEGs.AGING), function(i) {
  data.frame(exp = 'Aging',
             cell = names(DEGs.AGING)[i],
             data.frame(QC = DEGs.AGING[[i]] %>% mcols %>% .[, 'deviance'] %>% sum %>% `+`(5),
                        raw = DEGs.AGING.noqc[[i]] %>% mcols %>% .[, 'deviance'] %>% sum %>% `+`(1))
  ) %>% rbind(
    data.frame(exp = 'Stress',
               cell = names(DEGs.AGING)[i],
               data.frame(QC = DEGs.STRESS[[i]] %>% mcols %>% .[, 'deviance'] %>% sum %>% `+`(5),
                          raw = DEGs.STRESS.noqc[[i]] %>% mcols %>% .[, 'deviance'] %>% sum %>% `+`(1))
    )
  )
}) %>% rbindlist

lapply(1:length(DEGs.AGING), function(i) {
  data.frame(exp = 'Aging',
             cell = names(DEGs.AGING)[i],
             DEGs.AGING[[i]] %>% mcols %>% .[, 'deviance'] %>% `+`(5) %>%
               cbind(
                 DEGs.AGING.noqc[[i]] %>% mcols %>% .[, 'deviance'] %>% `+`(1)
               ) %>% { .[, 1] - .[, 2] } %>%
               { data.frame(better = mean(.[. < 0]), worse = mean(.[. > 0])) }
               # { data.frame(QC = sum(. < 0), raw = sum(. > 0)) }
  ) %>% rbind(
    data.frame(exp = 'Stress',
               cell = names(DEGs.AGING)[i],
               DEGs.STRESS[[i]] %>% mcols %>% .[, 'deviance'] %>% `+`(5) %>%
                 cbind(
                   DEGs.STRESS.noqc[[i]] %>% mcols %>% .[, 'deviance'] %>% `+`(1)
                 ) %>% { .[, 1] - .[, 2] } %>%
                 { data.frame(better = mean(.[. < 0]), worse = mean(.[. > 0])) }
                 # { data.frame(QC = sum(. < 0), raw = sum(. > 0)) }
    )
  )
}) %>% rbindlist %>% mutate(efficacy = QC/(raw+QC)) %>% group_by(exp) %>% summarize(sum(QC)/sum(QC+raw))