source('R/analysis_funs.R')
library(ggpubr)

## set gene list and GO file pathway
sperm_list
oocyte_list
GO_file <- 'data/msigdb/mh.all.v2023.1.Mm.symbols.gmt'

### sperm
sperm_GSEA1 <- my_GSEA(sperm_list[[1]], GO_file, min_size = 15, collapse = TRUE, pval = 0.05, lfc_cutoff = log(1.5, 2))
sperm_GSEA1$Results
sp1 <- sperm_GSEA1$Plot + labs(title = 'Sperm', caption = 'Hu-FS-400ppb vs. Hu-FD-400ppb')

sperm_GSEA2 <- my_GSEA(sperm_list[[2]], GO_file, min_size = 15, collapse = TRUE, pval = 0.05, lfc_cutoff = log(1.5, 2))
sperm_GSEA2$Results
sperm_GSEA2$Plot

sperm_GSEA3 <- my_GSEA(sperm_list[[3]], GO_file, min_size = 15, collapse = TRUE, pval = 0.05, lfc_cutoff = log(1.5, 2))
sperm_GSEA3$Results
sperm_GSEA3$Plot

sperm_GSEA4 <- my_GSEA(sperm_list[[4]], GO_file, min_size = 15, collapse = TRUE, pval = 0.05, lfc_cutoff = log(1.5, 2))
sperm_GSEA4$Results
sp4 <- sperm_GSEA4$Plot + labs(title = 'Sperm', caption = 'Hu-FA-400ppb vs. Hu-FD-400ppb')

sperm_GSEA5 <- my_GSEA(sperm_list[[5]], GO_file, min_size = 15, collapse = TRUE, pval = 0.05, lfc_cutoff = log(1.5, 2))
sperm_GSEA5$Results
sperm_GSEA5$Plot

sperm_GSEA6 <- my_GSEA(sperm_list[[6]], GO_file, min_size = 15, collapse = TRUE, pval = 0.05, lfc_cutoff = log(1.5, 2))
sperm_GSEA6$Results
sperm_GSEA6$Plot

### oocyte
oocyte_GSEA1 <- my_GSEA(oocyte_list[[1]], GO_file, min_size = 15, collapse = TRUE, pval = 0.2)
oocyte_GSEA1$Results
oocyte_GSEA1$Plot

oocyte_GSEA2 <- my_GSEA(oocyte_list[[2]], GO_file, min_size = 15, collapse = TRUE, pval = 0.05, lfc_cutoff = log(1.5, 2))
oocyte_GSEA2$Results
oocyte_GSEA2$Plot

oocyte_GSEA3 <- my_GSEA(oocyte_list[[3]], GO_file, min_size = 15, collapse = TRUE, pval = 0.05, lfc_cutoff = log(1.5, 2))
oocyte_GSEA3$Results
oocyte_GSEA3$Plot

oocyte_GSEA4 <- my_GSEA(oocyte_list[[4]], GO_file, min_size = 15, collapse = TRUE, pval = 0.05, lfc_cutoff = log(1.5, 2))
oocyte_GSEA4$Results
op4 <- oocyte_GSEA4$Plot + labs(title = 'Oocyte', caption = 'Hu-FA-400ppb vs. Hu-FD-400ppb')

oocyte_GSEA5 <- my_GSEA(oocyte_list[[5]], GO_file, min_size = 15, collapse = TRUE, pval = 0.05, lfc_cutoff = log(1.5, 2))
oocyte_GSEA5$Results
oocyte_GSEA5$Plot

oocyte_GSEA6 <- my_GSEA(oocyte_list[[6]], GO_file, min_size = 15, collapse = TRUE, pval = 0.05, lfc_cutoff = log(1.5, 2))
oocyte_GSEA6$Results
op6 <- oocyte_GSEA6$Plot + labs(title = 'Oocyte', caption = 'Hu-FA-400ppb vs. Hu-FD-0ppb')
ggarrange(sp1, sp4, op4, op6, nrow = 2, ncol = 2, common.legend = TRUE, legend = 'bottom')
ggsave('GSEA_plot.png', width = 10, height = 6)
