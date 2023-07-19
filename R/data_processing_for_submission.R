library(seqUtils)
library(tidyverse)
library(readxl)
library(janitor)

# Data processing, QC, and prefiltering ---------------------------------

## Raw data file name is omitted

## rename the 4 unmatched samples (for sperms)
sperm_sample_names <- colnames(fpfm_sperm_seq)
sperm_sample_names[which(sperm_sample_names == 'FPFM_76A_111_400ppb')] <- 'FPFM_76A_1111_400ppb'
sperm_sample_names[which(sperm_sample_names == 'FPFM_78_1059_0ppb')] <- 'FPFM_78A_1059_0ppb'
sperm_sample_names[which(sperm_sample_names == 'FPFM_49_N_400ppb')] <- 'FPFM_49_RN_400ppb'
sperm_sample_names[which(sperm_sample_names == 'FPFM_50_101X_400ppb')] <- 'FPFM_50_1016_400ppb'
colnames(fpfm_sperm_seq) <- sperm_sample_names

## load the handwritten data spreadsheet
fpfm_col_raw$Sex <- relevel(fpfm_col_raw$Sex, ref = 'M')
fpfm_col_raw$Strain <- relevel(fpfm_col_raw$Strain, ref = 'WT')
fpfm_col_raw$Diet <- relevel(fpfm_col_raw$Diet, ref = 'FA')
fpfm_col_raw$Dose <- relevel(fpfm_col_raw$Dose, ref = '0ppb')
fpfm_col_raw <- fpfm_col_raw %>% 
  filter(!duplicated(Sample))

## create column data for samples that are sequenced (n = 5 for each group)
fpfm_sperm_col <- tibble(Sample = colnames(fpfm_sperm_seq)[-1]) %>% 
  filter(Sample != 'SUM') %>% 
  left_join(fpfm_col_raw, by = 'Sample') %>% 
  distinct(Sample, .keep_all = TRUE) 
fpfm_sperm_col %>% 
  print(n = 60) # check if anything mouse doesn't have coldata

fpfm_sperm_col %>% 
  group_by(Strain, Dose, Diet) %>% 
  summarise(n = n())

## also read Bingzhen's coldata to take a look
fpfm_oocyte_col <- bz_col %>% 
  select(-`n=57`) %>% 
  dplyr::rename('Strain' = 'Genotype',
         'Dose' = 'Treatment',
         'Diet' = 'Folate') %>% 
  mutate(Cage = gsub("(.*_){1}(\\d+)(_.*){2}", "\\2", Sample)) %>% 
  mutate_at(2:6, as.factor) 
fpfm_oocyte_col$Dose <- recode_factor(fpfm_oocyte_col$Dose, 
                                      '400' = '400ppb',
                                      '40' = '40ppb',
                                      '0' = '0ppb')
fpfm_oocyte_col$Strain <- recode_factor(fpfm_oocyte_col$Strain, 
                                        'Hu' = 'HU')
fpfm_oocyte_col$Strain <- relevel(fpfm_oocyte_col$Strain, ref = 'WT')
fpfm_oocyte_col$Diet <- relevel(fpfm_oocyte_col$Diet, ref = 'FA')
fpfm_oocyte_col$Dose <- relevel(fpfm_oocyte_col$Dose, ref = '0ppb')

## inspecting the unexpected number of rows in the merged dataset; comment out when running the analysis
# setdiff(colnames(fpfm_seq)[-1], fpfm_sperm_col$Sample)
# tbl <- table(fpfm_sperm_col$Sample)
# fpfm_sperm_col %>% 
#   filter(Sample %in% names(tbl[tbl > 1]))

## make regular data
fpfm_sperm_d <- as.matrix(fpfm_sperm_seq[, -c(1, 61)])
rownames(fpfm_sperm_d) <- fpfm_sperm_seq$Gene
fpfm_oocyte_d <- as.matrix(fpfm_oocyte_seq[, -1])
rownames(fpfm_oocyte_d) <- fpfm_oocyte_seq$gene

## pre-filer
fpfm_sperm_d_filter <- man.filter(fpfm_sperm_d)
fpfm_oocyte_d_filter <- man.filter(fpfm_oocyte_d)

## realign oocyte samples with the coldata
fpfm_oocyte_d_filter <- fpfm_oocyte_d_filter[, fpfm_oocyte_col$Sample]

## make QC plots
## choose 2 samples from each combination of conditions to plot
fpfm_sperm_to_plt <- fpfm_sperm_col %>% 
  group_by(Strain, Diet, Dose) %>% 
  sample_n(1) %>% 
  ungroup()
fpfm_oocyte_to_plt <- fpfm_oocyte_col %>% 
  group_by(Strain, Diet, Dose) %>% 
  sample_n(1) %>% 
  ungroup()
  
## boxplots
fpfm_sperm_bp <- expr.bp(Y = fpfm_sperm_d_filter[, which(colnames(fpfm_sperm_d_filter) %in% fpfm_sperm_to_plt$Sample)],
                         n = nrow(fpfm_sperm_to_plt))
fpfm_oocyte_bp <- expr.bp(Y = fpfm_oocyte_d_filter[, fpfm_oocyte_to_plt$Sample],
                          n = nrow(fpfm_oocyte_to_plt))

## pairwise scatterplots
pw.scatter(Y = fpfm_sperm_d_filter[, which(colnames(fpfm_sperm_d_filter) %in% fpfm_sperm_to_plt$Sample)], 
           notes = 'FPFM_sperm', 
           dir_name = 'plots')
pw.scatter(Y = fpfm_oocyte_d_filter[, fpfm_oocyte_to_plt$Sample], 
           notes = 'FPFM_oocyte', 
           dir_name = 'plots')

## make a pairwise heatmap for all the mice
Y <- fpfm_sperm_d_filter
coldata <- fpfm_sperm_col %>% 
  arrange(Sex, Strain, Diet, Dose)
Y_rearrange <- Y[, match(coldata$Sample, colnames(Y))]
cor_mat <- cor(Y_rearrange)
rownames(cor_mat) <- colnames(cor_mat) <- colnames(Y_rearrange)
cor_df <- reshape2::melt(t(cor_mat)) %>% as_tibble()
cor_plot <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  # scale_fill_gradient(low = 'yellow', high = 'red')
  viridis::scale_fill_viridis()
pdf('plots/FPFM_sperm_cor_matrx_heatmap.pdf')
print(cor_plot)
dev.off()

## write the pairwise heatmap function
coldata <- fpfm_oocyte_col %>% 
  arrange(Sex, Strain, Diet, Dose)
Y <- fpfm_oocyte_d_filter[, coldata$Sample]
dir_name <- 'plots'
notes <- 'FPFM_oocyte_cor_matrix'
make.pw.heatmap <- function(Y, dir_name, notes) {
  cor_mat <- cor(Y)
  rownames(cor_mat) <- colnames(cor_mat) <- colnames(Y)
  cor_df <- reshape2::melt(t(cor_mat)) %>% as_tibble()
  cor_plot <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile() +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    # scale_fill_gradient(low = 'yellow', high = 'red')
    viridis::scale_fill_viridis()
  pdf(paste0(dir_name, '/', notes, '_heatmap.pdf'))
  print(cor_plot)
  dev.off()
}
make.pw.heatmap(Y, dir_name, notes)

## save the final data
## check the colnames of the sequencing data match those of the sample coldata
all.equal(colnames(fpfm_sperm_d_filter), fpfm_sperm_col$Sample) # proceed if TRUE
all.equal(colnames(fpfm_oocyte_d_filter), fpfm_oocyte_col$Sample) # proceed if TRUE
