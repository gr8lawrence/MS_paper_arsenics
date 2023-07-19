library(DESeq2)
library(tidyverse)
library(patchwork)
library(ggvenn)
library(ggVennDiagram)
source('R/analysis_funs.R')

## Processed data file and excluded sample names are omitted

### sperms
## combine all experimental conditions into different groups
sperm_col_fpfm <- sperm_col_fpfm %>% 
  mutate(Group = as.factor(paste(Strain, Diet, Dose, sep = '_'))) %>% 
  mutate(Group = fct_relevel(Group, 'WT_FD_0ppb', after = 0)) # set WT_FD_0ppb to the reference level

## fit the DEseq analysis
dds_sperm <- DESeqDataSetFromMatrix(countData = sperm_counts_fpfm,
                                    colData = sperm_col_fpfm,
                                    design = ~ Group)
dds_sperm <- DESeq(dds_sperm, test = 'Wald')
resultsNames(dds_sperm)

## pairwise comparisons (1-6)
sperm1 <- results(dds_sperm, contrast = c('Group', 'HU_FS_400ppb', 'HU_FD_400ppb')) %>% 
  as_tibble(rownames = 'geneSymbol') %>% 
  arrange(padj)
sperm2 <- results(dds_sperm, contrast = c('Group', 'HU_FS_0ppb', 'HU_FD_0ppb')) %>% 
  as_tibble(rownames = 'geneSymbol') %>% 
  arrange(padj)
sperm3 <- results(dds_sperm, contrast = c('Group', 'HU_FS_400ppb', 'HU_FD_0ppb')) %>% 
  as_tibble(rownames = 'geneSymbol') %>% 
  arrange(padj)
sperm4 <- results(dds_sperm, contrast = c('Group', 'HU_FA_400ppb', 'HU_FD_400ppb')) %>% 
  as_tibble(rownames = 'geneSymbol') %>% 
  arrange(padj)
sperm5 <- results(dds_sperm, contrast = c('Group', 'HU_FA_0ppb', 'HU_FD_0ppb')) %>% 
  as_tibble(rownames = 'geneSymbol') %>% 
  arrange(padj)
sperm6 <- results(dds_sperm, contrast = c('Group', 'HU_FA_400ppb', 'HU_FD_0ppb')) %>% 
  as_tibble(rownames = 'geneSymbol') %>% 
  arrange(padj)

## provide a lists of genes
sperm_list <- list(6)
sperm_list[[1]] <- get_gene_list(sperm1 %>% arrange(desc(log2FoldChange)))
sperm_list[[2]] <- get_gene_list(sperm2 %>% arrange(desc(log2FoldChange)))
sperm_list[[3]] <- get_gene_list(sperm3 %>% arrange(desc(log2FoldChange)))
sperm_list[[4]] <- get_gene_list(sperm4 %>% arrange(desc(log2FoldChange)))
sperm_list[[5]] <- get_gene_list(sperm5 %>% arrange(desc(log2FoldChange)))
sperm_list[[6]] <- get_gene_list(sperm6 %>% arrange(desc(log2FoldChange)))

### oocytes
## combine all experimental conditions
oocyte_col_fpfm <- oocyte_col_fpfm %>% 
  mutate(Group = as.factor(paste(Strain, Diet, Dose, sep = '_'))) %>% 
  mutate(Group = fct_relevel(Group, 'WT_FD_0ppb', after = 0)) # set WT_FD_0ppb to the reference level

## fit the DEseq analysis
dds_oocyte <- DESeqDataSetFromMatrix(countData = oocyte_counts_fpfm,
                                     colData = oocyte_col_fpfm,
                                     design = ~ Group)
dds_oocyte <- DESeq(dds_oocyte, test = 'Wald')
resultsNames(dds_oocyte)

## pairwise comparisons (1-6)
oocyte1 <- results(dds_oocyte, contrast = c('Group', 'HU_FS_400ppb', 'HU_FD_400ppb')) %>% 
  as_tibble(rownames = 'geneSymbol') %>% 
  arrange(padj)
oocyte2 <- results(dds_oocyte, contrast = c('Group', 'HU_FS_0ppb', 'HU_FD_0ppb')) %>% 
  as_tibble(rownames = 'geneSymbol') %>% 
  arrange(padj)
oocyte3 <- results(dds_oocyte, contrast = c('Group', 'HU_FS_400ppb', 'HU_FD_0ppb')) %>% 
  as_tibble(rownames = 'geneSymbol') %>% 
  arrange(padj)
oocyte4 <- results(dds_oocyte, contrast = c('Group', 'HU_FA_400ppb', 'HU_FD_400ppb')) %>% 
  as_tibble(rownames = 'geneSymbol') %>% 
  arrange(padj)
oocyte5 <- results(dds_oocyte, contrast = c('Group', 'HU_FA_0ppb', 'HU_FD_0ppb')) %>% 
  as_tibble(rownames = 'geneSymbol') %>% 
  arrange(padj)
oocyte6 <- results(dds_oocyte, contrast = c('Group', 'HU_FA_400ppb', 'HU_FD_0ppb')) %>% 
  as_tibble(rownames = 'geneSymbol') %>% 
  arrange(padj)

## provide a lists of genes
oocyte_list <- list(6)
oocyte_list[[1]] <- get_gene_list(oocyte1 %>% arrange(desc(log2FoldChange)))
oocyte_list[[2]] <- get_gene_list(oocyte2 %>% arrange(desc(log2FoldChange)))
oocyte_list[[3]] <- get_gene_list(oocyte3 %>% arrange(desc(log2FoldChange)))
oocyte_list[[4]] <- get_gene_list(oocyte4 %>% arrange(desc(log2FoldChange)))
oocyte_list[[5]] <- get_gene_list(oocyte5 %>% arrange(desc(log2FoldChange)))
oocyte_list[[6]] <- get_gene_list(oocyte6 %>% arrange(desc(log2FoldChange)))



