## get significant DE results in the tibble format
get_tidy_df <- function(res, gene_names, cutoff = 0.1) {
  res <- as_tibble(res) %>% 
    mutate(gene = gene_names) %>% 
    arrange(padj) %>% 
    filter(padj < cutoff) %>% 
    select(c(7, 1:6))
  res
}
get_tidy_df_top_genes <- function(res, gene_names, top = 20) {
  res <- as_tibble(res) %>% 
    mutate(gene = gene_names) %>% 
    arrange(padj) %>% 
    filter(row_number() <= top) %>% 
    select(c(7, 1:6))
  res
}

## function to get the stratified analysis
# counts = sperm_counts_fpfm
# cols = sperm_col_fpfm
# strain = 'HU'
# dose = '0ppb'
# tissue = 'sperm'

folate_strat_dose <- function(counts, cols, strain, dose, tissue = c('sperm', 'oocyte')) {
  cols1 = cols %>% 
    filter(Strain == strain & Dose == dose)
  counts1 <- counts[, cols1$Sample]
  # all.equal(colnames(sperm_counts1), sperm_cols1$Sample) # TRUE
  dds_s1 <- DESeqDataSetFromMatrix(countData = counts1,
                                   colData = cols1,
                                   design = ~ Diet)
  dds_s1 <- DESeq(dds_s1)
  gene_names <- rownames(counts1)
  fs_v_fa <- results(dds_s1, contrast = c('Diet', 'FS', 'FA'))
  fs_v_fd <- results(dds_s1, contrast = c('Diet', 'FS', 'FD'))
  fa_v_fd <- results(dds_s1, contrast = c('Diet', 'FA', 'FD'))
  
  ## extract mean normalized counts for each group
  nc <- counts(dds_s1, normalized = TRUE)
  nc_df <- as_tibble(nc) %>% 
    mutate(Gene = gene_names) %>% 
    pivot_longer(seq(ncol(nc)), names_to = 'Sample', values_to = 'Normalized_count') %>% 
    left_join(cols, by = 'Sample') %>% 
    select(-Sex, -Strain, -Dose, -Cage) %>% 
    group_by(Gene, Diet) %>% 
    summarise(Mean_nc = mean(Normalized_count)) %>% 
    ungroup %>% 
    pivot_wider(id_cols = 'Gene', names_from = 'Diet', values_from = 'Mean_nc')
  rownames(nc_df) <- nc_df$Gene
  nc_df <- nc_df[gene_names, ]
  # all.equal(gene_names, nc_df$Gene) # should be TRUE
  
  ## Output1: a list of normalized counts, adjusted, and unadjusted p-values for each gene
  p_val_df <- tibble(Gene = gene_names,
                     nc_FD = nc_df$FD,
                     nc_FA = nc_df$FA,
                     nc_FS = nc_df$FS,
                     FS_FA_lfc = fs_v_fa$log2FoldChange,
                     FS_FA_p_raw = fs_v_fa$pvalue,
                     FS_FA_p_adj = fs_v_fa$padj,
                     FS_FD_lfc = fs_v_fd$log2FoldChange,
                     FS_FD_p_raw = fs_v_fd$pvalue,
                     FS_FD_p_adj = fs_v_fd$padj,
                     FA_FD_lfc = fa_v_fd$log2FoldChange,
                     FA_FD_p_raw = fa_v_fd$pvalue,
                     FA_FD_p_adj = fa_v_fd$padj)
  write_csv(p_val_df, paste0('new_results/', strain, '_', tissue, '_', dose, '_pairwise_diet_p_values.csv'))
  
  ## get the significant genes
  fs_v_fa_sig <- get_tidy_df(fs_v_fa, gene_names)
  fs_v_fd_sig <- get_tidy_df(fs_v_fd, gene_names)
  fa_v_fd_sig <- get_tidy_df(fa_v_fd, gene_names)
  # write_csv(fs_v_fa_sig, paste('results/stratified', strain, dose, 'FS_v_FA_sig.csv', sep = '_'))
  # write_csv(fs_v_fd_sig, paste('results/stratified', strain, dose, 'FS_v_FD_sig.csv', sep = '_'))
  # write_csv(fa_v_fd_sig, paste('results/stratified', strain, dose, 'FA_v_FD_sig.csv', sep = '_'))
  res_ls <- list(FS_v_FA = fs_v_fa_sig,
                 FS_v_FD = fs_v_fd_sig,
                 FA_v_FD = fa_v_fd_sig,
                 dds = dds_s1,
                 col_data = cols1,
                 nc_df = nc_df,
                 nc_full = nc)
  res_ls
}

# counts = sperm_counts_fpfm
# cols = sperm_col_fpfm
# strain = 'HU'
# diet = 'FD'
# tissue = 'sperm'

folate_strat_diet <- function(counts, cols, strain, diet, tissue = c('sperm', 'oocyte')) {
  cols1 <- cols %>% 
    filter(Strain == strain & Diet == diet)
  counts1 <- counts[, cols1$Sample]
  # all.equal(colnames(sperm_counts1), sperm_cols1$Sample) # TRUE
  dds_s1 <- DESeqDataSetFromMatrix(countData = counts1,
                                   colData = cols1,
                                   design = ~ Dose)
  dds_s1 <- DESeq(dds_s1)
  gene_names <- rownames(counts1)
  res_400_v_0 <- results(dds_s1, contrast = c('Dose', '400ppb', '0ppb'))
  
  ## extract mean normalized counts for each group
  nc <- counts(dds_s1, normalized = TRUE)
  nc_df <- as_tibble(nc) %>% 
    mutate(Gene = gene_names) %>% 
    pivot_longer(seq(ncol(nc)), names_to = 'Sample', values_to = 'Normalized_count') %>% 
    left_join(cols, by = 'Sample') %>% 
    select(-Sex, -Strain, -Diet, -Cage) %>% 
    group_by(Gene, Dose) %>% 
    summarise(Mean_nc = mean(Normalized_count)) %>% 
    ungroup %>% 
    pivot_wider(id_cols = 'Gene', names_from = 'Dose', values_from = 'Mean_nc')
  rownames(nc_df) <- nc_df$Gene
  nc_df <- nc_df[gene_names, ]
  # all.equal(gene_names, nc_df$Gene) # should be TRUE
  
  ## Output1: a list of adjusted and unadjusted p-values for each gene
  p_val_df <- tibble(Gene = gene_names,
                     nc_0ppb = nc_df[['0ppb']],
                     nc_400ppb = nc_df[['400ppb']],
                     lfc = res_400_v_0$log2FoldChange,
                     p_raw = res_400_v_0$pvalue,
                     p_adj = res_400_v_0$padj)
  write_csv(p_val_df, paste0('new_results/', strain, '_', tissue, '_', diet, '_pairwise_dose_p_values.csv'))
  
  ## get the significant genes
  sig_400_v_0 <- get_tidy_df(res_400_v_0, gene_names)
  res_ls <- list(sig_400_v_0 = sig_400_v_0,
                 dds = dds_s1,
                 col_data = cols1,
                 nc_df = nc_df,
                 nc_full = nc)
  res_ls
}

# counts = sperm_counts_fpfm
# cols = sperm_col_fpfm
# dose = '400ppb'
# diet = 'FD'
# tissue = 'sperm'

folate_strat_strain <- function(counts, cols, dose, diet, tissue = c('sperm', 'oocyte')) {
  cols1 = cols %>% 
    filter(Dose == dose & Diet == diet)
  counts1 <- counts[, cols1$Sample]
  # all.equal(colnames(sperm_counts1), sperm_cols1$Sample) # TRUE
  dds_s1 <- DESeqDataSetFromMatrix(countData = counts1,
                                   colData = cols1,
                                   design = ~ Strain)
  dds_s1 <- DESeq(dds_s1)
  gene_names <- rownames(counts1)
  res_HU_v_WT <- results(dds_s1, contrast = c('Strain', 'HU', 'WT'))
  
  ## extract mean normalized counts for each group
  nc <- counts(dds_s1, normalized = TRUE)
  nc_df <- as_tibble(nc) %>% 
    mutate(Gene = gene_names) %>% 
    pivot_longer(seq(ncol(nc)), names_to = 'Sample', values_to = 'Normalized_count') %>% 
    left_join(cols, by = 'Sample') %>% 
    select(-Sex, -Dose, -Diet, -Cage) %>% 
    group_by(Gene, Strain) %>% 
    summarise(Mean_nc = mean(Normalized_count)) %>% 
    ungroup %>% 
    pivot_wider(id_cols = 'Gene', names_from = 'Strain', values_from = 'Mean_nc')
  rownames(nc_df) <- nc_df$Gene
  nc_df <- nc_df[gene_names, ]
  # all.equal(gene_names, nc_df$Gene) # should be TRUE
  
  ## Output1: a list of adjusted and unadjusted p-values for each gene
  p_val_df <- tibble(Gene = gene_names,
                     nc_WT = nc_df[['WT']],
                     nc_HU = nc_df[['HU']],
                     lfc = res_HU_v_WT$log2FoldChange,
                     p_raw = res_HU_v_WT$pvalue,
                     p_adj = res_HU_v_WT$padj)
  write_csv(p_val_df, paste0('new_results/', dose, '_', tissue, '_', diet, '_pairwise_strain_p_values.csv'))
  
  ## get the significant genes
  sig_HU_v_WT <- get_tidy_df(res_HU_v_WT, gene_names)
  res_ls <- list(sig_HU_v_WT = sig_HU_v_WT,
                 dds = dds_s1,
                 col_data = cols1,
                 nc_df = nc_df,
                 nc_full = nc)
  res_ls
}


## function for making the 


## function for getting the significant genes in One-way ANOVA
folate_one_way <- function(counts, coldata, strain, cutoff = 0.1) {
  
  ## subset coldata and counts to specified strain
  coldata2 <- coldata %>% 
    filter(Strain == strain)
  counts2 <- counts[, coldata2$Sample]
  if (!all.equal(colnames(counts2), coldata2$Sample)) message('Sample ordering does not match between count matrix and coldata.')
  
  ## get DESeq model fit
  dds_s2 <- DESeqDataSetFromMatrix(countData = counts2,
                                   colData = coldata2,
                                   design = ~ DDG)
  dds_s2 <- DESeq(dds_s2, test = 'LRT', reduced = ~ 1) # compare to a model with only the intercept using the LRT
  
  ## obtain a list of significant DE genes (cutoff padj = 0.1)
  gene_names <- rownames(counts2)
  res2 <- results(dds_s2) %>% 
    get_tidy_df(gene_names = gene_names,
                cutoff = cutoff)
  res_ls <- list(dds = dds_s2, cols = coldata2, res = res2)
  res_ls
}

## function to just get the group means
# nc = sp_nc
# cols = sperm_col_fpfm
get_means <- function(nc, cols) {
  gene_names <- rownames(nc)
  nc_df <- as_tibble(nc) %>% 
    mutate(Gene = gene_names) %>% 
    pivot_longer(seq(ncol(nc)), names_to = 'Sample', values_to = 'Normalized_count') %>% 
    left_join(cols, by = 'Sample') %>% 
    select(-Sex, -Strain, -Diet, -Dose, -Cage) %>% 
    group_by(Gene, DDG) %>% 
    summarise(Mean_nc= mean(Normalized_count)) %>% 
    ungroup %>% 
    pivot_wider(id_cols = 'Gene', names_from = 'DDG', values_from = 'Mean_nc')
  nc_df
}

## Make mean normalized count table for those genes
# strain = 'HU'
# dds2 <- HU_anova_ls$dds
# gene_names2 <- HU_anova_genes$gene
get_mean_table <- function(dds2, strain, cols2, gene_names2, padj2, cutoff = 0.1) {
  padj_df <- tibble(Gene = gene_names2, P_adj = padj2)
  normalized_counts2 <- counts(dds2, normalized = TRUE)
  normalized_sub2 <- normalized_counts2[which(rownames(normalized_counts2) %in% gene_names2), ] 
  if (is.null(dim(normalized_sub2))) {
    sample_names <- names(normalized_sub2)
    normalized_sub2 <- normalized_sub2 %>% 
      as_tibble %>% 
      mutate(Sample = sample_names) 
    normalized_sub2[[gene_names2]] <- normalized_sub2$value
    normalized_sub2 <- normalized_sub2 %>%
      select(-value) %>% 
      mutate(Sample = colnames(normalized_counts2)) %>% 
      pivot_longer(-Sample, names_to = 'Gene', values_to = 'Normalized_Count') # in the long format for merging
    # normalized_grp2 <- cols2 %>% 
    #   select(Sample, DDG) %>% 
    #   left_join(normalized_sub2, by = 'Sample')
  } else {
    normalized_sub2 <- normalized_sub2 %>% 
      t %>% 
      as_tibble %>% 
      mutate(Sample = colnames(normalized_counts2)) %>% 
      pivot_longer(-Sample, names_to = 'Gene', values_to = 'Normalized_Count') # in the long format for merging
  }
  normalized_grp2 <- cols2 %>% 
    select(Sample, DDG) %>% 
    left_join(normalized_sub2, by = 'Sample') 
  # normalized_grp2 %>% 
  #   filter(Gene == '4933406C10Rik' & DDG == 'FD_0ppb') %>% 
  #   pluck('Normalized_Count') %>% 
  #   mean # check the above results
  mean_normalized_tbl2 <- normalized_grp2 %>% 
    group_by(Gene, DDG) %>% 
    summarise(Mean = mean(Normalized_Count)) %>% 
    ungroup() # long format table
  final_tbl2 <- mean_normalized_tbl2 %>% 
    pivot_wider(Gene, names_from = 'DDG', values_from = 'Mean') %>% 
    left_join(padj_df, by = 'Gene') %>% 
    arrange(P_adj) %>% 
    filter(P_adj < cutoff)
  final_tbl2
}

## function for drawing a Venn diagram from a gene list
get_venn <- function(gene_lst) {
  vd <- ggvenn(
    gene_lst, 
    fill_color = c("#0073C2FF", "#EFC000FF"),
    stroke_size = 0.5, set_name_size = 4
  )
  vd
}

## function for getting Venn diagrams for more than two groups
get_venn2 <- function(sig_ls, show_percentage = FALSE) {
  name_ls <- lapply(sig_ls[-seq(length(sig_ls) - 3, length(sig_ls))], function(x) x$gene)
  vd <- ggvenn(
    name_ls, 
    show_percentage = show_percentage,
    stroke_size = 0.5,
    set_name_size = 4
  )
  vd
}

## get long formatted normalized counts from the normalized count matrix and the coldata
get_long_format <- function(nc, cols) {
  nc_df <- as_tibble(nc, rownames = 'Gene') %>% 
    pivot_longer(-Gene, names_to = 'Sample', values_to = 'Normalized_count') %>% 
    left_join(cols %>% select(Sample, Diet, Dose), by = 'Sample')
  nc_df
}

## get a list of genes from the resulting analysis
# sorted_res <- sperm1
get_gene_list <- function(sorted_res, p_cutoff = 1) {
  require(dplyr)
  sorted_sub <- sorted_res %>% 
    filter(padj <= p_cutoff)
  ls <- sorted_sub$log2FoldChange
  names(ls) <- sorted_sub$geneSymbol
  ls
}


## function for performing the pathway analysis
my_GSEA <- function(gene_list, GO_file, pval, min_size = 0, max_size = Inf, collapse = TRUE, seed = 100, lfc_cutoff = 1) {
  set.seed(seed)
  require(dplyr)
  require(data.table)
  require(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names; remove the duplicates")
    gene_list <- gene_list[!duplicated(names(gene_list))]
  }
  if( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted; sort the gene list based on decreasing order")
    gene_list <- sort(gene_list, decreasing = TRUE)
  }
  
  ## filter gene list based on the LFC cutoff
  gene_list <- gene_list[abs(gene_list) >= lfc_cutoff]
  
  myGO <- fgsea::gmtPathways(GO_file)
  names(myGO) <- substring(names(myGO), 10)
  fgRes <- fgsea::fgsea(pathways = myGO,
                        stats = gene_list,
                        minSize = min_size, ## minimum gene set size
                        maxSize = max_size ## maximum gene set size
                        # , nperm = 10000
                        ) %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval) %>% 
    arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))
  
  ## can choose if the pathways need to be collapsed
  if (collapse) {
    message("Collapsing Pathways -----")
    concise_pathways <- collapsePathways(data.table::as.data.table(fgRes),
                                         pathways = myGO,
                                         stats = gene_list)
    fgRes <- fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
    message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  }
  
  fgRes$Enrichment <- ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes <- rbind(head(fgRes, n = 10), tail(fgRes, n = 10 ))
  total_up <- sum(fgRes$Enrichment == "Up-regulated")
  total_down <- sum(fgRes$Enrichment == "Down-regulated")
  header <- paste0("Top 10 (Total pathways: Up=", total_up,", Down=",    total_down, ")")
  
  colos <- setNames(c("firebrick2", "dodgerblue2"),
                    c("Up-regulated", "Down-regulated"))
  
  ## control the pathway labels to be no longer than 15 characters per row
  filtRes$pathway2 <- stringr::str_wrap(filtRes$pathway, width = 15, whitespace_only = FALSE)
  g1 <- ggplot(filtRes, aes(reorder(pathway2, NES) , NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) +
    scale_fill_manual(values = colos ) +
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score (NES)",
         title=header) + 
    theme_minimal() +
    theme(aspect.ratio = 1)
  
  output <- list("Results" = fgRes, "Plot" = g1)
  return(output)
}

# rotatedAxisElementText = function(angle,position='x'){
#   angle     = angle[1]; 
#   position  = position[1]
#   positions = list(x=0,y=90,top=180,right=270)
#   if(!position %in% names(positions))
#     stop(sprintf("'position' must be one of [%s]",paste(names(positions),collapse=", ")),call.=FALSE)
#   if(!is.numeric(angle))
#     stop("'angle' must be numeric",call.=FALSE)
#   rads  = (angle - positions[[ position ]])*pi/180
#   hjust = 0.5*(1 - sin(rads))
#   vjust = 0.5*(1 + cos(rads))
#   element_text(angle=angle,vjust=vjust,hjust=hjust)
# }


