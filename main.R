library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('testthat')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
  
  counts_df <- read_tsv(counts_csv)
  meta_df <- read_csv(metafile_csv)
  
  meta_subset <- meta_df %>%
    dplyr::filter(timepoint %in% selected_times) %>%
    dplyr::select(samplename, timepoint) %>%
    dplyr::mutate(timepoint = factor(timepoint, levels = selected_times))
  
  sample_order <- meta_subset$samplename
  
  counts_subset <- counts_df %>%
    dplyr::select(gene, dplyr::all_of(sample_order))
  
  counts_matrix <- as.matrix(counts_subset[, -1])
  rownames(counts_matrix) <- counts_subset$gene
  
  meta_subset <- meta_subset %>%
    dplyr::slice(match(colnames(counts_matrix), samplename))
  
  rownames(meta_subset) <- meta_subset$samplename
  
  se <- SummarizedExperiment(
    assays = list(counts = counts_matrix),
    colData = meta_subset
  )
  
    return(se)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  dds <- DESeq2::DESeqDataSet(se, design = design)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds)
  res_df <- as.data.frame(res)
  
  return(list(dds = dds, results = res_df))
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`. Ensure
#' that the column name for your rownames is called "genes". 
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
  
  labeled_res <- deseq2_res %>%
    tibble::rownames_to_column(var = "genes") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      volc_plot_status = dplyr::case_when(
        padj < padj_threshold & log2FoldChange > 0 ~ "UP",
        padj < padj_threshold & log2FoldChange < 0  ~ "DOWN",
        TRUE ~ "NS"
      )
    )
    
  return(labeled_res)
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  
  my_plot <- labeled_results %>%
    dplyr::filter(!is.na(pvalue)) %>%
    ggplot2::ggplot(ggplot2::aes(x = pvalue)) +
    ggplot2::geom_histogram(bins = 50) +
    ggplot2::labs(
      title = "Raw p-values histogram",
      x = "p-value",
      y = "Frequency"
    ) +
    ggplot2::theme_minimal()
  
  return(my_plot)
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  
  my_plot <- labeled_results %>%
    dplyr::filter(
      !is.na(log2FoldChange),
      !is.na(padj),
      padj < padj_threshold
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = log2FoldChange)) +
    ggplot2::geom_histogram(bins = 50) +
    ggplot2::labs(
      title = paste0("log2 Fold Changes (padj < ", padj_threshold, ") histogram"),
      x = "log2 Fold Change",
      y = "Frequency"
    ) +
    ggplot2::theme_minimal()
    
  return(my_plot)
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  
  # Get top genes by smallest padj:
  
  top_genes <- labeled_results %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice(1:num_genes) %>%
    dplyr::pull(genes)
  
  # Get normalized counts:
  
  normalized_counts <- counts(dds_obj, normalized = TRUE)
  
  # Convert to long format: (Because ggplot wants one row per point it needs to draw on the plot)
  
  df_long <- as.data.frame(normalized_counts) %>%
    tibble::rownames_to_column("genes") %>%
    dplyr::filter(genes %in% top_genes) %>%
    tidyr::pivot_longer(
      cols = -genes,
      names_to = "samplename",
      values_to = "normalized_counts"
    )
  
  # Add metadata (timepoint) using left join:

  meta_df <- as.data.frame(colData(dds_obj))
  
  df_long <- df_long %>%
    dplyr::left_join(meta_df, by = "samplename")
  
  # Plot:
  
  my_plot <- ggplot2::ggplot(df_long, 
                       ggplot2::aes(x = timepoint, 
                                    y = normalized_counts, 
                                    color = timepoint)) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.1)) +
    ggplot2::facet_wrap(~ genes, scales = "free_y") +
    ggplot2::labs(
      title = paste0("Normalized Counts for Top ", num_genes, " Genes"),
      x = "Timepoint",
      y = "Normalized Counts"
    ) +
    ggplot2::theme_minimal()
  
  return(my_plot)
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  
  my_plot <- labeled_results %>%
    dplyr::filter(!is.na(log2FoldChange), !is.na(padj)) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = log2FoldChange,
        y = -log10(padj),
        color = volc_plot_status
      )
    ) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::labs(
      title = "Differential Expression Results Volcano Plot",
      x = "log2 Fold Change",
      y = "-log10 adjusted p-value",
      color = "Status"
    ) +
    ggplot2::theme_minimal()
  
  return(my_plot)
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  
  # NOTE: header = FALSE is required because the id2gene file has no header row.
  # Without it, read.delim treats the first gene (ENSMUSG00000102693.2 -> 4933401J01Rik)
  # as column names instead of data, causing it to return NA on the join.
  
  id2gene <- read.delim(id2gene_path, stringsAsFactors = FALSE, header = FALSE)
  
  # Need to make sure first two columns are genes and symbol:
  
  colnames(id2gene)[1:2] <- c("genes", "symbol")
  
  # Join gene symbol column to our results and clean up NA values:
  
  ranked_df <- labeled_results %>%
    dplyr::left_join(id2gene, by = "genes") %>%
    dplyr::filter(!is.na(log2FoldChange), !is.na(symbol), symbol != "") %>%
    dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
    dplyr::distinct(symbol, .keep_all = TRUE) %>%
    dplyr::arrange(dplyr::desc(log2FoldChange))
  
  # Create a named vector where the values are log2FoldChange numbers and the names are gene symbols:
  
  ranked_vector <- ranked_df$log2FoldChange
  names(ranked_vector) <- ranked_df$symbol
  
    return(ranked_vector)
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  
  pathways <- fgsea::gmtPathways(gmt_file_path)
  
  fgsea_results <- fgsea::fgsea(
    pathways = pathways,
    stats = rnk_list,
    minSize = min_size,
    maxSize = max_size
  )
  
    return(tibble::as_tibble(fgsea_results))
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  
  # Select top positively enriched pathways:
  
  top_pos <- fgsea_results %>%
    dplyr::filter(!is.na(NES)) %>%
    dplyr::arrange(dplyr::desc(NES)) %>%
    dplyr::slice(1:num_paths)
  
  # Select top negatively enriched pathways:
  
  top_neg <- fgsea_results %>%
    dplyr::filter(!is.na(NES)) %>%
    dplyr::arrange(NES) %>%
    dplyr::slice(1:num_paths)
  
  # Combine and label by NES direction:
  
  plot_df <- dplyr::bind_rows(top_pos, top_neg) %>%
    dplyr::mutate(
      nes_sign = dplyr::if_else(NES > 0, "Positive NES", "Negative NES"),
      nes_sign = factor(nes_sign, levels = c("Negative NES", "Positive NES")),
      pathway = factor(pathway, levels = pathway[order(NES)])
    )
  
  # Make bar plot:
  
  my_plot <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = pathway, y = NES, fill = nes_sign)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = c("Negative NES" = "salmon", "Positive NES" = "steelblue")
    ) +
    ggplot2::labs(
      title = paste0("Top ", num_paths, " Positive and Negative Enriched Pathways"),
      x = "Pathway",
      y = "Normalized Enrichment Score (NES)",
      fill = "NES Sign"
    ) +
    ggplot2::theme_minimal()
  
  return(my_plot)
}

