# PCA module helpers using SNPRelate
# Handles VCF → GDS conversion, LD pruning, and PCA calculation

suppressPackageStartupMessages({
  library(SNPRelate)
  library(plotly)
})

#' Compute PCA from filtered VCF
#'
#' @param vcf_path Path to VCF file (or use annotated VCF from /scratch)
#' @param filtered_variant_ids Vector of variant IDs to include
#' @param ld_prune Logical, perform LD pruning
#' @param sample_metadata_path Optional path to sample metadata TSV
#' @param progress_cb Optional progress callback function
#' @return List with pca_result, eigenvalues, sample_ids, and metadata
compute_pca <- function(vcf_path = "/scratch/annotated.vcf",
                       filtered_variant_ids = NULL,
                       ld_prune = TRUE,
                       sample_metadata_path = "/data/sample_metadata.tsv",
                       progress_cb = NULL) {
  
  report_progress <- function(message) {
    if (!is.null(progress_cb)) progress_cb(message)
  }
  
  # Convert VCF to GDS
  report_progress("Converting VCF to GDS...")
  gds_file <- "/scratch/pca.gds"
  
  tryCatch({
    SNPRelate::snpgdsVCF2GDS(vcf_path, gds_file, method = "biallelic.only", verbose = FALSE)
  }, error = function(e) {
    stop("Failed to convert VCF to GDS: ", e$message)
  })
  
  gdsobj <- SNPRelate::snpgdsOpen(gds_file)
  on.exit(SNPRelate::snpgdsClose(gdsobj), add = TRUE)
  
  # Filter to selected variant IDs if provided
  snp_ids <- NULL
  if (!is.null(filtered_variant_ids) && length(filtered_variant_ids) > 0) {
    all_snp_ids <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
    snp_ids <- all_snp_ids[all_snp_ids %in% filtered_variant_ids]
    if (length(snp_ids) == 0) snp_ids <- NULL
  }
  
  # LD pruning
  if (ld_prune) {
    report_progress("LD pruning...")
    set.seed(1000)
    snpset_pruned <- SNPRelate::snpgdsLDpruning(
      gdsobj,
      ld.threshold = sqrt(0.2),
      slide.max.bp = 500000,
      snp.id = snp_ids,
      verbose = FALSE
    )
    snp_ids <- unlist(snpset_pruned, use.names = FALSE)
  }
  
  # Run PCA
  report_progress("Running PCA...")
  pca <- SNPRelate::snpgdsPCA(gdsobj, snp.id = snp_ids, num.thread = 2, verbose = FALSE)
  
  # Load sample metadata if available
  sample_info <- NULL
  if (file.exists(sample_metadata_path)) {
    sample_info <- read.table(sample_metadata_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }
  
  list(
    pca = pca,
    eigenvalues = pca$eigenval,
    sample_ids = pca$sample.id,
    sample_info = sample_info
  )
}

#' Plot PCA results with plotly
#'
#' @param pca_result Result from compute_pca
#' @param pc_x PC for X axis (default 1)
#' @param pc_y PC for Y axis (default 2)
#' @return Plotly object
plot_pca <- function(pca_result, pc_x = 1, pc_y = 2) {
  
  pca <- pca_result$pca
  sample_info <- pca_result$sample_info
  
  # Build data frame
  df <- data.frame(
    sample_id = pca$sample.id,
    pc1 = pca$eigenvect[, 1],
    pc2 = pca$eigenvect[, 2],
    pc3 = if (ncol(pca$eigenvect) >= 3) pca$eigenvect[, 3] else NA,
    pc4 = if (ncol(pca$eigenvect) >= 4) pca$eigenvect[, 4] else NA,
    pc5 = if (ncol(pca$eigenvect) >= 5) pca$eigenvect[, 5] else NA,
    pc6 = if (ncol(pca$eigenvect) >= 6) pca$eigenvect[, 6] else NA,
    pc7 = if (ncol(pca$eigenvect) >= 7) pca$eigenvect[, 7] else NA,
    pc8 = if (ncol(pca$eigenvect) >= 8) pca$eigenvect[, 8] else NA,
    pc9 = if (ncol(pca$eigenvect) >= 9) pca$eigenvect[, 9] else NA,
    pc10 = if (ncol(pca$eigenvect) >= 10) pca$eigenvect[, 10] else NA,
    stringsAsFactors = FALSE
  )
  
  # Merge with sample metadata for coloring
  color_col <- NULL
  if (!is.null(sample_info) && "group" %in% names(sample_info)) {
    df <- merge(df, sample_info[, c("sample_id", "group"), drop = FALSE], by = "sample_id", all.x = TRUE)
    color_col <- ~group
  }
  
  # Select PCs for axes
  x_col <- paste0("pc", pc_x)
  y_col <- paste0("pc", pc_y)
  
  if (!x_col %in% names(df) || !y_col %in% names(df)) {
    return(plot_ly() %>% layout(title = "Selected PCs not available"))
  }
  
  # Variance explained
  var_explained <- pca$varprop * 100
  x_lab <- sprintf("PC%d (%.2f%%)", pc_x, var_explained[pc_x])
  y_lab <- sprintf("PC%d (%.2f%%)", pc_y, var_explained[pc_y])
  
  # Plot
  p <- plot_ly(df, x = as.formula(paste0("~", x_col)), y = as.formula(paste0("~", y_col)),
               type = "scatter", mode = "markers",
               text = ~sample_id, hoverinfo = "text",
               marker = list(size = 8))
  
  if (!is.null(color_col)) {
    p <- plot_ly(df, x = as.formula(paste0("~", x_col)), y = as.formula(paste0("~", y_col)),
                 color = color_col,
                 type = "scatter", mode = "markers",
                 text = ~sample_id, hoverinfo = "text",
                 marker = list(size = 8))
  }
  
  p %>% layout(
    xaxis = list(title = x_lab),
    yaxis = list(title = y_lab)
  )
}
