# Refactored preprocess.R as callable function for Shiny server
# Converts VCF to annotated Parquet with progress callbacks

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(digest)
})

#' Preprocess VCF file with SnpEff annotation
#'
#' @param vcf_path Path to input VCF file
#' @param snpeff_data_dir Path to SnpEff data directory
#' @param snpeff_jar Path to SnpEff JAR file
#' @param reference_dir Path to reference genome directory
#' @param cache_dir Directory for cached outputs
#' @param force Force re-processing even if cache is valid
#' @param progress_cb Progress callback function(fraction, message, detail)
#' @return List with parquet_path, genes_rds_path, n_variants, from_cache
preprocess_vcf <- function(vcf_path,
                           snpeff_data_dir = "/data/snpeff/data",
                           snpeff_jar = "/data/snpeff/snpEff.jar",
                           reference_dir = "/data/reference",
                           cache_dir = "/scratch",
                           force = FALSE,
                           progress_cb = NULL) {
  
  start_time <- Sys.time()
  
  # Helper to call progress callback safely
  report_progress <- function(fraction, message, detail = "") {
    if (!is.null(progress_cb)) {
      tryCatch(progress_cb(fraction, message, detail), error = function(e) NULL)
    }
  }
  
  # File paths
  parquet_file <- file.path(cache_dir, "variants.parquet")
  sha_file <- file.path(cache_dir, "source.sha256")
  genes_file <- file.path(cache_dir, "genes.rds")
  annotated_vcf <- file.path(cache_dir, "annotated.vcf")
  
  # Compute SHA256
  compute_sha <- function(file) {
    digest(file, algo = "sha256", file = TRUE)
  }
  
  # Check cache validity
  cache_valid <- function() {
    if (force) return(FALSE)
    if (!file.exists(parquet_file) || !file.exists(sha_file)) return(FALSE)
    cached_sha <- readLines(sha_file, warn = FALSE)[1]
    current_sha <- compute_sha(vcf_path)
    identical(cached_sha, current_sha)
  }
  
  # Cache hit path
  if (cache_valid()) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    n_variants <- nrow(read_parquet(parquet_file))
    return(list(
      parquet_path = parquet_file,
      genes_rds_path = genes_file,
      n_variants = n_variants,
      from_cache = TRUE,
      elapsed = elapsed
    ))
  }
  
  # Full processing path
  
  # Stage 1: Count input variants (5%)
  report_progress(0, "Reading VCF header and counting records...")
  
  # Quick count via bcftools if available
  n_total <- tryCatch({
    bcftools_out <- system2("bcftools", args = c("view", "-H", shQuote(vcf_path), "|", "wc", "-l"),
                            stdout = TRUE, stderr = FALSE)
    as.integer(bcftools_out)
  }, error = function(e) {
    # Fallback: parse VCF header
    con <- file(vcf_path, "r")
    count <- 0
    while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
      if (!startsWith(line, "#")) count <- count + 1
    }
    close(con)
    count
  })
  
  report_progress(0.05, sprintf("Found %d variants", n_total), "")
  
  # Stage 2: Run SnpEff (70%)
  report_progress(0.05, "Annotating with SnpEff...", sprintf("0 / %d variants", n_total))
  
  if (!file.exists(snpeff_jar)) {
    stop("SnpEff JAR not found at ", snpeff_jar)
  }
  
  # Run SnpEff synchronously
  exit_code <- system2(
    "java",
    args = c(
      "-Xmx8g",
      "-jar", snpeff_jar,
      "-dataDir", snpeff_data_dir,
      "GRCh38.105",
      vcf_path
    ),
    stdout = annotated_vcf,
    stderr = FALSE
  )
  
  if (exit_code != 0 || !file.exists(annotated_vcf) || file.size(annotated_vcf) == 0) {
    stop("SnpEff annotation failed. Check that /data/snpeff contains GRCh38.105 database.")
  }
  
  report_progress(0.75, "SnpEff complete", sprintf("%d / %d variants", n_total, n_total))
  
  # Stage 3: Parse annotated VCF (15%)
  report_progress(0.75, "Parsing annotated VCF...")
  
  ref_genome <- list.files(reference_dir, pattern = "\\.fa$|\\.fasta$",
                           full.names = TRUE, recursive = FALSE)[1]
  
  if (is.na(ref_genome) || !file.exists(ref_genome)) {
    stop("Reference genome FASTA not found in ", reference_dir,
         ". Attach the hg38 reference Data Asset to /data/reference/ before processing.")
  }
  if (!file.exists(paste0(ref_genome, ".fai"))) {
    stop("Reference index (.fai) not found alongside ", ref_genome,
         ". The reference Data Asset must include a .fai index.")
  }
  
  vcf <- readVcf(annotated_vcf, genome = ref_genome)
  
  # Extract fields
  df <- data.frame(
    CHROM = as.character(seqnames(vcf)),
    POS = start(vcf),
    ID = names(vcf),
    REF = as.character(ref(vcf)),
    ALT = sapply(alt(vcf), function(x) paste(as.character(x), collapse = ",")),
    QUAL = qual(vcf),
    FILTER = sapply(filt(vcf), paste, collapse = ","),
    stringsAsFactors = FALSE
  )
  
  # INFO fields
  info_fields <- info(vcf)
  if ("AF" %in% names(info_fields)) df$AF <- info_fields$AF
  if ("AC" %in% names(info_fields)) df$AC <- info_fields$AC
  if ("AN" %in% names(info_fields)) df$AN <- info_fields$AN
  
  # SnpEff ANN field — pick the most-severe annotation per variant
  if ("ANN" %in% names(info_fields)) {
    ann <- info_fields$ANN
    impact_rank <- c(HIGH = 4, MODERATE = 3, LOW = 2, MODIFIER = 1)
    
    pick_top_annotation <- function(ann_field) {
      if (is.null(ann_field) || is.na(ann_field) || ann_field == "") {
        return(rep(NA_character_, 7))
      }
      annotations <- unlist(strsplit(as.character(ann_field), ","))
      parts_list <- lapply(annotations, function(a) unlist(strsplit(a, "\\|")))
      impacts <- sapply(parts_list, function(p) if (length(p) >= 3) p[3] else NA)
      ranks <- impact_rank[impacts]
      ranks[is.na(ranks)] <- 0
      top <- parts_list[[which.max(ranks)]]
      get_field <- function(x, i) if (length(x) >= i) x[i] else NA_character_
      c(
        Effect        = get_field(top, 2),
        Impact        = get_field(top, 3),
        Gene_Name     = get_field(top, 4),
        Gene_ID       = get_field(top, 5),
        Transcript_ID = get_field(top, 7),
        HGVS_c        = get_field(top, 10),
        HGVS_p        = get_field(top, 11)
      )
    }
    
    parsed <- t(sapply(as.character(ann), pick_top_annotation))
    df$Effect        <- parsed[, "Effect"]
    df$Impact        <- parsed[, "Impact"]
    df$Gene_Name     <- parsed[, "Gene_Name"]
    df$Gene_ID       <- parsed[, "Gene_ID"]
    df$Transcript_ID <- parsed[, "Transcript_ID"]
    df$HGVS_c        <- parsed[, "HGVS_c"]
    df$HGVS_p        <- parsed[, "HGVS_p"]
  }
  
  # Genotype fields
  geno_list <- geno(vcf)
  if ("DP" %in% names(geno_list)) df$mean_DP <- rowMeans(geno_list$DP, na.rm = TRUE)
  if ("GQ" %in% names(geno_list)) df$mean_GQ <- rowMeans(geno_list$GQ, na.rm = TRUE)
  
  report_progress(0.90, "Parsing complete", "")
  
  # Stage 4: Write Parquet cache (5%)
  report_progress(0.90, "Caching results...")
  write_parquet(df, parquet_file)
  writeLines(compute_sha(vcf_path), sha_file)
  
  # Stage 5: Build gene index (5%)
  report_progress(0.95, "Building gene autocomplete index...")
  if ("Gene_Name" %in% names(df)) {
    genes <- unique(na.omit(df$Gene_Name))
    genes <- sort(genes)
    saveRDS(genes, genes_file)
  }
  
  report_progress(1.0, "Complete", "")
  
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  list(
    parquet_path = parquet_file,
    genes_rds_path = genes_file,
    n_variants = nrow(df),
    from_cache = FALSE,
    elapsed = elapsed
  )
}
