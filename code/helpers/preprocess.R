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
  
  # Run SnpEff with stderr capture
  stderr_file <- tempfile()
  cmd <- sprintf(
    "java -Xmx8g -jar %s -dataDir %s GRCh38.105 %s 2> %s > %s",
    shQuote(snpeff_jar), shQuote(snpeff_data_dir),
    shQuote(vcf_path), shQuote(stderr_file), shQuote(annotated_vcf)
  )
  
  # Launch SnpEff in background
  system(cmd, wait = FALSE, intern = FALSE)
  
  # Poll stderr for progress (every 2 seconds)
  snpeff_complete <- FALSE
  processed <- 0
  
  while (!snpeff_complete) {
    Sys.sleep(2)
    
    # Check if SnpEff finished
    if (file.exists(annotated_vcf)) {
      snpeff_complete <- TRUE
      report_progress(0.75, "SnpEff complete", sprintf("%d / %d variants", n_total, n_total))
      break
    }
    
    # Parse stderr for progress
    if (file.exists(stderr_file)) {
      stderr_lines <- readLines(stderr_file, warn = FALSE)
      # SnpEff writes lines like "Done: 12345 variants"
      progress_lines <- grep("Done:|variants processed", stderr_lines, value = TRUE, ignore.case = TRUE)
      if (length(progress_lines) > 0) {
        last_line <- tail(progress_lines, 1)
        numbers <- as.integer(regmatches(last_line, gregexpr("[0-9]+", last_line))[[1]])
        if (length(numbers) > 0) {
          processed <- max(numbers)
          frac <- min(processed / max(n_total, 1), 1.0)
          progress_val <- 0.05 + 0.70 * frac
          report_progress(progress_val, "Annotating with SnpEff...", 
                          sprintf("%d / %d variants", processed, n_total))
        }
      }
    }
    
    # Timeout after 30 minutes
    elapsed_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    if (elapsed_secs > 1800) {
      stop("SnpEff timed out after 30 minutes")
    }
  }
  
  unlink(stderr_file)
  
  # Stage 3: Parse annotated VCF (15%)
  report_progress(0.75, "Parsing annotated VCF...")
  
  ref_genome <- list.files(reference_dir, pattern = "\\.fa$|\\.fasta$",
                           full.names = TRUE, recursive = FALSE)[1]
  
  vcf <- if (!is.na(ref_genome) && file.exists(ref_genome)) {
    readVcf(annotated_vcf, genome = ref_genome)
  } else {
    readVcf(annotated_vcf, genome = "GRCh38")
  }
  
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
  
  # SnpEff ANN field
  if ("ANN" %in% names(info_fields)) {
    ann <- info_fields$ANN
    ann_split <- strsplit(as.character(ann), "\\|")
    df$Effect <- sapply(ann_split, function(x) if(length(x) > 1) x[2] else NA)
    df$Impact <- sapply(ann_split, function(x) if(length(x) > 2) x[3] else NA)
    df$Gene_Name <- sapply(ann_split, function(x) if(length(x) > 3) x[4] else NA)
    df$Gene_ID <- sapply(ann_split, function(x) if(length(x) > 4) x[5] else NA)
    df$Transcript_ID <- sapply(ann_split, function(x) if(length(x) > 6) x[7] else NA)
    df$HGVS_c <- sapply(ann_split, function(x) if(length(x) > 9) x[10] else NA)
    df$HGVS_p <- sapply(ann_split, function(x) if(length(x) > 10) x[11] else NA)
  }
  
  # Genotype fields
  geno_list <- geno(vcf)
  if ("GT" %in% names(geno_list)) df$n_samples <- ncol(geno_list$GT)
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
