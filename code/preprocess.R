#!/usr/bin/env Rscript
# Stage 1: SnpEff annotation + VCF parsing → cached Parquet

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(digest)
})

CACHE_DIR <- "/scratch"
PARQUET_FILE <- file.path(CACHE_DIR, "variants.parquet")
SHA_FILE <- file.path(CACHE_DIR, "source.sha256")
GENES_FILE <- file.path(CACHE_DIR, "genes.rds")

# Discover VCF in /data
find_vcf <- function() {
  vcf_files <- list.files("/data", pattern = "\\.vcf(\\.gz)?$", 
                          recursive = TRUE, full.names = TRUE)
  if (length(vcf_files) == 0) {
    stop("No VCF file found in /data")
  }
  if (length(vcf_files) > 1) {
    warning("Multiple VCF files found; using first: ", vcf_files[1])
  }
  vcf_files[1]
}

# Compute SHA256 of VCF
compute_sha <- function(file) {
  digest(file, algo = "sha256", file = TRUE)
}

# Check if cache is valid
cache_valid <- function(vcf_path) {
  if (!file.exists(PARQUET_FILE) || !file.exists(SHA_FILE)) {
    return(FALSE)
  }
  cached_sha <- readLines(SHA_FILE, warn = FALSE)[1]
  current_sha <- compute_sha(vcf_path)
  identical(cached_sha, current_sha)
}

# Run SnpEff annotation
run_snpeff <- function(vcf_path, output_path) {
  cat("Running SnpEff annotation...\n")
  snpeff_jar <- "/data/snpeff/snpEff.jar"
  snpeff_data <- "/data/snpeff/data"
  
  if (!file.exists(snpeff_jar)) {
    stop("SnpEff jar not found at ", snpeff_jar, 
         "\nEnsure snpeff-grch38-105 Data Asset is mounted at /data/snpeff/")
  }
  
  cmd <- sprintf(
    "java -Xmx8g -jar %s -dataDir %s GRCh38.105 %s > %s",
    shQuote(snpeff_jar), shQuote(snpeff_data), 
    shQuote(vcf_path), shQuote(output_path)
  )
  
  system(cmd, intern = FALSE)
  
  if (!file.exists(output_path)) {
    stop("SnpEff annotation failed; output VCF not created")
  }
}

# Parse annotated VCF
parse_vcf <- function(vcf_path) {
  cat("Parsing VCF...\n")
  # Find reference genome (should be in /data/reference/)
  ref_genome <- list.files("/data/reference", pattern = "\\.fa$|\\.fasta$", 
                           full.names = TRUE, recursive = FALSE)[1]
  if (is.na(ref_genome) || !file.exists(ref_genome)) {
    warning("Reference genome not found in /data/reference/. Using default GRCh38.")
    vcf <- readVcf(vcf_path, genome = "GRCh38")
  } else {
    cat("Using reference:", ref_genome, "\n")
    vcf <- readVcf(vcf_path, genome = ref_genome)
  }
  
  # Extract fixed fields
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
  
  # Extract INFO fields
  info_fields <- info(vcf)
  if ("AF" %in% names(info_fields)) df$AF <- info_fields$AF
  if ("AC" %in% names(info_fields)) df$AC <- info_fields$AC
  if ("AN" %in% names(info_fields)) df$AN <- info_fields$AN
  
  # Parse SnpEff ANN field (simplified - first annotation per variant)
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
  
  # Extract genotype fields (GT, DP, GQ) - simplified
  geno_list <- geno(vcf)
  if ("GT" %in% names(geno_list)) {
    df$n_samples <- ncol(geno_list$GT)
  }
  if ("DP" %in% names(geno_list)) {
    df$mean_DP <- rowMeans(geno_list$DP, na.rm = TRUE)
  }
  if ("GQ" %in% names(geno_list)) {
    df$mean_GQ <- rowMeans(geno_list$GQ, na.rm = TRUE)
  }
  
  df
}

# Main preprocessing
main <- function() {
  cat("=== Clinical Variant Explorer Preprocessing ===\n")
  
  vcf_path <- find_vcf()
  cat("Found VCF:", vcf_path, "\n")
  
  if (cache_valid(vcf_path)) {
    cat("Cache is valid; skipping preprocessing.\n")
    return(invisible(NULL))
  }
  
  cat("Cache invalid or missing; running full preprocessing...\n")
  
  # Annotate with SnpEff
  annotated_vcf <- file.path(CACHE_DIR, "annotated.vcf")
  run_snpeff(vcf_path, annotated_vcf)
  
  # Parse VCF
  variants <- parse_vcf(annotated_vcf)
  
  # Write to Parquet
  cat("Writing to Parquet...\n")
  write_parquet(variants, PARQUET_FILE)
  
  # Save SHA
  writeLines(compute_sha(vcf_path), SHA_FILE)
  
  # Build gene autocomplete index
  if ("Gene_Name" %in% names(variants)) {
    genes <- unique(na.omit(variants$Gene_Name))
    genes <- sort(genes)
    saveRDS(genes, GENES_FILE)
    cat("Saved", length(genes), "unique gene names for autocomplete.\n")
  }
  
  cat("Preprocessing complete. Variants:", nrow(variants), "\n")
}

main()
