#!/usr/bin/env Rscript
# Validate environment setup before launching app

cat("=== Clinical Variant Explorer Environment Check ===\n\n")

errors <- c()
warnings <- c()

# Check for VCF file
cat("Checking for VCF file in /data...\n")
vcf_files <- list.files("/data", pattern = "\\.vcf(\\.gz)?$", 
                        recursive = TRUE, full.names = TRUE)
if (length(vcf_files) == 0) {
  errors <- c(errors, "No VCF file found in /data. Please attach a VCF Data Asset.")
} else {
  cat("  ✓ Found VCF:", vcf_files[1], "\n")
  if (length(vcf_files) > 1) {
    warnings <- c(warnings, paste("Multiple VCF files found; will use:", vcf_files[1]))
  }
}

# Check for SnpEff database
cat("\nChecking for SnpEff database at /data/snpeff/...\n")
if (!file.exists("/data/snpeff/snpEff.jar")) {
  errors <- c(errors, "SnpEff.jar not found at /data/snpeff/snpEff.jar. Ensure snpeff-grch38-105 Data Asset is mounted to 'snpeff'.")
} else {
  cat("  ✓ Found SnpEff jar\n")
}

if (!file.exists("/data/snpeff/data/GRCh38.105/snpEffectPredictor.bin")) {
  errors <- c(errors, "GRCh38.105 database not found. Ensure snpeff-grch38-105 Data Asset is complete.")
} else {
  cat("  ✓ Found GRCh38.105 database\n")
}

# Check for FASTA (optional but recommended)
cat("\nChecking for reference FASTA (optional)...\n")
fasta_files <- list.files("/data", pattern = "\\.fa(\\.gz)?$|.fasta(\\.gz)?$", 
                          recursive = TRUE, full.names = TRUE)
if (length(fasta_files) == 0) {
  warnings <- c(warnings, "No reference FASTA found. IGV viewer may not display sequences.")
} else {
  cat("  ✓ Found FASTA:", fasta_files[1], "\n")
}

# Check Java
cat("\nChecking Java installation...\n")
java_version <- tryCatch({
  system("java -version 2>&1 | head -n 1", intern = TRUE)
}, error = function(e) NULL)

if (is.null(java_version)) {
  errors <- c(errors, "Java not found. SnpEff requires Java.")
} else {
  cat("  ✓ Java installed:", java_version, "\n")
}

# Summary
cat("\n=== Summary ===\n")
if (length(errors) > 0) {
  cat("\n❌ ERRORS (must fix before running):\n")
  for (e in errors) cat("  -", e, "\n")
}

if (length(warnings) > 0) {
  cat("\n⚠️  WARNINGS (recommended to fix):\n")
  for (w in warnings) cat("  -", w, "\n")
}

if (length(errors) == 0 && length(warnings) == 0) {
  cat("\n✓ All checks passed! You're ready to launch the app.\n")
}

if (length(errors) > 0) {
  cat("\n❌ Please fix errors before proceeding.\n")
  quit(status = 1)
} else {
  cat("\n✓ Environment validation complete.\n")
}
