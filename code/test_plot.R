#!/usr/bin/env Rscript
# Standalone test script - generates a simple plot to validate environment

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

cat("=== Clinical Variant Explorer Test ===\n")
cat("Generating test plot...\n")

# Create synthetic variant data
set.seed(123)
test_data <- data.frame(
  CHROM = rep(paste0("chr", 1:22), each = 50),
  POS = sample(1:248956422, 1100, replace = TRUE),
  QUAL = runif(1100, 20, 100),
  AF = rbeta(1100, 2, 8),
  Gene_Name = sample(c("BRCA1", "TP53", "EGFR", "KRAS", "MYC"), 1100, replace = TRUE),
  Impact = sample(c("HIGH", "MODERATE", "LOW"), 1100, replace = TRUE, prob = c(0.1, 0.4, 0.5))
)

# Generate Manhattan-style plot
p <- ggplot(test_data, aes(x = POS, y = -log10(1 - QUAL/100), color = CHROM)) +
  geom_point(alpha = 0.6, size = 2) +
  facet_wrap(~CHROM, ncol = 6, scales = "free_x") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Test Variant Distribution Across Chromosomes",
    subtitle = "Clinical Variant Explorer - Environment Validation",
    x = "Genomic Position",
    y = "-log10(P-value proxy)"
  )

# Save plot
output_file <- "/results/test_manhattan_plot.png"
ggsave(output_file, p, width = 14, height = 10, dpi = 150)

# Verify file size
file_info <- file.info(output_file)
file_size_kb <- round(file_info$size / 1024, 2)

cat("\n✓ Plot generated successfully!\n")
cat("  File:", output_file, "\n")
cat("  Size:", file_size_kb, "KB\n")

if (file_size_kb > 0) {
  cat("\n✓ SUCCESS: Plot is >0KB\n")
  quit(status = 0)
} else {
  cat("\n✗ FAIL: Plot is 0KB\n")
  quit(status = 1)
}
