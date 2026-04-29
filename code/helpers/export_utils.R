# Export utilities for Clinical Variant Explorer
# Handles VCF, CSV, PNG, and PDF report exports

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

#' Export filtered variants to VCF
#'
#' @param filtered_df Data frame of filtered variants
#' @param output_path Output VCF file path
#' @param original_vcf Path to original annotated VCF
export_vcf <- function(filtered_df, output_path, original_vcf = "/scratch/annotated.vcf") {
  
  if (!file.exists(original_vcf)) {
    stop("Annotated VCF not found at ", original_vcf)
  }
  
  # Write variant IDs to temp file
  temp_ids <- tempfile()
  writeLines(filtered_df$ID, temp_ids)
  
  # Use bcftools to extract matching variants
  cmd <- sprintf(
    "bcftools view -i 'ID=@%s' %s | bgzip > %s",
    shQuote(temp_ids), shQuote(original_vcf), shQuote(output_path)
  )
  
  system(cmd)
  
  # Create tabix index
  system(sprintf("tabix -p vcf %s", shQuote(output_path)))
  
  unlink(temp_ids)
  
  if (file.exists(output_path)) {
    return(output_path)
  } else {
    stop("VCF export failed")
  }
}

#' Export filtered variants to CSV
#'
#' @param filtered_df Data frame of filtered variants
#' @param output_path Output CSV file path
export_csv <- function(filtered_df, output_path) {
  arrow::write_csv_arrow(filtered_df, output_path)
  
  if (file.exists(output_path)) {
    return(output_path)
  } else {
    stop("CSV export failed")
  }
}

#' Export plot to PNG
#'
#' @param plot_obj Plotly or ggplot object
#' @param output_path Output PNG file path
#' @param width Width in pixels
#' @param height Height in pixels
export_plot_png <- function(plot_obj, output_path, width = 1200, height = 800) {
  if (!requireNamespace("htmlwidgets", quietly = TRUE) ||
      !requireNamespace("webshot2", quietly = TRUE)) {
    stop("Plot export requires htmlwidgets and webshot2 packages")
  }
  
  temp_html <- tempfile(fileext = ".html")
  on.exit(unlink(temp_html), add = TRUE)
  
  if (inherits(plot_obj, "plotly")) {
    htmlwidgets::saveWidget(plot_obj, temp_html, selfcontained = TRUE)
  } else if (inherits(plot_obj, "ggplot")) {
    htmlwidgets::saveWidget(plotly::ggplotly(plot_obj), temp_html, selfcontained = TRUE)
  } else {
    stop("Unsupported plot object type: ", paste(class(plot_obj), collapse = ", "))
  }
  
  webshot2::webshot(temp_html, output_path, vwidth = width, vheight = height)
  
  if (!file.exists(output_path)) {
    stop("Plot export failed: webshot2 did not produce ", output_path)
  }
  output_path
}

#' Export PDF report
#'
#' @param filtered_df Data frame of filtered variants
#' @param filter_params List of active filter parameters
#' @param plots List of plot objects (manhattan, maf, sfs, pca, karyotype)
#' @param output_path Output PDF file path
export_pdf_report <- function(filtered_df, filter_params, plots, output_path) {
  
  template_path <- "/code/report_template.Rmd"
  
  if (!file.exists(template_path)) {
    stop("Report template not found at ", template_path)
  }
  
  # Render with parameters
  rmarkdown::render(
    template_path,
    output_file = output_path,
    params = list(
      filtered_df = filtered_df,
      filter_params = filter_params,
      plots = plots
    ),
    envir = new.env()
  )
  
  if (file.exists(output_path)) {
    return(output_path)
  } else {
    stop("PDF report export failed")
  }
}
