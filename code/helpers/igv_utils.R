# IGV module using igvShiny or iframe fallback
# Provides genome browser functionality with custom reference

#' Initialize IGV viewer
#'
#' @param session Shiny session object
#' @param reference_fasta Path to reference FASTA file
#' @param reference_fai Path to reference FAI index
#' @param vcf_paths Vector of VCF paths to load
#' @return igvShiny widget or fallback HTML
init_igv <- function(session, reference_fasta, reference_fai, vcf_paths = c()) {
  
  # Try igvShiny first
  if (requireNamespace("igvShiny", quietly = TRUE)) {
    tryCatch({
      
      # Configure custom genome
      genome_config <- list(
        id = "custom_genome",
        name = "Custom Reference",
        fastaURL = reference_fasta,
        indexURL = reference_fai
      )
      
      # Initialize igvShiny
      igv_widget <- igvShiny::igvShiny(
        genomeName = "custom",
        genomeOptions = genome_config,
        displayMode = "SQUISHED"
      )
      
      # Load VCF tracks
      for (vcf_path in vcf_paths) {
        if (file.exists(vcf_path)) {
          igvShiny::loadVcfTrack(session, trackName = basename(vcf_path), vcfURL = vcf_path)
        }
      }
      
      return(igv_widget)
      
    }, error = function(e) {
      message("igvShiny failed: ", e$message)
      # Fall through to iframe
    })
  }
  
  # Fallback: iframe with self-hosted igv.js
  # This requires an HTML file serving igv.js
  fallback_iframe()
}

#' Fallback IGV iframe
fallback_iframe <- function() {
  htmltools::tags$iframe(
    src = "https://igv.org/app/",
    width = "100%",
    height = "600px",
    style = "border: 1px solid #ccc;"
  )
}

#' Navigate IGV to locus
#'
#' @param session Shiny session object
#' @param chr Chromosome
#' @param pos Position
#' @param window Window size around position (default 200bp)
navigate_igv <- function(session, chr, pos, window = 200) {
  
  if (requireNamespace("igvShiny", quietly = TRUE)) {
    tryCatch({
      locus <- sprintf("%s:%d-%d", chr, max(1, pos - window), pos + window)
      igvShiny::searchTrackBy(session, locus = locus)
    }, error = function(e) {
      message("IGV navigation failed: ", e$message)
    })
  }
}
