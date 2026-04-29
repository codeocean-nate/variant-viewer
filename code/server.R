suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(arrow)
  library(dplyr)
  library(DT)
  library(plotly)
  library(karyoploteR)
  library(jsonlite)
})

# Source helpers
source("/code/helpers/export_utils.R")
source("/code/helpers/igv_utils.R")
source("/code/helpers/ui_helpers.R")

server <- function(input, output, session) {
  
  # Data-loaded flag
  data_loaded <- reactiveVal(FALSE)
  
  `%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a
  
  available_datasets <- reactive({
    dirs <- list.dirs("/data", recursive = FALSE, full.names = TRUE)
    dirs <- dirs[vapply(dirs, function(d) file.exists(file.path(d, "variants.parquet")), logical(1))]
    if (!length(dirs)) return(character(0))
    setNames(dirs, basename(dirs))
  })
  
  observe({
    ch <- available_datasets()
    updateSelectInput(session, "dataset_picker", choices = ch,
                      selected = if (length(ch)) ch[[1]] else character(0))
  })
  
  selected_dataset <- reactive({
    d <- input$dataset_picker
    if (is.null(d) || d == "" || !dir.exists(d)) NULL else d
  })
  
  observe({
    d <- selected_dataset()
    data_loaded(!is.null(d) && file.exists(file.path(d, "variants.parquet")))
  })
  
  output$dataset_meta <- renderText({
    d <- selected_dataset()
    if (is.null(d)) return("No annotated-variants dataset found in /data/. Attach one.")
    mp <- file.path(d, "vcf_meta.json")
    if (!file.exists(mp)) return(basename(d))
    m <- jsonlite::read_json(mp)
    sprintf("%s • %s variants • %s samples", m$source_vcf %||% "?",
            m$n_variants %||% "?", m$n_samples %||% "?")
  })
  
  # Disable the entire filter panel until data is loaded
  observe({
    shinyjs::toggleState("filter_panel", condition = data_loaded())
  })
  
  # Load variants data from selected dataset
  variants_data <- reactive({
    req(data_loaded())
    d <- selected_dataset()
    req(d)
    tryCatch({
      read_parquet(file.path(d, "variants.parquet"))
    }, error = function(e) {
      showNotification(paste("Error loading variants:", e$message), type = "error", duration = NULL)
      NULL
    })
  })
  
  # Load gene list for autocomplete
  observe({
    d <- selected_dataset(); req(d)
    gp <- file.path(d, "genes.rds")
    if (file.exists(gp))
      updateSelectizeInput(session, "gene_search", choices = readRDS(gp), server = TRUE)
  })
  
  # Update consequence filter choices dynamically
  observe({
    df <- variants_data()
    if ("Effect" %in% names(df)) {
      consequences <- unique(na.omit(df$Effect))
      updateSelectInput(session, "consequence_filter", 
                        choices = consequences)
    }
  })
  
  # Filtered data
  filtered_variants <- reactive({
    df <- variants_data()
    
    # Chromosome filter
    if (!is.null(input$chr_filter) && !"All" %in% input$chr_filter) {
      df <- df %>% filter(CHROM %in% input$chr_filter)
    }
    
    # Position filter (if single chromosome selected)
    if (length(input$chr_filter) == 1 && input$chr_filter != "All") {
      if (!is.null(input$start_pos) && !is.na(input$start_pos)) {
        df <- df %>% filter(POS >= input$start_pos)
      }
      if (!is.null(input$end_pos) && !is.na(input$end_pos)) {
        df <- df %>% filter(POS <= input$end_pos)
      }
    }
    
    # Variant type filter (simplified heuristic)
    if (!is.null(input$var_type) && length(input$var_type) > 0) {
      df <- df %>%
        mutate(
          var_type = case_when(
            nchar(REF) == 1 & nchar(ALT) == 1 ~ "SNP",
            nchar(REF) != nchar(ALT) ~ "Indel",
            TRUE ~ "MNP"
          )
        ) %>%
        filter(var_type %in% input$var_type)
    }
    
    # QUAL filter (handle NA values based on checkbox)
    if ("QUAL" %in% names(df)) {
      if (input$include_na_qual) {
        df <- df %>% filter(is.na(QUAL) | (QUAL >= input$qual_filter[1] & QUAL <= input$qual_filter[2]))
      } else {
        df <- df %>% filter(!is.na(QUAL) & QUAL >= input$qual_filter[1] & QUAL <= input$qual_filter[2])
      }
    }
    
    # PASS filter
    if (input$pass_only) {
      df <- df %>% filter(FILTER == "PASS" | FILTER == ".")
    }
    
    # Allele Frequency filter
    if ("AF" %in% names(df)) {
      df <- df %>% filter(AF >= input$af_filter[1] & AF <= input$af_filter[2])
    }
    
    # MAF filter
    if ("AF" %in% names(df)) {
      df <- df %>%
        mutate(MAF = pmin(AF, 1 - AF)) %>%
        filter(MAF >= input$maf_filter[1] & MAF <= input$maf_filter[2])
    }
    
    # Gene filter
    if (!is.null(input$gene_search) && length(input$gene_search) > 0) {
      df <- df %>% filter(Gene_Name %in% input$gene_search)
    }
    
    # Consequence filter
    if (!is.null(input$consequence_filter) && length(input$consequence_filter) > 0) {
      df <- df %>% filter(Effect %in% input$consequence_filter)
    }
    
    # Impact filter
    if (!is.null(input$impact_filter) && length(input$impact_filter) > 0) {
      df <- df %>% filter(Impact %in% input$impact_filter)
    }
    
    df
  })
  
  # Reset filters
  observeEvent(input$reset_filters, {
    updateSelectInput(session, "chr_filter", selected = "All")
    updateCheckboxGroupInput(session, "var_type", 
                              selected = c("SNP", "Indel", "MNP"))
    updateSliderInput(session, "qual_filter", value = c(20, 1000))
    updateCheckboxInput(session, "include_na_qual", value = TRUE)
    updateCheckboxInput(session, "pass_only", value = TRUE)
    updateSliderInput(session, "af_filter", value = c(0, 1))
    updateSliderInput(session, "maf_filter", value = c(0, 0.5))
    updateSelectizeInput(session, "gene_search", selected = character(0))
    updateSelectInput(session, "consequence_filter", selected = NULL)
    updateCheckboxGroupInput(session, "impact_filter", 
                              selected = c("HIGH", "MODERATE", "LOW", "MODIFIER"))
  })
  
  # Table output
  output$variants_table <- renderDT({
    df <- filtered_variants()
    
    datatable(
      df,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        order = list(list(5, 'desc'))  # Sort by QUAL descending
      ),
      selection = 'single',
      rownames = FALSE
    )
  })
  
  # Track selected variant for IGV
  selected_variant <- reactive({
    s <- input$variants_table_rows_selected
    if (!is.null(s) && length(s) > 0) {
      filtered_variants()[s, ]
    } else {
      NULL
    }
  })
  

  
  # Reactive plot builder — used by both render output and export handler
  build_maf_histogram <- reactive({
    df <- filtered_variants()
    if (nrow(df) == 0 || !"AF" %in% names(df)) return(NULL)
    df <- df %>% mutate(MAF = pmin(AF, 1 - AF))
    plot_ly(df, x = ~MAF, type = "histogram", nbinsx = 50) %>%
      layout(xaxis = list(title = "Minor Allele Frequency"),
             yaxis = list(title = "Count"))
  })
  
  # MAF histogram
  output$maf_histogram <- renderPlotly({
    if (!data_loaded()) return(empty_state_plotly())
    p <- build_maf_histogram()
    if (is.null(p)) return(empty_state_plotly("AF field not available or no variants match filters."))
    p
  })
  

  

  # Karyotype plot
  output$karyotype_plot <- renderPlot({
    if (!data_loaded()) {
      empty_state_baseplot()
      return()
    }
    df <- filtered_variants()
    if (nrow(df) == 0) {
      empty_state_baseplot("No variants match current filters.")
      return()
    }
    
    # Convert to GRanges
    gr <- GRanges(
      seqnames = df$CHROM,
      ranges = IRanges(start = df$POS, width = 1),
      impact = df$Impact
    )
    
    # Plot karyotype with variant density
    kp <- plotKaryotype(genome = "hg38", plot.type = 2)
    kpPlotDensity(kp, data = gr, window.size = 1e6)
    
  }, height = 700)
  
  # IGV viewer
  output$igv_ui <- renderUI({
    if (!data_loaded()) return(empty_state_html())
    
    # Initialize IGV with reference and VCF tracks
    ref_fasta <- list.files("/data/reference", pattern = "\\.fa$|\\.fasta$", 
                            full.names = TRUE, recursive = FALSE)[1]
    ref_fai <- paste0(ref_fasta, ".fai")
    
    if (is.na(ref_fasta) || !file.exists(ref_fasta) || !file.exists(ref_fai)) {
      return(p("Reference genome or index not found in /data/reference/. IGV requires FASTA and .fai index files.",
               style = "color: red;"))
    }
    
    d <- selected_dataset()
    vcf_tracks <- if (!is.null(d)) file.path(d, "annotated.vcf.gz") else character(0)
    vcf_tracks <- vcf_tracks[file.exists(vcf_tracks)]
    
    tryCatch({
      init_igv(session, ref_fasta, ref_fai, vcf_tracks)
    }, error = function(e) {
      p(paste("IGV initialization failed:", e$message), style = "color: red;")
    })
  })
  
  # Navigate IGV on table row click
  observe({
    var <- selected_variant()
    if (!is.null(var)) {
      navigate_igv(session, chr = var$CHROM, pos = var$POS, window = 100)
    }
  })
  
  # Export handlers
  observeEvent(input$export_btn, {
    showModal(modalDialog(
      title = "Export Options",
      checkboxInput("export_vcf", "Filtered VCF", value = TRUE),
      checkboxInput("export_csv", "Table CSV", value = TRUE),
      checkboxInput("export_plots", "Plots (PNG)", value = FALSE),
      checkboxInput("export_pdf", "PDF Report", value = FALSE),
      footer = tagList(
        actionButton("do_export", "Export"),
        modalButton("Cancel")
      )
    ))
  })
  
  observeEvent(input$do_export, {
    req(data_loaded())
    
    df <- filtered_variants()
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    exported_files <- c()
    
    tryCatch({
      
      # Export VCF
      if (input$export_vcf) {
        vcf_path <- file.path("/results", paste0("filtered_", timestamp, ".vcf.gz"))
        export_vcf(df, vcf_path)
        exported_files <- c(exported_files, vcf_path)
      }
      
      # Export CSV
      if (input$export_csv) {
        csv_path <- file.path("/results", paste0("filtered_variants_", timestamp, ".csv"))
        export_csv(df, csv_path)
        exported_files <- c(exported_files, csv_path)
      }
      
      # Export plots — use the build_* reactive, NOT the output$ render expression
      if (input$export_plots) {
        maf_plot <- build_maf_histogram()
        if (!is.null(maf_plot)) {
          maf_path <- file.path("/results", paste0("maf_histogram_", timestamp, ".png"))
          export_plot_png(maf_plot, maf_path)
          exported_files <- c(exported_files, maf_path)
        }
      }
      
      # Export PDF report
      if (input$export_pdf) {
        pdf_path <- file.path("/results", paste0("variant_report_", timestamp, ".pdf"))
        
        filter_params <- list(
          chr = input$chr_filter,
          qual_range = input$qual_filter,
          pass_only = input$pass_only,
          maf_range = input$maf_filter,
          genes = input$gene_search,
          consequences = input$consequence_filter,
          impacts = input$impact_filter
        )
        
        # Collect plots (placeholder for now - would need actual plot objects)
        plots <- list()
        
        export_pdf_report(df, filter_params, plots, pdf_path)
        exported_files <- c(exported_files, pdf_path)
      }
      
      removeModal()
      showNotification(
        sprintf("Export complete! %d file(s) saved to /results", length(exported_files)),
        type = "message", duration = 8
      )
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Export failed:", e$message), type = "error", duration = NULL)
    })
  })
}
