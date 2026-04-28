suppressPackageStartupMessages({
  library(shiny)
  library(shinyFiles)
  library(shinyjs)
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(DT)
  library(plotly)
  library(qqman)
  library(SNPRelate)
  library(karyoploteR)
})

# Source helpers
source("/code/helpers/preprocess.R")
source("/code/helpers/export_utils.R")
source("/code/helpers/pca_utils.R")
source("/code/helpers/igv_utils.R")

server <- function(input, output, session) {
  
  # File browser for /data/
  volumes <- c(data = "/data")
  shinyFileChoose(input, "vcf_browser", roots = volumes, 
                  filetypes = c("vcf", "vcf.gz"))
  
  # Selected VCF path
  selected_vcf <- reactiveVal(NULL)
  
  observe({
    req(input$vcf_browser)
    if (is.integer(input$vcf_browser)) return()
    
    file_selected <- parseFilePaths(volumes, input$vcf_browser)
    if (nrow(file_selected) > 0) {
      selected_vcf(as.character(file_selected$datapath[1]))
    }
  })
  
  # Display selected path
  output$selected_vcf_path <- renderText({
    vcf <- selected_vcf()
    if (is.null(vcf)) {
      "No file selected"
    } else {
      paste("Selected:", vcf)
    }
  })
  
  # Enable/disable Process button
  observe({
    if (is.null(selected_vcf())) {
      shinyjs::disable("process_vcf")
    } else {
      shinyjs::enable("process_vcf")
    }
  })
  
  # Data loaded flag
  data_loaded <- reactiveVal(FALSE)
  
  # Check if data exists at startup
  observe({
    if (file.exists("/scratch/variants.parquet")) {
      data_loaded(TRUE)
    }
  })
  
  # Disable filters until data is loaded
  observe({
    if (!data_loaded()) {
      shinyjs::disable("chr_filter")
      shinyjs::disable("start_pos")
      shinyjs::disable("end_pos")
      shinyjs::disable("var_type")
      shinyjs::disable("qual_filter")
      shinyjs::disable("pass_only")
      shinyjs::disable("maf_filter")
      shinyjs::disable("gene_search")
      shinyjs::disable("consequence_filter")
      shinyjs::disable("impact_filter")
      shinyjs::disable("reset_filters")
    } else {
      shinyjs::enable("chr_filter")
      shinyjs::enable("start_pos")
      shinyjs::enable("end_pos")
      shinyjs::enable("var_type")
      shinyjs::enable("qual_filter")
      shinyjs::enable("pass_only")
      shinyjs::enable("maf_filter")
      shinyjs::enable("gene_search")
      shinyjs::enable("consequence_filter")
      shinyjs::enable("impact_filter")
      shinyjs::enable("reset_filters")
    }
  })
  
  # Process VCF button handler
  observeEvent(input$process_vcf, {
    req(selected_vcf())
    
    vcf_path <- selected_vcf()
    force <- input$force_reprocess
    
    withProgress(message = "Processing VCF", value = 0, {
      
      progress_cb <- function(fraction, message, detail = "") {
        setProgress(value = fraction, message = message, detail = detail)
      }
      
      result <- tryCatch({
        preprocess_vcf(
          vcf_path = vcf_path,
          force = force,
          progress_cb = progress_cb
        )
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error", duration = NULL)
        return(NULL)
      })
      
      if (!is.null(result)) {
        data_loaded(TRUE)
        
        if (result$from_cache) {
          showNotification(
            sprintf("Cached results loaded in %.1f seconds", result$elapsed),
            type = "message", duration = 5
          )
        } else {
          showNotification(
            sprintf("Processing complete: %d variants in %.1f seconds", 
                    result$n_variants, result$elapsed),
            type = "message", duration = 8
          )
        }
      }
    })
  })
  
  # Load cached data
  variants_data <- reactive({
    req(data_loaded())
    req(file.exists("/scratch/variants.parquet"))
    read_parquet("/scratch/variants.parquet")
  })
  
  # Load gene list for autocomplete
  observe({
    if (file.exists("/scratch/genes.rds")) {
      genes <- readRDS("/scratch/genes.rds")
      updateSelectizeInput(session, "gene_search", 
                           choices = genes, server = TRUE)
    }
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
    
    # QUAL filter
    df <- df %>% filter(QUAL >= input$qual_filter[1] & QUAL <= input$qual_filter[2])
    
    # PASS filter
    if (input$pass_only) {
      df <- df %>% filter(FILTER == "PASS" | FILTER == ".")
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
    updateSliderInput(session, "qual_filter", value = c(20, 100))
    updateCheckboxInput(session, "pass_only", value = TRUE)
    updateSliderInput(session, "maf_filter", value = c(0.01, 0.5))
    updateSelectizeInput(session, "gene_search", selected = character(0))
    updateSelectInput(session, "consequence_filter", selected = NULL)
    updateCheckboxGroupInput(session, "impact_filter", 
                              selected = c("HIGH", "MODERATE"))
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
  
  # Manhattan plot
  output$manhattan_plot <- renderPlotly({
    df <- filtered_variants()
    req(nrow(df) > 0)
    
    # Prepare data for Manhattan plot
    df_manhattan <- df %>%
      mutate(
        CHR = as.integer(gsub("chr", "", gsub("X", "23", gsub("Y", "24", gsub("M", "25", CHROM))))),
        BP = POS,
        P = 10^(-QUAL / 10)  # Convert QUAL to pseudo p-value
      ) %>%
      filter(!is.na(CHR) & !is.na(BP) & !is.na(P))
    
    plot_ly(df_manhattan, x = ~BP, y = ~-log10(P), color = ~as.factor(CHR),
            type = "scatter", mode = "markers",
            text = ~paste("ID:", ID, "<br>Gene:", Gene_Name, "<br>Effect:", Effect),
            hoverinfo = "text") %>%
      layout(
        xaxis = list(title = "Genomic Position"),
        yaxis = list(title = "-log10(P)"),
        showlegend = FALSE
      )
  })
  
  # MAF histogram
  output$maf_histogram <- renderPlotly({
    if (!data_loaded()) {
      return(plot_ly() %>% 
               layout(title = list(text = "No variant data loaded.", x = 0.5, xanchor = "center")))
    }
    
    df <- filtered_variants()
    req(nrow(df) > 0, "AF" %in% names(df))
    
    df <- df %>% mutate(MAF = pmin(AF, 1 - AF))
    
    plot_ly(df, x = ~MAF, type = "histogram", nbinsx = 50) %>%
      layout(
        xaxis = list(title = "Minor Allele Frequency"),
        yaxis = list(title = "Count")
      )
  })
  
  # Site frequency spectrum
  output$sfs_plot <- renderPlotly({
    df <- filtered_variants()
    req(nrow(df) > 0, "AC" %in% names(df), "AN" %in% names(df))
    
    df_sfs <- df %>%
      mutate(freq_bin = cut(AC / AN, breaks = seq(0, 1, 0.05), include.lowest = TRUE)) %>%
      count(freq_bin) %>%
      filter(!is.na(freq_bin))
    
    plot_ly(df_sfs, x = ~freq_bin, y = ~n, type = "bar") %>%
      layout(
        xaxis = list(title = "Allele Frequency Bin"),
        yaxis = list(title = "Count")
      )
  })
  
  # PCA computation (reactive)
  pca_result <- reactive({
    if (!data_loaded()) return(NULL)
    
    df <- filtered_variants()
    if (nrow(df) == 0) return(NULL)
    
    # Check sample count
    if (!"n_samples" %in% names(df) || is.na(df$n_samples[1]) || df$n_samples[1] < 3) {
      return(NULL)
    }
    
    # Debounce: only recompute if filter changed > 500ms ago
    # For now, compute on demand
    tryCatch({
      compute_pca(
        vcf_path = "/scratch/annotated.vcf",
        filtered_variant_ids = df$ID,
        ld_prune = input$ld_prune,
        progress_cb = function(msg) message(msg)
      )
    }, error = function(e) {
      message("PCA computation failed: ", e$message)
      return(NULL)
    })
  })
  
  # PCA plot
  output$pca_plot <- renderPlotly({
    if (!data_loaded()) {
      return(plot_ly() %>% 
               layout(title = list(text = "No variant data loaded. Use the Data panel in the sidebar to select a VCF file from /data/, then click Process VCF.", 
                                   x = 0.5, xanchor = "center")))
    }
    
    result <- pca_result()
    
    if (is.null(result)) {
      df <- variants_data()
      n_samples <- if ("n_samples" %in% names(df)) df$n_samples[1] else 0
      return(plot_ly() %>% 
               layout(title = list(text = sprintf("PCA requires at least 3 samples — current data has %d sample(s)", n_samples),
                                   x = 0.5, xanchor = "center")))
    }
    
    plot_pca(result, pc_x = 1, pc_y = 2)
  })
  
  # Karyotype plot
  output$karyotype_plot <- renderPlot({
    if (!data_loaded()) {
      plot.new()
      text(0.5, 0.5, "No variant data loaded. Use the Data panel in the sidebar to select a VCF file from /data/, then click Process VCF.", 
           cex = 1.2, col = "gray40")
      return()
    }
    
    df <- filtered_variants()
    req(nrow(df) > 0)
    
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
    if (!data_loaded()) {
      return(p("No variant data loaded. Use the Data panel in the sidebar to select a VCF file from /data/, then click Process VCF.",
               style = "color: gray; font-size: 14px;"))
    }
    
    # Initialize IGV with reference and VCF tracks
    ref_fasta <- list.files("/data/reference", pattern = "\\.fa$|\\.fasta$", 
                            full.names = TRUE, recursive = FALSE)[1]
    ref_fai <- paste0(ref_fasta, ".fai")
    
    if (is.na(ref_fasta) || !file.exists(ref_fasta) || !file.exists(ref_fai)) {
      return(p("Reference genome or index not found in /data/reference/. IGV requires FASTA and .fai index files.",
               style = "color: red;"))
    }
    
    vcf_tracks <- c(selected_vcf(), "/scratch/annotated.vcf")
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
      
      # Export plots
      if (input$export_plots) {
        # Manhattan
        manhattan_path <- file.path("/results", paste0("manhattan_", timestamp, ".png"))
        export_plot_png(output$manhattan_plot(), manhattan_path)
        exported_files <- c(exported_files, manhattan_path)
        
        # MAF histogram
        maf_path <- file.path("/results", paste0("maf_histogram_", timestamp, ".png"))
        export_plot_png(output$maf_histogram(), maf_path)
        exported_files <- c(exported_files, maf_path)
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
  })
}
   removeModal()
    showNotification("Export complete! Check /results", type = "message")
  })
}
ts", type = "message")
  })
}
