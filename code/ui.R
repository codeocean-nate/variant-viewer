suppressPackageStartupMessages({
  library(shinydashboard)
  library(shinyjs)
  library(DT)
  library(plotly)
})

ui <- dashboardPage(
  skin = "blue",
  
  dashboardHeader(
    title = "Clinical Variant Explorer",
    tags$li(
      class = "dropdown",
      actionButton("export_btn", "Export", icon = icon("download"))
    )
  ),
  
  dashboardSidebar(
    useShinyjs(),
    
    # Data panel (collapsible)
    tags$div(
      id = "data_panel",
      style = "padding: 15px; border-bottom: 1px solid #ddd;",
      h4("Data", style = "margin-top: 0;"),
      shinyFilesButton("vcf_browser", "Browse VCF Files", 
                       "Select a VCF file from /data/", 
                       multiple = FALSE, icon = icon("folder-open")),
      br(), br(),
      textOutput("selected_vcf_path"),
      br(),
      checkboxInput("force_reprocess", "Force re-process (ignore cache)", value = FALSE),
      br(),
      actionButton("process_vcf", "Process VCF", icon = icon("play"), 
                   class = "btn-primary", width = "100%")
    ),
    
    hr(),
    
    h4("Filters", style = "padding-left: 15px;"),
    
    selectInput("chr_filter", "Chromosome",
                choices = c("All", paste0("chr", c(1:22, "X", "Y", "M"))),
                selected = "All", multiple = TRUE),
    
    conditionalPanel(
      condition = "input.chr_filter.length == 1 && input.chr_filter != 'All'",
      numericInput("start_pos", "Start Position", value = NULL),
      numericInput("end_pos", "End Position", value = NULL)
    ),
    
    checkboxGroupInput("var_type", "Variant Type",
                       choices = c("SNP", "Indel", "MNP"),
                       selected = c("SNP", "Indel", "MNP")),
    
    sliderInput("qual_filter", "QUAL Score",
                min = 0, max = 100, value = c(20, 100)),
    
    checkboxInput("pass_only", "PASS Variants Only", value = TRUE),
    
    sliderInput("maf_filter", "Minor Allele Frequency",
                min = 0, max = 0.5, value = c(0.01, 0.5), step = 0.01),
    
    selectizeInput("gene_search", "Gene Name",
                   choices = NULL, multiple = TRUE,
                   options = list(placeholder = "Search genes...")),
    
    selectInput("consequence_filter", "Consequence",
                choices = NULL, multiple = TRUE),
    
    checkboxGroupInput("impact_filter", "SnpEff Impact",
                       choices = c("HIGH", "MODERATE", "LOW", "MODIFIER"),
                       selected = c("HIGH", "MODERATE")),
    
    hr(),
    actionButton("reset_filters", "Reset Filters", icon = icon("refresh"))
  ),
  
  dashboardBody(
    useShinyjs(),
    
    tabsetPanel(
      id = "main_tabs",
      
      tabPanel("Table", icon = icon("table"),
               fluidRow(
                 column(12,
                        h3("Filtered Variants"),
                        DTOutput("variants_table")
                 )
               )
      ),
      
      tabPanel("Manhattan Plot", icon = icon("chart-area"),
               fluidRow(
                 column(12,
                        h3("Manhattan Plot"),
                        plotlyOutput("manhattan_plot", height = "600px")
                 )
               )
      ),
      
      tabPanel("Allele Frequency", icon = icon("chart-bar"),
               fluidRow(
                 column(6,
                        h4("MAF Distribution"),
                        plotlyOutput("maf_histogram", height = "400px")
                 ),
                 column(6,
                        h4("Site Frequency Spectrum"),
                        plotlyOutput("sfs_plot", height = "400px")
                 )
               )
      ),
      
      tabPanel("PCA", icon = icon("project-diagram"),
               fluidRow(
                 column(12,
                        h3("Principal Component Analysis"),
                        checkboxInput("ld_prune", "LD Pruning (r² < 0.2)", value = TRUE),
                        plotlyOutput("pca_plot", height = "600px")
                 )
               )
      ),
      
      tabPanel("Karyotype", icon = icon("dna"),
               fluidRow(
                 column(12,
                        h3("Variant Density by Chromosome"),
                        plotOutput("karyotype_plot", height = "700px")
                 )
               )
      ),
      
      tabPanel("IGV", icon = icon("eye"),
               fluidRow(
                 column(12,
                        h3("IGV Genome Browser"),
                        p("Click a variant in the Table tab to navigate here."),
                        uiOutput("igv_ui")
                 )
               )
      )
    )
  )
)
