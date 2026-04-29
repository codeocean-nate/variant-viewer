suppressPackageStartupMessages({
  library(shinydashboard)
  library(shinyjs)
  library(DT)
  library(plotly)
})

ui <- dashboardPage(
  skin = "blue",
  
  dashboardHeader(
    title = "Variant Viewer",
    tags$li(
      class = "dropdown",
      actionButton("export_btn", "Export", icon = icon("download"))
    )
  ),
  
  dashboardSidebar(
    useShinyjs(),
    
    tags$div(
      id = "data_panel",
      style = "padding: 15px; border-bottom: 1px solid #ddd;",
      h4("Variants Dataset", style = "margin-top: 0;"),
      selectInput("dataset_picker", NULL, choices = NULL),
      tags$small(textOutput("dataset_meta"), style = "color: #666;")
    ),
    hr(),
    
    tags$div(
      id = "filter_panel",
      
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
                  min = 0, max = 1000, value = c(20, 1000), step = 10),
      
      checkboxInput("pass_only", "PASS Variants Only", value = TRUE),
      
      sliderInput("maf_filter", "Minor Allele Frequency",
                  min = 0, max = 0.5, value = c(0.01, 0.5), step = 0.01),
      
      selectizeInput("gene_search", "Gene Name (HGNC)",
                     choices = NULL, multiple = TRUE,
                     options = list(placeholder = "Search genes...")),
      
      selectInput("consequence_filter", "Consequence",
                  choices = NULL, multiple = TRUE),
      
      checkboxGroupInput("impact_filter", "SnpEff Impact",
                         choices = c("HIGH", "MODERATE", "LOW", "MODIFIER"),
                         selected = c("HIGH", "MODERATE", "LOW", "MODIFIER")),
      
      hr(),
      actionButton("reset_filters", "Reset Filters", icon = icon("refresh"))
    )
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
      
      tabPanel("Allele Frequency", icon = icon("chart-bar"),
               fluidRow(
                 column(12,
                        h3("Minor Allele Frequency Distribution"),
                        plotlyOutput("maf_histogram", height = "500px")
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
