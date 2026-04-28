suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyjs)
})

# Load UI and server
source("ui.R")
source("server.R")

# Run app
shinyApp(ui = ui, server = server)
