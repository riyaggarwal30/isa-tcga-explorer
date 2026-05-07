
# app.R
library(shiny)
library(ggplot2)
library(dplyr)
library(IsoformSwitchAnalyzeR)

# Source our module
source("R/ui_modules.R")

# 1. Load the pre-computed data using a relative path
# (Ensure you ran the scripts/01_data_prep.R first!)
switchData <- readRDS("data/switchData.rds")

# 2. Define Server Logic
server <- function(input, output, session) {
  
  output$switchSummary <- renderTable({
    extractSwitchSummary(switchData, filterForConsequences = TRUE)
  })
  
  output$volcanoPlot <- renderPlot({
    ggplot(data = switchData$isoformFeatures,
           aes(x = dIF, y = -log10(isoform_switch_q_value))) +
      geom_point(aes(color = abs(dIF) > 0.1 & isoform_switch_q_value < 0.05), size = 1.5) +
      scale_color_manual('Significant', values = c('grey40', '#e05c5c')) +
      theme_bw() +
      facet_wrap(~ condition_2)
  })
  
  output$switchPlotTitle <- renderText({
    paste("Isoform Switch:", input$gene, "in", input$comparison)
  })
  
  output$switchPlot <- renderPlot({
    switchPlot(switchData, gene = input$gene, 
               condition1 = paste0(input$comparison, "_ctrl"),
               condition2 = paste0(input$comparison, "_cancer"))
  })
  
  output$consequenceSummary <- renderPlot({
    extractConsequenceSummary(switchData)
  })
  
  output$consequenceEnrichment <- renderPlot({
    extractConsequenceEnrichment(switchData)
  })
  
  output$splicingEnrichment <- renderPlot({
    extractSplicingEnrichment(switchData)
  })
}

# 3. Launch App
shinyApp(ui = create_ui(switchData), server = server)

