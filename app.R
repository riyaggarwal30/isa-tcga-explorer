# ISA TCGA EXPLORER- An interactive isoform switching analysis dashboard
# Author: Riya Aggarwal
# Data: TCGA example dataset (IsoformSwitchAnalyzeR, Vitting-Seerup et al 2019)



#importing libraries
library(shiny)
library(ggplot2)
library(dplyr)
library(bslib)
library(IsoformSwitchAnalyzeR)


#original data we are using
data("exampleSwitchListAnalyzed")

#saving this data as RDS file on the folder, to avoid loading the entire IsoformSwitchAnalyzeR on server
saveRDS(
  exampleSwitchListAnalyzed, 
  "/Users/riyaaggarwal/Desktop/isa-tcga-explorer/switchData.rds"
)
# loading pre-computed data
exampleSwitchListAnalyzed <- readRDS("switchData.rds")

# Pre-extract top genes for dropdown
topGenes <- exampleSwitchListAnalyzed$isoformFeatures %>%
  filter(abs(dIF) > 0.1 & isoform_switch_q_value < 0.05) %>%
  filter(switchConsequencesGene == TRUE) %>%
  arrange(gene_switch_q_value) %>%
  pull(gene_name) %>%
  unique() %>%
  head(30)


#UI
ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "flatly"),
  
  titlePanel("🧬 ISA TCGA Explorer — Isoform Switch Analysis in Cancer"),
  
  tabsetPanel(
    tabPanel("Overview",
             br(),
             fluidRow(
               column(12,
                      h4("About this Dashboard"),
                      p("This dashboard explores isoform switching in two cancer types from The Cancer Genome Atlas (TCGA):
                      Colon Adenocarcinoma (COAD) and Lung Adenocarcinoma (LUAD), analyzed using the 
                      IsoformSwitchAnalyzeR R package (Vitting-Seerup et al., 2019)."),
                      p("Isoform switching occurs when the dominant transcript of a gene changes between conditions.
                      Unlike standard differential expression analysis, this captures functional changes 
                      that are invisible at the gene level."),
               )
             ),
             br(),
             fluidRow(
               column(6,
                      h4("Switch Summary"),
                      tableOutput("switchSummary")
               ),
               column(6,
                      h4("What is an Isoform Switch?"),
                      p("When total gene expression stays flat but isoform usage shifts dramatically,
                      standard RNA-seq analysis misses it entirely. The plots in this dashboard
                      reveal these hidden regulatory events and their functional consequences."),
                      p("Key metric: dIF (difference in Isoform Fraction) measures the effect size
                      of each switch, equivalent to fold change in standard expression analysis.")
               )
             ),
             br(),
             fluidRow(
               column(12,
                      h4("Volcano Plot — All Isoform Switches"),
                      plotOutput("volcanoPlot", height = "450px")
               )
             )
    ),
    
    tabPanel("Gene Explorer",
             br(),
             fluidRow(
               column(4,
                      selectInput("gene", 
                                  label = "Select Gene:",
                                  choices = extractTopSwitches(
                                    exampleSwitchListAnalyzed,
                                    filterForConsequences = TRUE,
                                    n = 20,
                                    inEachComparison = TRUE
                                  ) %>% pull(gene_name) %>% unique()
                      ),
                      selectInput("comparison",
                                  label = "Select Comparison:",
                                  choices = list(
                                    "COAD (Colon Cancer)" = "COAD",
                                    "LUAD (Lung Cancer)" = "LUAD"
                                  )
                      ),
                      br(),
                      p("The switch plot shows transcript structure (top),
                      gene and isoform expression (bottom left/middle),
                      and isoform fraction — IF (bottom right)."),
                      p("A significant switch occurs when IF changes between
                      conditions while total gene expression may remain stable.")
               ),
               column(8,
                      h4(textOutput("switchPlotTitle")),
                      plotOutput("switchPlot", height = "500px")
               )
             )
    ),
    tabPanel("Consequence Analysis",
             br(),
             fluidRow(
               column(12,
                      h4("Functional Consequences of Isoform Switches"),
                      p("Shows how many isoforms gain or lose specific biological features 
                      as a result of switching. Domain loss is the most frequent consequence 
                      in both cancer types — cancer-associated isoforms often lack functional 
                      motifs present in normal tissue.")
               )
             ),
             br(),
             fluidRow(
               column(12,
                      h4("Consequence Summary"),
                      plotOutput("consequenceSummary", height = "450px")
               )
             ),
             br(),
             fluidRow(
               column(12,
                      h4("Consequence Enrichment — Are Losses More Common Than Gains?"),
                      p("Tests whether losses are systematically more frequent than gains 
                      for each consequence type. A significant bias toward loss suggests 
                      cancer isoforms are functionally truncated versions of normal isoforms."),
                      plotOutput("consequenceEnrichment", height = "450px")
               )
             )
    ),
    tabPanel("Splicing Analysis",
             br(),
             fluidRow(
               column(12,
                      h4("Splicing Mechanisms Underlying Isoform Switches"),
                      p("Identifies which alternative splicing events drive the switches detected above.
                      Each dot represents a splicing type — dot size indicates number of genes affected,
                      position indicates whether gains or losses dominate."),
                      br(),
                      p(strong("Key splicing types:")),
                      tags$ul(
                        tags$li(strong("ATSS"), " — Alternative Transcription Start Site: different promoter usage"),
                        tags$li(strong("ATTS"), " — Alternative Transcription Termination Site: different 3' end usage"),
                        tags$li(strong("A3"), " — Alternative 3' acceptor site"),
                        tags$li(strong("A5"), " — Alternative 5' donor site"),
                        tags$li(strong("SE"), " — Skipped exon"),
                        tags$li(strong("MES"), " — Mutually exclusive exons")
                      )
               )
             ),
             br(),
             fluidRow(
               column(12,
                      h4("Splicing Enrichment Analysis"),
                      plotOutput("splicingEnrichment", height = "500px")
               )
             )
    )
  )
)

#server
server <- function(input, output, session) {
  output$switchSummary <- renderTable({
    extractSwitchSummary(
      exampleSwitchListAnalyzed,
      filterForConsequences = TRUE
    )
  })
  
  output$volcanoPlot <- renderPlot({
    ggplot(data = exampleSwitchListAnalyzed$isoformFeatures,
           aes(x = dIF, y = -log10(isoform_switch_q_value))) +
      geom_point(
        aes(color = abs(dIF) > 0.1 & isoform_switch_q_value < 0.05),
        size = 1.5
      ) +
      geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
      geom_vline(xintercept = c(-0.1, 0.1), linetype = 'dashed') +
      facet_wrap(~ condition_2) +
      scale_color_manual('Significant\nIsoform Switch', values = c('grey40', '#e05c5c')) + 
      labs(x = 'dIF', y = '-Log10 ( Isoform Switch Q Value )') +
      theme_dark() +
      theme_bw() +
      theme(strip.background = element_rect(fill = '#f0f0f0'))
  })
  
  output$switchPlotTitle <- renderText({
    paste("Isoform Switch:", input$gene, 
          ifelse(input$comparison == "COAD", 
                 "(COAD_ctrl vs COAD_cancer)", 
                 "(LUAD_ctrl vs LUAD_cancer)"))
  })
  #for gene comparison
  output$switchPlot <- renderPlot({
    cond1 <- paste0(input$comparison, "_ctrl")
    cond2 <- paste0(input$comparison, "_cancer")
    
    switchPlot(
      exampleSwitchListAnalyzed,
      gene = input$gene,
      condition1 = cond1,
      condition2 = cond2,
      localTheme = theme_bw(base_size = 12)
    )
  })
  
  output$consequenceSummary <- renderPlot({
    extractConsequenceSummary(
      exampleSwitchListAnalyzed,
      consequencesToAnalyze = 'all',
      plotGenes = FALSE,
      asFractionTotal = FALSE
    )
  })
  
  output$consequenceEnrichment <- renderPlot({
    extractConsequenceEnrichment(
      exampleSwitchListAnalyzed,
      consequencesToAnalyze = 'all',
      analysisOppositeConsequence = TRUE,
      returnResult = FALSE
    )
  })
  
  output$splicingEnrichment <- renderPlot({
    extractSplicingEnrichment(
      exampleSwitchListAnalyzed,
      returnResult = FALSE
    )
  })
}

shinyApp(ui = ui, server = server)
