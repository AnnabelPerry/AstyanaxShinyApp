source("functions.R")

ui <- fluidPage(
  tabsetPanel(
    sidebarLayout(
      sidebarPanel(
        radioButtons(
          inputId = "QTLsub_mod",
          label = "I want to...",
          selected = character(0),
          choiceNames = c("... search for markers and genes within range of a 
                          known marker or gene.", 
                          "... search for markers and genes in a specific region
                          of the genome.",
                          "... search for markers associated with a quantitative
                          trait"),
          choiceValues = c("MR", "GR", "TM")
        ),
        conditionalPanel(
          condition = "input.QTLsub_mod == 'MR'",
          textInput("marker_term",
                    label = "Search for a marker:",
                    width = "500px"),
          textInput("marker_position",
                    label = "Search for position around marker:",
                    width = "500px"),
          actionButton(
            # Unique identifier for the action button
            inputId = "MR_action_button",
            # Label to be displayed on the action button
            label = "Find Genes & Markers"),
          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                           tags$div("Crawling through the data...",
                                    id="loadmessage"))
        ),
        conditionalPanel(
          condition = "input.QTLsub_mod == 'GR'",
          fluidRow(
            # Any data placed within the same "column()" function call will appear
            # in the same column of that row
            column(
              # We must assign each column a width between 1 and 12
              width = 5,
              br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
              tags$b("Select a chromosome: "),
              br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
              tags$b("Starting position: "),
              br(),br(),
              tags$b("Ending position: "),
              br(),br(),
              actionButton("GR_action_button","Find Genes & Markers"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Crawling through the data...",id="GR_loadmessage"))
            ),
            column(
              # We must assign each column a width between 1 and 12
              width = 7,
              radioButtons(
                # A unique identifier which will be stored as an entry in a list called "input"
                inputId = "GR_chr",
                # The text to be displayed above the radiobuttons
                label = NULL,
                choices = 1:25,
                # Tell shiny whether a default value should show as "selected"
                selected = F
              ),
              textInput("GR_start",label = NULL,width = "500px"),
              textInput("GR_end",label = NULL,width = "500px")
            )
          )
        ),
        conditionalPanel(
          condition = "input.QTLsub_mod == 'TM'",
          uiOutput("TM_checkboxes"),
          actionButton("TM_enter", "Find Markers"),
          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                           tags$div("Crawling through the data...",id="loadmessage"))
          
        ),
        "Chromosome lengths obtained from assembly GCA_000372685.2 via the 
        European Nucleotide Archive Browser (Accessed July 2022)"
      ),
      mainPanel(
        tableOutput("QTL_Marker_table"),
        tableOutput("QTL_Gene_table"),
        textOutput("QTL_wrnings"),
        #tableOutput("QTLmarker_table"),
        #tableOutput("QTLgene_table"),
        # TODO Copy-and-paste
        #chromoMapOutput("QTLplot"),
        #textOutput("QTLwrnings")
        textOutput("test")
      )
    )
  )
)

server <- function(input, output, session) {
  
  QTLoutput <- eventReactive(input$QTLsub_mod, valueExpr = {
    QTL(chr_table, position_table, QTL_table,
        GR.bool = F, GR.chr = NA, GR.start = NA, GR.end = NA,
        MR.bool = F, MR.search_term = NA, MR.bp = NA,
        TM.bool = F, TM.QT = NA)
    })
  
  #output$QTLmarker_table <- renderTable({QTLoutput()[[1]]})
  #output$QTLgene_table <- renderTable({QTLoutput()[[2]]})
  # TODO copy-and-paste
  #output$QTLplot <- renderChromoMap({QTLoutput()[[3]]})
  #output$QTLwrnings <- renderText({QTLoutput()[[4]]})
  output$test <- renderText({c(paste("marker_term == empty string?", 
                                     input$marker_term == ""),
                                     paste("marker_position == empty string?", 
                                           input$marker_position == ""))})
}

shinyApp(ui, server)