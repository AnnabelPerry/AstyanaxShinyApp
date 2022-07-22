# Annabel Perry
# July 2022
# When it was created in Summer 2022, this script was originally intended as a 
# fill-in-the-blank style exercise for Zelun Liu & Daniel Nguyen to learn how to
# code the QTL module user interface in CaveCrawler.

# TODO: Throughout this script, there are comments marked with "TODO". These 
# comments describe portions of code which need to be added (aka "pseudocode").
# After the code described in a "TODO" comment has been successfully added and
# troubleshooted, please erase the "TODO" part of the comment, but don't erase 
# the whole comment! Remember that you will have to describe this code to a NEW 
# undergrad at some point (maybe 1-3 years from now!). So, rephrase the comment
# in a manner which will help describe the code to an unfamiliar user. If a 
# comment does NOT start with "TODO", there is no need to mess with it.

source("functions.R")

ui <- fluidPage(
  tabsetPanel(
    sidebarLayout(
      sidebarPanel(
        # With this radiobutton, the user specifies if they would like to...
        # 1. Use the Marker Range sub-module.
        # 2. Use the Genomic Range sub-module.
        # 3. Use the Trait-to-Marker sub-module.
        radioButtons(
          inputId = "QTLsub_mod",
          label = "I want to...",
          selected = character(0),
          # The choicenames are user-friendly descriptions of the actions of 
          # different sub-modules
          choiceNames = c("... search for markers and genes within range of a 
                          known marker or gene.", 
                          "... search for markers and genes in a specific region
                          of the genome.",
                          "... search for markers associated with a quantitative
                          trait"),
          # The choice values are abbreviations for the sub-modules which the 
          # user-friendly names refer to.
          choiceValues = c("MR", "GR", "TM")
        ),
        conditionalPanel(
          condition = "input.QTLsub_mod == 'MR'",
          # TODO MR:
          #   Based on the slides describing the "Marker Range" sub-module and 
          #   the 'MR' inputs specified in the QTL function, add widgets for the 
          #   Marker Range sub-module.
          textInput("marker_term",label = "Search for a marker:",width = "500px"),
          textInput("marker_position",label = "Search for position around marker:",width = "500px"),
          actionButton(
            # Unique identifier for the action button
            inputId = "MR_action_button",
            # Label to be displayed on the action button
            label = "Find Genes & Markers"),
          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                           tags$div("Crawling through the data...",id="loadmessage"))
        ),
        conditionalPanel(
          condition = "input.QTLsub_mod == 'GR'",
          # TODO GR:
          #   Based on the slides describing the "Genomic Range" sub-module and 
          #   the 'GR' inputs specified in the QTL function, add widgets for the 
          #   Genomic Range sub-module.
          # Hint: You may need to load additional libraries.
          # Hint: Use your example apps!
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
                selected = NULL
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

            
          #),
          # TODO TM:
          #   Based on the slides describing the "Trait-to-Marker" sub-module and 
          #   the 'TM' inputs specified in the QTL function, add widgets for the 
          #   Trait-to-Marker sub-module.
          # Hint: You may need to load additional libraries.
          # Hint: Use your example apps!
        ),
        "Chromosome lengths obtained from assembly GCA_000372685.2 via the 
        European Nucleotide Archive Browser (Accessed July 2022)"
      ),
      mainPanel(
        # Note: Later, Annabel will add plot output here
        conditionalPanel(
          condition = "input.GR_action_button | input.MR_action_button",
          tableOutput("QTL_Marker_table"),
          tableOutput("QTL_Gene_table")
        ),
        conditionalPanel(
          condition = "input.TM_enter",
          tableOutput("QTLmarker_table")
        ),
        textOutput("chr_plot"),
        tableOutput("QTL_wrnings")
      )
    )
  )
)

server <- function(input, output, session) {
  
  QTLoutput <- eventReactive(input$QTLsub_mod, valueExpr = {
    # TODO P: To functions.R, read in the CSVs needed to make this work
    QTL(chr_table, position_table, QTL_table,
        # TODO GR: To minimize errors in simultaneously troubleshooting multiple
        #          sub-modules, these arguments currently are set to not have
        #          values. Fill in the arguments with the reactive inputs from
        #          the widgets you made in the UI section.
        # Hint: Control-F for "GeneSearch(" in the official CaveCrawler app.R. 
        # Hint: Use your example apps!
        # GR.bool = F, GR.chr = NA, GR.start = NA, GR.end = NA,
        if (input$QTLsub_mod == 'GR') {
          GR.bool = T
        },
        if (input$GR_chr != "") {
          GR.chr = as.numeric(input$GR_chr)
        },
        if (input$GR_start != "") {
          GR.start = as.numeric(input$GR_start)
        },
        if (input$GR_end != "") {
          GR.end = as.numeric(input$GR_end)
        },
        # TODO MR: To minimize errors in simultaneously troubleshooting multiple
        #          sub-modules, these arguments currently are set to not have
        #          values. Fill in the arguments with the reactive inputs from
        #          the widgets you made in the UI section.
        # Hint: Control-F for "GeneSearch(" in the official CaveCrawler app.R. 
        # Hint: Use your example apps!
        # MR.bool = F, MR.search_term = NA, MR.bp = NA,
        if (input$QTLsub_mod == 'MR') {
          MR.bool = T
        } else {
          MR.bool = F
        },
        if (length(input$marker_term) > 0) {
          MR.search_term = input$marker_term
        },
        if (length(input$marker_position) > 0) {
          MR.bp = as.numeric(input$marker_position)
        },
        # TODO TM: To minimize errors in simultaneously troubleshooting multiple
        #          sub-modules, these arguments currently are set to not have
        #          values. Fill in the arguments with the reactive inputs from
        #          the widgets you made in the UI section.
        # Hint: Control-F for "GeneSearch(" in the official CaveCrawler app.R. 
        # Hint: Use your example apps!
        if(input$QTLsub_mod == 'TM'){
          TM.bool = T
        },
        if(input$Trait_search > 0){
          TM.QT = input$Trait_search
        }
    )
  })
  
  output$QTL_Marker_table <- renderTable({QTLoutput()[[2]]})
  output$QTL_Gene_table <- renderTable({QTLoutput()[[3]]})
  # Note: Later, Annabel will add plot output here
  output$QTLwrnings <- renderText({QTLoutput()[[4]]})
  
  output$QTLmarker_table <- renderTable({
    if(('TM' %in% input$QTLsub_mod) & (length(input$Trait_search) > 0)){
      QTLoutput()[[1]]
    } else{
      data.frame()
    }
  })
  
  output$TM_checkboxes <- renderUI({
    #populate checkboxes with all traits in trait column
    traits <- QTL_table$Quantitative_Trait
    #remove duplicates
    traits <- traits[!duplicated(traits) &!grepl('NA', traits) &!is.na(traits)]
    
    checkboxGroupInput(
      inputId = "Trait_search",
      label = "Select trait(s)",
      choices = traits
    )
    
  })
  # Note: Later, Annabel will add plot output here
  
}

shinyApp(ui, server)