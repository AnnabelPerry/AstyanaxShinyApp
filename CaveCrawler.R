#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

########################### Load Required Libraries ############################

library(shinyWidgets)
library(shiny)
library(plotly)
library(WVPlots)
library(stringr)
library(tibble)
library(ggplot2)
library(gridExtra)
library(dplyr)

################################## Load Data ###################################

position_table <- read.csv("data/AmexPositionTable.csv", fill = TRUE)

condition_control <- read.csv("data/Morph_Control_TranscData.csv")

# EDIT: This data is fake. We currently lack between-morph comparisons, so this
# CSV file is simply used as a filler to ensure the function works
morph1.morph2 <- read.csv("data/Toy_RioChoyPachon.csv")

GeneToGO <- read.csv("data/AMexGOTerms.csv", fill = T)

GoIDToNames <- read.table("data/GOIDs_and_Names.txt", fill = T, sep = "\t", header = T)
# When you first read in GoIDToNames, some entire lines are, for whatever reason,
# combined into a single GO term cell

# This section of the script...
# 1. Identifies all GO term cells containing information from multiple lines
messed_cells <- GoIDToNames$GO.Term[grepl("\n", GoIDToNames$GO.Term)]

# 2. Splices the cells into strings based on newline characters and adds strings
#    to a vector
newline_strings <- c()
sites.of.errors <- c()
for(c in 1:length(messed_cells)){
  newline_strings <- append(newline_strings, str_split(string = messed_cells[c],
                                                       pattern = "\n"))
  sites.of.errors <- append(sites.of.errors, 
                            which(GoIDToNames$GO.Term == messed_cells[c]))
}


error.replacements <- c()
for(i in 1:length(newline_strings)){
  # 3. Records sites at which an erroneous GO term must be replaced with a
  #    corrected GO term
  error.replacements <- append(error.replacements, newline_strings[[i]][1])
  newline_strings[[i]] <- newline_strings[[i]][-1]
  # 4. Splices the strings of the vector into substrings based on \t and adds 
  #    substring pairs as columns to a temporary dataframe
  temp_df <- data.frame(matrix(ncol = 2, nrow = 1))
  new_rows <- str_split(string = newline_strings[[i]], pattern = "\t")
  for(r in 1:length(new_rows)){
    temp_df <- rbind(temp_df, new_rows[[r]])
  }
  temp_df <- temp_df[-1,]
  names(temp_df) <- names(GoIDToNames)
  # 5. Inserts the temporary data frame into the master data frame
  GoIDToNames <- rbind(GoIDToNames, temp_df)
  next
}
# 6. Replaces the erroneous GO term with a corrected GO term
GoIDToNames$GO.Term[sites.of.errors] <- error.replacements

UpperLower <- read.table("data/GOTermAssociations.txt", fill = T, sep = "\t", header = T)

# Read in both outlier and non-outlier datasets and integrate into a single data
# frame, adding a column describing the publication from which the data came
outlier_table <- read.csv("data/AMexicanus_Outlier_Stats.csv")
outlier_table <- outlier_table[,(names(outlier_table) != "X")]
outlier_table$Publication_Name <- rep("Outlier paper (Provisional)", 
                                      nrow(outlier_table))
chica_table <- read.csv("data/AMexicanus_iScienceS4R1_Stats.csv")
names(chica_table) <- names(outlier_table)
chica_table$Publication_Name <- rep("Chica paper (Provisional)", 
                                    nrow(chica_table))
chica_table <- chica_table[,1:ncol(outlier_table)]
stat_table <- rbind(outlier_table, chica_table)
stat_table <- stat_table[((stat_table$Gene_Name != "") & 
                                      (stat_table$Gene_Name != " ")),]

GO_classes <- read.table("data/GOID_Namespaces.txt", fill = T, sep = "\t", header = T)

# Obtain complete dataframe of all possible genes and corresponding IDs across
# the statistic and transcription data
all.genes_IDs <- data.frame(
  all_genes = c(position_table$Gene_Name, stat_table$Gene_Name,
                condition_control$Gene_name
  ),
  all_IDs = c(position_table$Gene_ID, stat_table$Stable_Gene_ID,
              condition_control$Gene_stable_ID
  )
)
all.genes_IDs <- all.genes_IDs[!duplicated(all.genes_IDs[,2]),]
all.genes_IDs <- all.genes_IDs[!duplicated(all.genes_IDs[,1]),]

# Obtain a complete vector of all GO IDs
all.GO_IDs <- c(GO_classes$GO_ID)
all.GO_IDs <- all.GO_IDs[!duplicated(all.GO_IDs)]

############################### Source Functions ###############################
source("functions/CaveCrawler_functions.R")

  ui = fluidPage(
    chooseSliderSkin("Flat", color = "#e8c4c2"),
    theme = "style.css",
    tags$head(tags$link(rel="shortcut icon", href="favicon.ico")),
    tags$style(type="text/css",
               ".shiny-output-error { visibility: hidden; }",
               ".shiny-output-error:before { visibility: hidden; }"),
    # Change background color of tabs
    tags$style(HTML("
    .tabbable > .nav > li > a                  {background-color: black; border-color: black;}
    .tabbable > .nav > li[class=active]    > a {background-color: #e8c4c2; border-color: #e4867e}
    .tabbable > .nav > li > a:hover {background-color: #e8c4c2; border-color: #e4867e}
  ")),
    tabsetPanel(
      tabPanel(h1("Home"), fluid = TRUE,
               h1("Welcome to CaveCrawler"),
               br(),
               textOutput("home_text2"),
               br(),
               textOutput("home_text3"),
               br(),
               textOutput("home_text4"),
               br(),
               textOutput("home_text5"),
               br(),
               plotOutput("home_plot")
      ),
      tabPanel(h2("Gene Search"), fluid = TRUE,
               sidebarLayout(
                 sidebarPanel(id = "sidebar",
                   searchInput(
                     inputId = "Gene_search",
                     label = "Gene name, gene stable ID, or phrase",
                     placeholder = "mtnr1al, ENSAMXG00000010894, melatonin receptor, etc...",
                     btnSearch = icon("search"),
                     btnReset = icon("remove"),
                     width = "450px"
                   )
                   ),
                 mainPanel(id = "main",
                   textOutput("GeneCent_warnings"),
                   tableOutput("GeneCent_table")
                 )
               )
      ),
      tabPanel(h2("Transcription"), fluid = TRUE,
               sidebarLayout(
                 sidebarPanel(id = "sidebar",
                   selectInput(
                     inputId = "morph1",
                     label = "Compare...",
                     choices = c("Pachon","Molino","Tinaja","Rascon","Rio Choy")
                   ),
                   selectInput(
                     inputId = "morph2",
                     label = "to...",
                     choices = c("")
                   ),

                   # If morph2 IS set to control, populate a drop-down of available conditions
                   conditionalPanel(
                     condition = "input.morph2 == 'Control'",
                     selectInput(
                       inputId = "condition",
                       label = "Condition?",
                       choices = c("Sleep deprivation","Other")
                     )
                   ),
                   radioButtons("direction",
                                label = textOutput("dir_label"),
                                choices = c("Upregulated", "Downregulated")
                   ),
                   sliderInput(inputId = "percent_change",
                               label = textOutput("per_label"),
                               min = 0, max = 100, value = 50),
                   actionButton("Transc_enter","Find Genes")
                   ),

                 mainPanel(fluidRow(
                   conditionalPanel(condition = "Transc_enter",
                                    tableOutput("transc_table_out")
                   )
                 )
                 )
               )
      ),
        tabPanel(h2("Population Genetics"), fluid = TRUE,
                 sidebarLayout(
                   sidebarPanel(id = "sidebar2",
                     radioButtons("which_function",
                                  label = "Would you like to search for genes within a range of statistic values or search for genes using GO terms?",
                                  choices = c("Range of Statistic Values" = "distr_func", 
                                              "GO Terms" = "stat_by_chr_func")
                     ),
                   ),
                   mainPanel()
                 ),
                 conditionalPanel(
                   condition = "input.which_function == 'distr_func'",
                   sidebarLayout(
                     sidebarPanel(id = "sidebar",
                       radioButtons("type",
                                    label = "Search for Top/Bottom Number of Genes or Genes Above/Below a Statistical Value?",
                                    choices = c("Number of Genes" = "Gene Count", "Statistic Value")
                       ),
                       radioButtons("dist_statist",
                                    label = "Statistic of Interest",
                                    choices = c("Fst","Dxy","Tajima's D" = "TajimasD","Pi")),
                       checkboxGroupInput("dist_pops",
                                          label = "Population(s) of Interest",
                                          choices = c("Molino", "Pachon", "Rascon", "Rio Choy", "Tinaja")),
                       
                       # Only show this panel if the user wants to find the number of genes
                       # with the greatest or smallest values for the stat of interest
                       conditionalPanel(
                         condition = "input.type == 'Gene Count'",
                         sliderInput("gene_count", "How many genes would you like to see?",
                                     min = 1, max = 1000, value = 10),
                         radioButtons("TB",
                                      label = "Output genes with largest or smallest values of the desired statistic?",
                                      choices = c("Largest" = "top", "Smallest" ="bottom")),
                         actionButton("GCDistTable_enter","Find Genes"),
                       ),
                       # Only show this panel if the user wants to find genes with stat values
                       # above or below a specific threshhold
                       conditionalPanel(
                         condition = "input.type == 'Statistic Value'",
                         uiOutput("thrsh_slider"),
                         radioButtons("TB",
                                      label = "Output genes whose value is above or below the threshhold?",
                                      choices = c("Above" = "top", "Below" ="bottom")),
                         actionButton("SVDistTable_enter","Find Genes"),
                         actionButton("SVDistPlot_enter", "Visualize"),
                       )
                     ),
                     mainPanel(
                       conditionalPanel(
                         condition = "input.type == 'Gene Count'",
                         tableOutput("GCdist_tab"),
                         textOutput("GCdist_wrnings"),
                       ),
                       # Only show this panel if the user wants to find genes with stat values
                       # above or below a specific threshhold
                       conditionalPanel(
                         condition = "input.type == 'Statistic Value'",
                         plotOutput("SVdist_plot"),
                         textOutput("SVdist_plot_wrnings"),
                         tableOutput("SVdist_tab"),
                         textOutput("SVdist_wrnings"),
                       )
                     )
                   )
                 ),
                 conditionalPanel(
                   condition = "input.which_function == 'stat_by_chr_func'",
                   sidebarLayout(
                     sidebarPanel(id = "sidebar",
                       checkboxGroupInput("sbc_statist", 
                                          label = "Statistic of Interest",
                                          choices = c("Fst","Dxy","Tajima's D" = "TajimasD","Pi")),
                       checkboxGroupInput("sbc_pops", 
                                          label = "Population(s) of Interest",
                                          choices = c("Molino", "Pachon", "Rascon", "Rio Choy", 
                                                      "Tinaja")),
                       searchInput(
                         inputId = "GO_search", label = "GO ID or Term",
                         placeholder = "GO:0000001, mitochondrion, etc...",
                         btnSearch = icon("search"),
                         btnReset = icon("remove"),
                         width = "450px"
                       ),
                       selectInput(
                         inputId = "stat_PlotSelect",
                         label = "Visualize statistic...",
                         choices = "",
                         selected = NULL,
                         multiple = FALSE
                       ),
                       selectInput(
                         inputId = "scaff_PlotSelect",
                         label = "...plotted along scaffold:",
                         choices = "",
                         selected = NULL,
                         multiple = FALSE
                       ),
                       actionButton("SBCP_enter","Visualize")
                     ),
                     mainPanel(
                        tableOutput("SBC_table"),
                        textOutput("SBC_wrnings"),
                        plotOutput("SBC_plot")
                     )
                   )
                 )
        ),
      tabPanel(h2("GO Term Info"), fluid = TRUE,
               sidebarLayout(
                 sidebarPanel(id = "sidebar",
                   searchInput(
                     inputId = "GO_info_search",
                     label = "Phrase or comma-separated list of GO IDs",
                     placeholder = "GO:0000001, mitochondrion, etc...",
                     btnSearch = icon("search"),
                     btnReset = icon("remove"),
                     width = "450px"
                   )
                 ),
                 mainPanel(
                   tableOutput("GOinfo_table"),
                   textOutput("GOinfo_wrnings")
                 )
               )
      )
    )
  )

  server = function(input, output) {
    observe({
      transc_morph_choices <- c("Control", "Pachon","Molino","Tinaja","Rascon",
                                "Rio Choy")
      # Transcription Page: Update morph-selection widget to only enable 
      # comparisons between current morph and morph which is NOT morph1
      updateSelectInput(session = getDefaultReactiveDomain(),
                        "morph2",
                        choices = transc_morph_choices[transc_morph_choices != input$morph1],
                        selected = tail(transc_morph_choices, 1)
      )
      # Population Genetics (Stat-By-Chr Suppage): Update widget for selecting 
      # which statistic to plot
      updateSelectInput(session = getDefaultReactiveDomain(),
                        inputId = "stat_PlotSelect",
                        label = "Visualize statistic...",
                        choices = input$sbc_statist)
      # Population Genetics (Stat-By-Chr Suppage): If a table has been created, 
      # update the widget for selecting which scaffold to plot
      if(length(SBCT()) == 2){
        all_plots <- StatByChrGraph(SBCT()[[2]], stat_vec = input$sbc_statist)
        available_scaffs <- c()
        for(i in 1:length(all_plots)){
          available_scaffs <- append(available_scaffs,
                                     str_split(names(all_plots)[[i]], ":")[[1]][2])
        }
        updateSelectInput(session = getDefaultReactiveDomain(),
                          inputId = "scaff_PlotSelect",
                          label = "...plotted along scaffold:",
                          choices = available_scaffs)
      }else{
        updateSelectInput(session = getDefaultReactiveDomain(),
                          inputId = "scaff_PlotSelect",
                          label = "...plotted along scaffold:",
                          choices = "")
      }
      # Population Genetics (StatDist Subpage): Update min and max values for 
      # threshold slider 
      output$thrsh_slider <- renderUI({
        if(((length(input$dist_pops) >= 2) & 
         ((input$dist_statist == "Fst") | (input$dist_statist == "Dxy"))) | 
         ((length(input$dist_pops) < 2) & 
         ((input$dist_statist != "Fst") & (input$dist_statist != "Dxy")))){
        min_max_list <- MinMax(mm_pops = input$dist_pops,
                               mm_stat = input$dist_statist,
                               stat_table)
        sliderInput("thrsh", "Threshhold statistical value: ",
                    min = round(min_max_list[[1]],2), 
                    max = round(min_max_list[[2]],2), value = 0, step = 0.0005)
      }else{
        sliderInput("thrsh", 
                    "Threshold selector will appear here once you have inputted enough populations for the selected statistic.",
                    min = 0, max = 0, value = 0, 
                    step = 0.05)
      }
      })
    })
    
    
    # Home Page: Output text describing website, functions, data contribution,
    # and Astyanax mexicanus
    output$home_text2 <- renderText("CaveCrawler is a reactive web interface for bioinformatic analysis of data in the Mexican tetra (Astyanax mexicanus), an emerging evolutionary model organism.")
    output$home_text3 <- renderText("CaveCrawler consists of 4 subpages: a Gene Search page for querying data about specific genes, a Transcription page for finding genes whose transcriptional levels differ between samples, a Population Genetics page for investigating statistics on diversity and selection, and a GO Term Info page for identifying and obtaining information on GO terms-of-interest.")
    output$home_text4 <- renderText("To request that new data be integrated into CaveCrawler, please email Annabel Perry at annabelperry@tamu.edu")
    output$home_text5 <- renderText("To cite CaveCrawler, please cite our paper <link>")
    
    # Home Page: Plot map of all populations
    # Plot map of Astyanax populations
    pop_map <- ggplot(data = world_map) + 
      geom_polygon(aes(x = long,
                       y = lat, 
                       group = group),
                   fill = "antiquewhite",
                   color = "black") + 
      # get the portion of the map where your data is located at
      coord_fixed(xlim = c(min(Latit_Longit$Longitude) - 10, max(Latit_Longit$Longitude) + 10),
                  ylim = c(min(Latit_Longit$Latitude) - 5, max(Latit_Longit$Latitude) + 5),
                  ratio = 1.3)+
      # plot your points
      geom_point(data = Latit_Longit, 
                 aes(y = Latitude,
                     x = Longitude,
                     colour = factor(Population))) +
      # this will change colour to viridis colour palette
      scale_color_viridis_d("Population",alpha = .7,) +
      # change the axis labels
      xlab("Longitude") +
      ylab("Latitude") +
      # change the size of the points in the legend
      guides(colour = guide_legend(override.aes = list(size=5))) +
      # change the theme according to your taste
      annotate("text", y = 28, x = -98.5, label = "Texas") +
      annotate("text", y = 28, x = -102.5, label = "Mexico") +
      ggtitle("Geographic Locations of Astyanax mexicanus Populations") +
      theme_bw() + 
      theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
    
    output$home_plot <- renderPlot(pop_map)
    # Gene Search Page: Output a table of all statistics associated with the 
    # entered gene
    output$GeneCent_table <- renderTable({
      if(typeof(GeneCentOutput()) == "list"){
        GeneCentOutput()
      }else if(typeof(GeneCentOutput()) == "character"){
        data.frame()
      }
    })
    
    # Gene Search Page: Output any warnings associated with the searched gene
    output$GeneCent_warnings <- renderText({
      if(typeof(GeneCentOutput()) == "list"){
        "No warnings or errors"
      }else if(typeof(GeneCentOutput()) == "character"){
        GeneCentOutput()
      }
    })

    # Transcription Page: Create a label which changes based on the morphs to be 
    # searched for
    output$dir_label <- renderText({
      paste(c("Search for genes which are UP or DOWNregulated in ",
              input$morph1, " relative to ", input$morph2, "?"), collapse = "")
    })
    
    # Transcription Page: Create a label which changes based on whether up- or
    # downregulation was specified
    output$per_label <- renderText({
      if(input$direction == "Upregulated"){
        "Percent upregulation:"
      }else{
        "Percent downregulation:"
      }
    })
    
    # Transcription Page: Output a table with specified transcription data
    transc_table <- eventReactive(input$Transc_enter, valueExpr = {
      if(input$morph2 == "Control"){
        condit <- input$condition
      }else{
        condit <- "Between morph"
      }
      TranscTable(morph1 = input$morph1,
                  morph2 = input$morph2,
                  condition = condit,
                  direction = input$direction,
                  percent = input$percent_change,
                  GOTable = GeneToGO)
    })
    output$transc_table_out <- renderTable({
      transc_table()
    })
    
    # Population Genetics (Distribution Suppage): If Statistic Value was 
    # specified, output a table of all genes within the specified range of 
    # statistic values
    SVDT <- eventReactive(input$SVDistTable_enter, valueExpr = {
      StatDistTable(input$type,
                    input$TB,
                    input$dist_statist,
                    thresh = input$thrsh,
                    stat_table,
                    input$dist_pops)
    }
    )
    # Render table with 5 units per statistic value
    output$SVdist_tab <- renderTable(
      if(length(SVDT()) == 2){
        SVtemp_df <- data.frame(
          SVDT()[[2]][,1:2],
          format(SVDT()[[2]][,3], digits = 5),
          SVDT()[[2]][,4:7]
        )
        names(SVtemp_df) <- names(SVDT()[[2]])
        SVtemp_df
      }
    )
    output$SVdist_wrnings <- renderText(SVDT()[[1]])

    # Population Genetics (Distribution Suppage): If Visualize was pressed, 
    # output a plot of the number of genes with each value of the specified 
    # statistic 
    SVDP <- eventReactive(input$SVDistPlot_enter, valueExpr = {
      StatDistPlot(stat = input$dist_statist,
                   UL = input$TB,
                   thresh = input$thrsh,
                   stat_table,
                   pops = input$dist_pops)
    }
    )
    output$SVdist_plot_wrnings <- renderText(SVDP()[[1]])
    output$SVdist_plot <- renderPlot(SVDP()[[2]])

    # Population Genetics (Distribution Suppage): If Gene Count was specified,
    # output the statistics associated with the indicated number of genes
    GCDT <- eventReactive(input$GCDistTable_enter, valueExpr = {
      StatDistTable(input$type,
                    input$TB,
                    input$dist_statist,
                    thresh = input$gene_count,
                    stat_table,
                    input$dist_pops)
    }
    )
    # Output table with 5 decimals for the statistic column
    output$GCdist_tab <- renderTable(
      if(length(GCDT()) == 2){
        GCtemp_df <- data.frame(
          GCDT()[[2]][,1:2],
          format(GCDT()[[2]][,3], digits = 5),
          GCDT()[[2]][,4:7]
        )
        names(GCtemp_df) <- names(GCDT()[[2]])
        GCtemp_df
      }
    )
    output$GCdist_wrnings <- renderText(GCDT()[[1]])
    
    # Population Genetics (Stat-By-Chr Subpage): If population, statistic, and 
    # GO term have been entered, output a table of associated genes
    SBCT <- eventReactive(input$GO_search, valueExpr = {
      if((length(input$sbc_statist) == 0) & (length(input$sbc_pops) != 0)){
        "Please choose at least one statistic of interest"
      }else if((length(input$sbc_statist) != 0) & (length(input$sbc_pops) == 0)){
        "Please choose at least one population of interest"
      }else if((length(input$sbc_statist) == 0) & (length(input$sbc_pops) == 0)){
        "Please choose at least one statistic and population of interest"
      }else if((length(input$sbc_statist) != 0) & (length(input$sbc_pops) != 0)){
        StatByChrTable(GOTerm = input$GO_search,
                       GeneToGO,
                       GoIDToNames, 
                       UpperLower, 
                       stat_vec = input$sbc_statist, 
                       position_table, 
                       stat_table, 
                       pops = input$sbc_pops
        )
      }
    }
    )
    
    # Population Genetics (Stat-By-Chr Subpage): If visualize was pressed and 
    # input table is NOT full of NAs, output a plot of the appropriate statistic x 
    # scaffold pair
    SBCP <- eventReactive(input$SBCP_enter, valueExpr = {
      if((!is.na(SBCT()[[2]][1,1])) & (nrow(SBCT()[[2]]) != 1)){
        SBC_plots <- StatByChrGraph(SBCT()[[2]], stat_vec = input$sbc_statist)
        for(i in 1:length(SBC_plots)){
          if((str_split(names(SBC_plots)[[i]], ":")[[1]][1] == input$stat_PlotSelect)
             & (str_split(names(SBC_plots)[[i]], ":")[[1]][2] == input$scaff_PlotSelect)){
            plot_out <- SBC_plots[[i]]
          }
        }
        plot_out
      }
    })
    output$SBC_table <- renderTable(
      if(length(SBCT()) == 2){
        temp_df <- data.frame(
          SBCT()[[2]][,1:7],
          format(SBCT()[[2]][,8], digits = 5),
          SBCT()[[2]][,9]
        )
        names(temp_df) <- names(SBCT()[[2]])
        temp_df
      }else{
        data.frame(Data = c("No data to display; see explanation below"))
      })
    output$SBC_wrnings <- renderText(SBCT()[[1]])
    output$SBC_plot <- renderPlot(SBCP())
    
    # GO Term Info: If GO ID or phrase was inputted, output class, lower-level
    # GO IDs, and GO terms associated with all relevant GO IDs
    GOInfoOutWarnings <- eventReactive(input$GO_info_search, valueExpr = {
      if(input$GO_info_search == ""){
        "Warnings will populate here once GO ID or phrase is inputted."
      }else{
        GOInfo(
          GO_input = input$GO_info_search, 
           GO_classes, 
           GOIDToNames, 
           UpperLower, 
           all.GO_IDs
          )[[1]]
      }
    })
    GOInfoOutTable <- eventReactive(input$GO_info_search, valueExpr = {
      if(input$GO_info_search == ""){
        data.frame(`Column Name` = "Data will populate here once GO ID or phrase is inputted.")
      }else{
        GOInfo(
          GO_input = input$GO_info_search, 
          GO_classes, 
          GOIDToNames, 
          UpperLower, 
          all.GO_IDs
        )[[2]]
      }
    })
    output$GOinfo_table <- renderTable(
        GOInfoOutTable()
    )
    output$GOinfo_wrnings <- renderText(
        GOInfoOutWarnings()
    )
    
    
  }


# Run the application
shinyApp(ui = ui, server = server)
