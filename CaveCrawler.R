#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

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
                   downloadButton("GeneSearchDL", "Download", class = "download"),
                   textOutput("GeneCent_warnings"),
                   tags$head(tags$style(".download{background-color:#c8feca;} .download{color: #71c596 !important;} .download{border-color: #71c596 !important;}")),
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
                                    downloadButton("TranscDL", "Download", class = "download"),
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
                                  label = "Would you like to find genes which are outliers with respect to a statistic-of-interest or find the statistic values for all genes related to a GO-term-of-interest?",
                                  choices = c("Outliers" = "distr_func", 
                                              "GO-term-of-interest" = "stat_by_chr_func")
                     ),
                   ),
                   mainPanel(
                     br(),
                     br(),
                     br(),
                     br(),
                     br(),
                     br(),
                     br(),
                     br(),
                     # Display download button appropriate to specific conditions
                     conditionalPanel(
                       condition = "input.which_function == 'stat_by_chr_func'",
                       downloadButton("SBCDL", "Download", class = "download"),
                     ),
                     conditionalPanel(
                       condition = "input.which_function == 'distr_func'",
                       conditionalPanel(
                         condition = "input.type == 'Gene Count'",
                         downloadButton("GCDistDL", "Download", class = "download"),
                       ),
                       conditionalPanel(
                         condition = "input.type == 'Statistic Value'",
                         downloadButton("SVDistDL", "Download", class = "download"),
                       )
                     )
                   )
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
                                          choices = c("Molino", "Pachon", 
                                                      "Rascon", "Rio Choy", 
                                                      "Tinaja", 
                                                      "Chica 1" = "Chica1",
                                                      "Chica 2" = "Chica2")),
                       
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
                                          choices = c("Molino", "Pachon", 
                                                      "Rascon", "Rio Choy", 
                                                      "Tinaja", 
                                                      "Chica 1" = "Chica1",
                                                      "Chica 2" = "Chica2")),
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
                       plotOutput("SBC_plot"),
                        tableOutput("SBC_table"),
                        textOutput("SBC_wrnings")
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
                   downloadButton("GOInfoDL", "Download", class = "download"),
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
      # Population Genetics (Stat-By-Chr Suppage): If a valid table has been 
      # created, update widget for selecting which statistic to plot
      if(length(SBCT()) == 2){
        if(sum(is.na(SBCT()[[2]])) != 9){
          updateSelectInput(session = getDefaultReactiveDomain(),
                            inputId = "stat_PlotSelect",
                            label = "Visualize statistic...",
                            choices = input$sbc_statist)
        }else{
          updateSelectInput(session = getDefaultReactiveDomain(),
                            inputId = "stat_PlotSelect",
                            label = "Visualize statistic...",
                            choices = "")
      }
      }else{
        updateSelectInput(session = getDefaultReactiveDomain(),
                          inputId = "stat_PlotSelect",
                          label = "Visualize statistic...",
                          choices = "")
      }
      # Population Genetics (Stat-By-Chr Suppage): If a valid table has been 
      # created, update the widget for selecting which scaffold to plot
      if(length(SBCT()) == 2){
          if(sum(is.na(SBCT()[[2]])) != 9){
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
    GeneCentOutput <- eventReactive(input$Gene_search, valueExpr = {
      if(input$Gene_search == ""){
        data.frame()
      }else{
        GeneCentered(input$Gene_search,
                     stat_table,
                     GeneToGO,
                     condition_control,
                     position_table)
      }
    })
    
    output$GeneCent_table <- renderTable({
      if(typeof(GeneCentOutput()) == "list"){
        # Format each statistic section to have appropriate number of decimals
        reformattedGeneCent <- data.frame(
          GeneCentOutput()[,1:7],
          format(GeneCentOutput()[,8:37], digits = 5),
          GeneCentOutput()[,38]
        )
        names(reformattedGeneCent) <- names(GeneCentOutput())
        reformattedGeneCent
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

    # Gene Search Page: Enable downloading of Gene Search table
    output$GeneSearchDL <- downloadHandler(
      filename = function() {
        paste("CaveCrawler-GeneSearch-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        # If a valid table was outputted. enable downloading of that table
        if(typeof(GeneCentOutput()) == "list"){
          reformattedGeneCent <- data.frame(
            GeneCentOutput()[,1:7],
            format(GeneCentOutput()[,8:37], digits = 5),
            GeneCentOutput()[,38]
          )
          names(reformattedGeneCent) <- names(GeneCentOutput())
        # If a valid table was not outputted, enable downloading of an empty df
        }else{
          reformattedGeneCent <- data.frame()
        }
        
        write.csv(reformattedGeneCent, file, row.names = F)
      }
    )
    
    
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
      # Change log fold changes to have 5 decimal points
      reformattedTranscT <- data.frame(
        transc_table()[,1:4],
        format(transc_table()[,5], digits = 5),
        transc_table()[,6]
      )
      names(reformattedTranscT) <- names(transc_table())
      reformattedTranscT
    })
    
    # Transcription Page: Enable downloading of Transcription table
    output$TranscDL <- downloadHandler(
      filename = function() {
        paste("CaveCrawler-Transcription-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(transc_table(), file, row.names = F)
      }
    )
    
    
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

    # Statistic Value SubPage: Enable downloading of Statistic Value table
    output$SVDistDL <- downloadHandler(
      filename = function() {
        paste("CaveCrawler-Outliers-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(SVDT()[[2]], file, row.names = F)
      }
    )
    
    # Population Genetics (Distribution Subpage): If Gene Count was specified,
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
    
    # Gene Count SubPage: Enable downloading of Gene Count table
    output$GCDistDL <- downloadHandler(
      filename = function() {
        paste("CaveCrawler-Outliers-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(GCDT()[[2]], file, row.names = F)
      }
    )
    
    
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
      if(sum(is.na(SBCT()[[2]])) != 9){
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
      # First, check if anything has been inputted into function
      if(length(SBCT()) == 2){
        # If function has inputs, check if inputs yielded a valid table
        if(sum(is.na(SBCT()[[2]])) != 9){
          # If inputs did yield a valid table, output the table with adjusted 
          # decimal places
          temp_df <- data.frame(
            SBCT()[[2]][,1:7],
            format(SBCT()[[2]][,8], digits = 5),
            SBCT()[[2]][,9]
          )
          names(temp_df) <- names(SBCT()[[2]])
          temp_df
        }else{
          # If inputs yielded an erroneous table, simply output an empty table
          data.frame(Data = "NA")
        }
      }else{
        data.frame(Data = c("No data to display; see explanation below"))
      })
    output$SBC_wrnings <- renderText(SBCT()[[1]])
    output$SBC_plot <- renderPlot(SBCP())
    
    # Statistic-by-Chromosome SubPage: Enable downloading of SBC table
    output$SBCDL <- downloadHandler(
      filename = function() {
        paste("CaveCrawler-Outliers-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(SBCT()[[2]], file, row.names = F)
      }
    )
    
    
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
    # GO Info Page: Enable downloading of GO Info table
    output$GOInfoDL <- downloadHandler(
      filename = function() {
        paste("CaveCrawler-GOInfo-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        # Do not write row numbers to CSV
        write.csv(GOInfoOutTable(), file, row.names = F)
      }
    )
    
  }


# Run the application
shinyApp(ui = ui, server = server)
