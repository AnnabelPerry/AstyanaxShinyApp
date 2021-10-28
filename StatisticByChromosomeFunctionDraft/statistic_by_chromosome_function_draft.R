#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyWidgets)
library(ggplot2)

ui <- fluidPage(
  checkboxGroupInput("statist", 
               label = "Statistic of Interest",
               choices = c("Fst","Dxy","Tajima's D" = "TajimasD","Pi")),
  checkboxGroupInput("pops", 
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
    inputId = "SBCP_select",
    label = "Visualize Scaffold...",
    choices = "",
    selected = NULL,
    multiple = FALSE
  ),
  tableOutput("test2"),
  textOutput("test1")
)


server <- function(input, output) {
  pos_table <- read.csv("../data/AmexPositionTable.csv", fill = TRUE)
  
  GeneToGO <- read.csv("../data/AMexGOTerms.csv", fill = T)
  GeneToGO <- GeneToGO[GeneToGO$Gene.names != "",]
  # Convert gene names in GeneToGO to lowercase to make compatible with other tables
  GeneToGO$Gene.names <- tolower(GeneToGO$Gene.names)
  
  GoID_Names <- read.table("../data/GOIDs_and_Names.txt", fill = T, sep = "\t", header = T)
  Up_Low <- read.table("../data/GOTermAssociations.txt", fill = T, sep = "\t", header = T)
  
  s_table <- read.csv("../data/AMexicanus_Genes_and_Stats.csv")
  s_table <- s_table[,(names(s_table) != "X")]
  
  StatByChrTable <- function(GOTerm, GeneToGo, GoIDToNames, UpperLower, stat_vec, position_table, stat_table, all_pops){
    # Initialize vector for all GO terms of interest
    GOs <- c()
    # Initialize vector in which to store warnings
    wrnings <- c("Notes: ")
    
    # Check if user inputted a word/phrase or a GO ID. If the user inputted a 
    # word/phrase, find the name(s) which contains that word/phrase and set the GO
    # vector equal to the corresponding GO IDs.
    if(!(GOTerm %in% GoIDToNames$GO.ID)){
      GOs <- GoIDToNames$GO.ID[which(grepl(GOTerm, GoIDToNames$GO.Term, 
                                           ignore.case = T))]
      # If user inputted a GO ID, add that ID 
    }else if(GOTerm %in% GoIDToNames$GO.ID){
      GOs <- GOTerm
    }
    
    # Find all GO terms which are LOWER than the GO terms in the GO vector.
    # First, loop through each GO ID in the vector-of-GOs. As more GO IDs, are 
    # added to the vector, the number of remaining iterations will increase.
    for(g in 1:length(GOs)){
      # Add all "Lower" GO IDs which occur on a row where the current vector entry
      # is an "Upper" to the vector of GO IDs, then move to the next GO ID
      if(GOs[g] %in% UpperLower$Upper){
        GOs <- append(GOs, UpperLower$Lower[UpperLower$Upper == GOs[g]])
        # If the current GO ID does NOT occur anywhere in the "Upper" column, skip it
      }else{
        next
      }
    }
    # Find all genes associated with the current GO IDs and add to vector of genes
    gene_vec <- c()
    found_GOs <- c()
    for(g in 1:length(GOs)){
      # If the GO term appears in the data frame of names AND the corresponding gene
      # occurs in the statistics vector, add the GO term and gene name
      if(GOs[g] %in% GeneToGO$Gene.ontology.IDs){
        gene_vec <- append(gene_vec, GeneToGO$Gene.names[grepl(GOs[g], GeneToGO$Gene.ontology.IDs)])
        found_GOs <- append(found_GOs, rep(GOs[g], length(gene_vec)))
        # If the GO term does NOT appear in the dataframe of names, skip it
      }else{
        next
      }
    }
    geneGOs <- data.frame(
      Gene = gene_vec,
      GO_ID = found_GOs
    )
    geneGOs <- geneGOs[!duplicated(geneGOs), ]
    
    # Extract all possible combinations of populations from populations of interest
    if((is.vector(all_pops)) & (length(all_pops) > 1)){
      all_pops <- combn(all_pops,2) 
    }
    
    # Create vectors in which to store values for later dataframe
    Statistic_Type_prelim <- c()
    Population_prelim <- c()
    Statistic_Value_prelim <- numeric()
    
    for(s in 1:length(stat_vec)){
      # Check if statistic of interest makes comparisons between TWO populations
      if((stat_vec[s] == "Fst") | (stat_vec[s] == "Dxy")){
        # Check if pops is a matrix
        if(is.matrix(all_pops)){
          # If so, iterate through each combination and output stat value for that
          # combination
          for(pair in 1:ncol(all_pops)){
            
            # If pops is a matrix, read the strings, find the column corresponding 
            # to the stat of interest for the populations in the vector, and set 
            # "index" equal to the column housing this statistic
            val <- which(grepl(all_pops[1, pair], names(stat_table))
                         & grepl(all_pops[2, pair], names(stat_table))
                         & grepl(stat_vec[s], names(stat_table)))
            # If statistic for populations-of-interest is not present, return error
            if(length(val) == 0){
              wrnings <- append(wrnings, paste(c("Statistic ",stat_vec[s],
                                                 " is not present for the populations ",
                                                 all_pops[1, pair]," and ",
                                                 all_pops[2, pair], " | "),collapse = ""))
            }else{
              # Create a row to add to the indices dataframe
              temp_str <- paste(c(all_pops[1, pair],"-",all_pops[2, pair]), 
                                collapse = "")
              Statistic_Type_prelim <- append(Statistic_Type_prelim,stat_vec[s])
              Population_prelim <- append(Population_prelim,temp_str)
              Statistic_Value_prelim <- append(Statistic_Value_prelim,val)
            }
            # If pops is NOT a matrix, return an error
          }
        }else if(!is.matrix(all_pops)){
          null.df <- data.frame(matrix(nrow = 1, ncol = 8))
          names(null.df) <- c("Gene",
                              "Scaffold",
                              "Start_Position",
                              "End_Position",
                              "GO_IDs",
                              "Statistic_Type",
                              "Population",
                              "Statistic_Value")
          return(list(paste(c("Warning: Only one population, ", all_pops, 
                              ", supplied for the two-population statistic ", 
                              stat_vec[s], "."),
                            collapse = ""), null.df))
          
          # Check if statistic of interest makes comparisons between ONE population
        }
      }else if((stat_vec[s] == "TajimasD") | (stat_vec[s] == "Pi")){
        # Iterate through each individual population
        for(p in 1:length(all_pops)){
          # If pops is a string, set indx equal to the column housing the stat of
          # interest for this population
          val <- which(grepl(all_pops[p], names(stat_table))
                       & grepl(stat_vec[s], names(stat_table)))
          # If statistic for populations-of-interest is not present, return a 
          # warning 
          if(length(val) == 0){
            wrnings <- append(wrnings, paste(c("Statistic ",stat_vec[s],
                                               " is not present for the population ",
                                               all_pops[p], " | "), collapse = ""))
          }else{
            # Create a row to add to the indices dataframe
            Statistic_Type_prelim <- append(Statistic_Type_prelim,stat_vec[s])
            Population_prelim <- append(Population_prelim,all_pops[p])
            Statistic_Value_prelim <- append(Statistic_Value_prelim,val)
          }
        }
      }
    }
    
    # For each GO Term in the vector of GO Terms-of-interest, find...
    # 1. All gene names which occur in the positions table AND in the statistics
    #    table
    # 2. The scaffolds of those genes
    # 3. The starting positions of those genes
    # 4. The ending positions of those genes
    # 5. The GO terms associated with those genes
    # 6. Each of the 4 statistic types
    # 7. Each of the populations/population combinations
    # 8. Each of the statistical values
    
    # For each statistic-population pair, iterate through each gene and find the 
    # statistic value, scaffold, starting position, and ending position and
    # output to a dataframe
    Gene <- c()
    Scaffold <- c()
    Start_Position <- c()
    End_Position <- c()
    GO_IDs <- c()
    Statistic_Type <- c()
    Population <- c()
    Statistic_Value <- c()
    
    for(s in 1:length(Statistic_Type_prelim)){
      for(g in 1:length(geneGOs$Gene)){
        # Check if current gene is in stat AND position table
        # If so...
        if((geneGOs$Gene[g] %in% stat_table$Gene_Name) & 
           (geneGOs$Gene[g] %in% position_table$Gene_Name)){
          # Output the current gene
          Gene <- append(Gene, geneGOs$Gene[g])
          # Output scaffold of current gene
          Scaffold <- append(Scaffold, 
                             position_table$Scaffold[position_table$Gene_Name == geneGOs$Gene[g]])
          # Output starting position of the current gene
          Start_Position <- append(Start_Position, 
                                   position_table$Start_Locus[position_table$Gene_Name == geneGOs$Gene[g]])
          # Output the ending position of the current gene
          End_Position <- append(End_Position, 
                                 position_table$End_Locus[position_table$Gene_Name == geneGOs$Gene[g]])
          # Output ALL GO terms associated with the current gene
          GO_IDs <- append(GO_IDs,
                           paste(geneGOs$GO_ID[geneGOs$Gene == geneGOs$Gene[g]], collapse = "; "))
          # Output the current statistic type
          Statistic_Type <- append(Statistic_Type, Statistic_Type_prelim[s])
          # Output the population(s)
          Population <- append(Population, Population_prelim[s])
          # Output the statistic value
          Statistic_Value <- append(Statistic_Value, 
                                    stat_table[stat_table$Gene_Name == geneGOs$Gene[g], 
                                               Statistic_Value_prelim[s]])
          # If not, skip the gene
        }else{
          next
        }
      }
    }
    output_df <- data.frame(Gene,
                            Scaffold,
                            Start_Position,
                            End_Position,
                            GO_IDs,
                            Statistic_Type,
                            Population,
                            Statistic_Value
    )
    return(list(wrnings, output_df))
  }
  StatByChrGraph <- function(Full_Table, stat_vec){
    
    # Count the number of unique scaffolds in the table
    unique_scaffs <- levels(as.factor(Full_Table$Scaffold))
    num_scaff <- length(unique_scaffs)
    
    # Find number of unique statistics so you can create a distinct plot for each
    # statistic-scaffold pair
    unique_stats <- levels(as.factor(Full_Table$Statistic_Type))
    num_stats <- length(unique_stats)
    
    # Initialize a list in which to store all plots
    plist <- vector(mode = "list", length = num_scaff*num_stats)
    
    # Create vector of names
    names_vec <- c()
    # Start at the 0th entry in the list
    p = 0
    # For each scaffold-statistic combo, create a new plot
    for(s in 1:num_stats){
      # For each scaffold represented by the selected genes, create a new plot
      for(c in 1:num_scaff){
        # Each time you reach a new scaffold-statistic combination, move to a new 
        # position in the plot list 
        p = p + 1
        # Create a new name from this combination
        names_vec <- append(names_vec, paste(c(stat_vec[s],unique_scaffs[c]), 
                                             collapse = ":"))
        # Find all rows in the table whose basic statistic and scaffold match the
        # current statistic and scaffold
        plot_dat <- data.frame(
          Locus = Full_Table$Start_Position[(Full_Table$Scaffold == unique_scaffs[c])
                                            & (Full_Table$Statistic_Type == unique_stats[s])],
          Value = Full_Table$Statistic_Value[(Full_Table$Scaffold == unique_scaffs[c])
                                             & (Full_Table$Statistic_Type == unique_stats[s])],
          Populations = Full_Table$Population[(Full_Table$Scaffold == unique_scaffs[c])
                                              & (Full_Table$Statistic_Type == unique_stats[s])]
        )
        
        # Plot the statistic values for all genes which occur on the scaffold
        # of interest, colored by population or population pair
        plist[[p]] <- ggplot(plot_dat, aes(x=Locus, y=Value, color=Populations)) + 
          ylab(stat_vec[s]) +
          xlab("Locus") +
          ggtitle(paste(c("Scaffold:", unique_scaffs[c]), collapse = " ")) +
          geom_point() +
          theme_bw()
        #Store the plot in a vector
      }
    }
    # Add names to plot list
    names(plist) <- names_vec
    return(plist)
  }
  
  
  SBCT <- eventReactive(input$GO_search, valueExpr = {
    if((length(input$statist) == 0) & (length(input$pops) != 0)){
      "Please choose at least one statistic of interest"
    }else if((length(input$statist) != 0) & (length(input$pops) == 0)){
      "Please choose at least one population of interest"
    }else if((length(input$statist) == 0) & (length(input$pops) == 0)){
      "Please choose at least one statistic and population of interest"
    }else if((length(input$statist) != 0) & (length(input$pops) != 0)){
      StatByChrTable(GOTerm = input$GO_search,
                     GeneToGO,
                     GoIDToNames = GoID_Names, 
                     UpperLower = Up_Low, 
                     stat_vec = input$statist, 
                     position_table = pos_table, 
                     stat_table = s_table, 
                     all_pops = input$pops
                     )
    }
  }
  )
  SBCP <- eventReactive(input$GO_search, valueExpr = {
    StatByChrTable_input <- SBCT()[[2]]
    updateSelectInput(session = getDefaultReactiveDomain(),
                      inputId = "SBCP_select",
                      label = "Visualize Scaffold...",
                      choices = levels(as.factor(StatByChrTable_input$Scaffold)))
  }
  )
  output$test2 <- renderTable(SBCT()[[2]])
  output$test1 <- renderText(SBCT()[[1]])
  
}

# Run the application 
shinyApp(ui = ui, server = server)
