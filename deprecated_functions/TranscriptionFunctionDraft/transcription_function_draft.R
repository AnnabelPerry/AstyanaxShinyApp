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

ui <- fluidPage(
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
    actionButton("Transc_enter","Find Genes"),
    conditionalPanel(condition = "Transc_enter",
                     tableOutput("test")
                     )
)

server <- function(input, output) {
    observe({
        all_choices <- c("Control", "Pachon","Molino","Tinaja","Rascon",
                         "Rio Choy")
        # Update widget to only enable comparisons between current morph and
        # morph which is NOT morph1
        updateSelectInput(session = getDefaultReactiveDomain(), 
                          "morph2",
                          choices = all_choices[all_choices != input$morph1],
                          selected = tail(all_choices, 1)
        )
    })
    
    condition_control <- read.csv("data/Morph_Control_TranscData.csv")
    colnames(condition_control)[1] <- gsub('^...','',colnames(condition_control)[1])
    morph1.morph2 <- read.csv("data/Toy_RioChoyPachon.csv")
    
    GeneToGO <- read.csv("data/AMexGOTerms.csv", fill = T)
    GeneToGO <- GeneToGO[GeneToGO$Gene.names != "",]
    GeneToGO$Gene.names <- tolower(GeneToGO$Gene.names)
    
    
    TranscTable <- function(morph1, morph2, condition, direction, percent, 
                            GOTable){
        # If condition is NOT "Between morph"...
        if(condition != "Between morph"){
            # Use transcription data of morph-control comparisons
            in_table <- condition_control
            # Store comparison
            comp <- paste(c(morph1, "Control"), collapse = "-")
            # If upregulated genes were requested...
            if(direction == "Upregulated"){
                # Find [percent]% of genes falling on morph(s)-of-interest with HIGHEST 
                # logFC scores AND p-value < 0.05
                
                # Find all rows-of-interest (ROIs) for morph-of-interest where genes are 
                # upregulated, condition matches the input specification, and p-value is 
                # less than 0.05
                ROIs <- in_table[(grepl(morph1, in_table$Class) & (in_table$logFC > 0) & 
                                      (in_table$Condition == condition) & (in_table$PValue < 0.05)), ]
                # Sort candidate rows with highest logFC values on top
                ROIs <- ROIs[order(ROIs[,3], decreasing = T),]
                # Find the number of rows corresponding to the specified percent
                n.rows <- as.integer((percent/100)*nrow(ROIs))
                # If downregulated genes were requested...
            }else if(direction == "Downregulated"){
                # Find [percent]% of genes falling on morph(s)-of-interest with LOWEST 
                # logFC scores AND p-value < 0.05
                
                # Find all rows-of-interest (ROIs) for morph-of-interest where genes are 
                # downregulated, condition matches the input specification, and p-value is 
                # less than 0.05
                ROIs <- in_table[(grepl(morph1, in_table$Class) & (in_table$logFC < 0) & 
                                      (in_table$Condition == condition) & (in_table$PValue < 0.05)), ]
                # Sort candidate rows with LOWEST logFC values on top
                ROIs <- ROIs[order(ROIs[,3], decreasing = F),]
                # Find the number of rows corresponding to the specified percent
                n.rows <- as.integer((percent/100)*nrow(ROIs))
            }
            # Obtain GO terms for ROIs
            GOTerms <- character(length = n.rows)
            for(i in 1:n.rows){
                if(tolower(ROIs$Gene_name)[i] %in% GeneToGO$Gene.names){
                    GOTerms[i] = GeneToGO$Gene.ontology.IDs[GeneToGO$Gene.names == tolower(ROIs$Gene_name)[i]]
                }else{
                    GOTerms[i] = NA
                }
            }
            # Output gene names, gene stable IDs, GO terms, morph-of-comparison, logFC, 
            # p-value, and Ensembl family information to a dataframe
            output.df <- data.frame(
                tolower(ROIs$Gene_name)[1:n.rows],
                ROIs[1:n.rows,1],
                GOTerms,
                rep(comp, n.rows)[1:n.rows],
                ROIs$logFC[1:n.rows],
                ROIs$PValue[1:n.rows],
                ROIs$Ensembl_Family_Description[1:n.rows]
            )
            names(output.df) <- c(
                "Gene Name",
                "Gene Stable ID",
                "GO Term(s)",
                "Comparison",
                "logFC",
                "p-value",
                "Ensembl Family Description"
            )
            
            # If condition is "Between morph"...
        }else if(condition == "Between morph"){
            # Use transcription data for between-morph comparisons
            in_table <- morph1.morph2
            
            # Find rows corresponding to morphs of interest
            morph1.rows <- in_table[in_table$Class == morph1,]
            morph2.rows <- in_table[in_table$Class == morph2,]
            # Find names of all genes present for morph1
            m1.genes <- morph1.rows$Gene_name[morph1.rows$Gene_name != ""]
            # Subtract all morph2 logFC values FROM all same-gene morph1 logFC value
            delta_FC <- c()
            genes <- c()
            G_IDs <- c()
            EF_IDs <- c()
            for(g in 1:length(m1.genes)){
                if(m1.genes[g] %in% morph2.rows$Gene_name){
                    diff <- morph1.rows$logFC[morph1.rows$Gene_name == m1.genes[g]] - 
                        morph2.rows$logFC[morph2.rows$Gene_name == m1.genes[g]]
                    delta_FC <- append(delta_FC, diff)
                    genes <- append(genes, m1.genes[g])
                    G_IDs <- append(G_IDs, morph1.rows$ï..Gene_stable_ID[
                        morph1.rows$Gene_name == m1.genes[g]])
                    EF_IDs <- append(EF_IDs, morph1.rows$Ensembl_Family_Description[
                        morph1.rows$Gene_name == m1.genes[g]])
                    # If morph1 and morph2 do NOT have matching genes, skip this gene
                }else{
                    next
                }
            }
            # Create dataframe to easily access the change in log FC values
            # corresponding to each gene pair
            ROIs <- data.frame(genes,G_IDs,delta_FC,EF_IDs)
            # If genes upregulated in morph1 were requested...
            if(direction == "Upregulated"){
                # Find genes with greatest POSITIVE change in expression between 
                # morphs
                ROIs <- ROIs[order(delta_FC, decreasing = T),]
                n.rows <- as.integer((percent/100)*nrow(ROIs))
                # If genes downregulated in morph1 were requested...
            }else if(direction == "Downregulated"){
                # Find indices of rows with greatest NEGATIVE change in expression between
                # morphs
                ROIs <- ROIs[order(delta_FC, decreasing = F),]
                n.rows <- as.integer((percent/100)*nrow(ROIs))
            }
            # Find GO terms for genes-of-interest
            GOTerms <- character(length = n.rows)
            for(i in 1:n.rows){
                if(tolower(ROIs$genes)[i] %in% GeneToGO$Gene.names){
                    GOTerms[i] = GeneToGO$Gene.ontology.IDs[GeneToGO$Gene.names == tolower(ROIs$genes)[i]]
                }else{
                    GOTerms[i] = NA
                }
            }
            # Output gene names, gene stable IDs, GO terms, direction, delta(logFC),
            # and Ensembl family description
            output.df <- data.frame(
                tolower(ROIs$genes)[1:n.rows],
                ROIs$G_IDs[1:n.rows],
                GOTerms,
                rep(paste(c(direction, " in ", morph1, " vs. ", morph2), collapse = ""),
                    n.rows),
                ROIs$delta_FC[1:n.rows],
                ROIs$EF_IDs[1:n.rows]
            )
            names(output.df) <- c(
                "Gene Name",
                "Gene Stable ID",
                "GO Term(s)",
                "Comparison",
                "delta(logFC)",
                "Ensembl Family Description"
            )
        }
        return(output.df)
    }
    
    output$dir_label <- renderText({
        paste(c("Search for genes which are UP or DOWNregulated in ", 
                input$morph1, " relative to ", input$morph2, "?"), collapse = "")
    })
    output$per_label <- renderText({
        if(input$direction == "Upregulated"){
            "Percent upregulation:"
        }else{
            "Percent downregulation:"
        }
    })
    testing <- eventReactive(input$Transc_enter, valueExpr = {
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
    output$test <- renderTable({
        testing()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
