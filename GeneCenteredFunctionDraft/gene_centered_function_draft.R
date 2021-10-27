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

setwd("~/Fall 2021/Capstone")
dat <- read.table("Astyanax_mexicanus.Astyanax_mexicanus-2.0.104.gtf", fill = TRUE, skip = 5)
position_table <- dat[dat$V3 == "gene",c(1,4,5,10,16)]

condition_control <- read.csv("ShinyInputData/Morph_Control_TranscData.csv")

GeneToGO <- read.csv("ShinyInputData/AMexGOTerms.csv", fill = T)
GeneToGO <- GeneToGO[GeneToGO$Gene.names != "",]
GeneToGO$Gene.names <- tolower(GeneToGO$Gene.names)

stat_table <- read.csv("ShinyInputData/AMexicanus_Genes_and_Stats.csv")
stat_table <- stat_table[,(names(stat_table) != "X")]


ui <- fluidPage(
    searchInput(
        inputId = "Gene_search", label = "Gene name, gene stable ID, or phrase",
        placeholder = "mtnr1al, ENSAMXG00000010894, melatonin receptor, etc...",
        btnSearch = icon("search"),
        btnReset = icon("remove"),
        width = "450px"
    ),
    textOutput("GeneCent_warnings"),
    tableOutput("GeneCent_table"),
)

server <- function(input, output) {
    GeneCentered <- function(input, stat_table, GeneToGO, condition_control,
                             position_table){
        library("stringr")
        comma <- ", "
        # Obtain complete dataframe of all possible genes and corresponding IDs across 
        # the statistic and transcription data
        all.genes_IDs <- data.frame(
            all_genes = c(position_table$V16, stat_table$Gene_Name, 
                          condition_control$Gene_name
            ),
            all_IDs = c(position_table$V10, stat_table$Stable_Gene_ID, 
                        condition_control$ï..Gene_stable_ID
            )
        )
        # Since some genes are found in one or more of the input dfs, remove all 
        # duplicate rows
        all.genes_IDs <- all.genes_IDs[!duplicated(all.genes_IDs[,1]),]
        all.genes_IDs <- all.genes_IDs[!duplicated(all.genes_IDs[,2]),]
        
        # If input is a comma-separated string, parse string and separate elements, 
        # then determine whether vector contains gene names or gene IDs
        if(grepl(comma, input)){
            input_vec <- str_split(input, pattern = comma)[[1]]
            # Dataframe in which to store inputs and associated values
            output.df <- data.frame(matrix(nrow = length(input_vec), ncol = 35))
            
            # Next, iterate through each element of vector
            for(i in 1:length(input_vec)){
                # If element is present in all_genes, consider vector a vector of gene
                # names and output all associated elements to output dataframe
                if(input_vec[i] %in% all.genes_IDs$all_genes){
                    output.df[i,1] = input_vec[i]
                    output.df[i,2] = all.genes_IDs$all_IDs[all.genes_IDs$all_genes == input_vec[i]]
                    
                    # If the current gene is present in the position table, output position 
                    # table info. If not, output all NA
                    if(input_vec[i] %in% position_table$V16){
                        output.df[i,3] = position_table$V1[position_table$V16 == input_vec[i]]
                        output.df[i,4] = position_table$V4[position_table$V16 == input_vec[i]]
                        output.df[i,5] = position_table$V5[position_table$V16 == input_vec[i]]
                    }else{
                        output.df[i,3] = NA
                        output.df[i,4] = NA
                        output.df[i,5] = NA
                    }
                    
                    # Check if gene is present in GO term table. If so, output GO terms. If 
                    # not, output NA
                    if(input_vec[i] %in% GeneToGO$Gene.names){
                        output.df[i,7] = GeneToGO$Gene.ontology.IDs[GeneToGO$Gene.names == input_vec[i]]
                    }else{
                        output.df[i,7] = NA
                    }
                    
                    # Check if gene is present in stat table. If so, output info. If not, 
                    # output NAs
                    if(input_vec[i] %in% stat_table$Gene_Name){
                        output.df[i,6] = stat_table$Gene_Description[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,8] = stat_table$Pi_RioChoy[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,9] = stat_table$Pi_Pachon[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,10] = stat_table$Pi_Molino[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,11] = stat_table$Pi_Tinaja[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,12] = stat_table$Pi_Rascon[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,13] = stat_table$Dxy_RioChoy.Pachon[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,14] = stat_table$Dxy_RioChoy.Molino[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,15] = stat_table$Dxy_RioChoy.Tinaja[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,16] = stat_table$Dxy_Rascon.Pachon[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,17] = stat_table$Dxy_Rascon.Tinaja[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,18] = stat_table$Fst_RioChoy.Pachon[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,19] = stat_table$Fst_RioChoy.Molino[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,20] = stat_table$Fst_RioChoy.Tinaja[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,21] = stat_table$Fst_Pachon.Rascon[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,22] = stat_table$Fst_Rascon.Tinaja[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,23] = stat_table$TajimasD_RioChoy[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,24] = stat_table$TajimasD_Pachon[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,25] = stat_table$TajimasD_Molino[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,26] = stat_table$TajimasD_Tinaja[stat_table$Gene_Name == input_vec[i]]
                        output.df[i,27] = stat_table$TajimasD_Rascon[stat_table$Gene_Name == input_vec[i]]
                    }else{
                        output.df[i,6] = NA
                        output.df[i,8] = NA
                        output.df[i,9] = NA
                        output.df[i,10] = NA
                        output.df[i,11] = NA
                        output.df[i,12] = NA
                        output.df[i,13] = NA
                        output.df[i,14] = NA
                        output.df[i,15] = NA
                        output.df[i,16] = NA
                        output.df[i,17] = NA
                        output.df[i,18] = NA
                        output.df[i,19] = NA
                        output.df[i,20] = NA
                        output.df[i,21] = NA
                        output.df[i,22] = NA
                        output.df[i,23] = NA
                        output.df[i,24] = NA
                        output.df[i,25] = NA
                        output.df[i,26] = NA
                        output.df[i,27] = NA
                    }
                    
                    # Check if current gene is found in transcription data. If so, output 
                    # associated information. If not, output NAs
                    if(input_vec[i] %in% condition_control$Gene_name){
                        if(length(condition_control$logFC[
                            (condition_control$Gene_name == input_vec[i]) &
                            (grepl("Choy",condition_control$Class))]) != 0){
                            output.df[i,28] = condition_control$logFC[
                                (condition_control$Gene_name == input_vec[i]) &
                                    (grepl("Choy",condition_control$Class))]
                            output.df[i,32] = condition_control$PValue[
                                (condition_control$Gene_name == input_vec[i]) &
                                    (grepl("Choy",condition_control$Class))]
                        }else{
                            output.df[i,28] = NA
                            output.df[i,32] = NA
                        }
                        if(length(condition_control$logFC[
                            (condition_control$Gene_name == input_vec[i]) &
                            (grepl("Pachon",condition_control$Class))]) != 0){
                            output.df[i,29] = condition_control$logFC[
                                (condition_control$Gene_name == input_vec[i]) &
                                    (grepl("Pachon",condition_control$Class))]
                            output.df[i,33] = condition_control$PValue[
                                (condition_control$Gene_name == input_vec[i]) &
                                    (grepl("Pachon",condition_control$Class))]
                        }else{
                            output.df[i,29] = NA
                            output.df[i,33] = NA
                        }
                        if(length(condition_control$logFC[
                            (condition_control$Gene_name == input_vec[i]) &
                            (grepl("Molino",condition_control$Class))]) != 0){
                            output.df[i,30] = condition_control$logFC[
                                (condition_control$Gene_name == input_vec[i]) &
                                    (grepl("Molino",condition_control$Class))]
                            output.df[i,34] = condition_control$PValue[
                                (condition_control$Gene_name == input_vec[i]) &
                                    (grepl("Molino",condition_control$Class))]
                        }else{
                            output.df[i,30] = NA
                            output.df[i,34] = NA
                        }
                        if(length(condition_control$logFC[
                            (condition_control$Gene_name == input_vec[i]) &
                            (grepl("Tinaja",condition_control$Class))]) != 0){
                            output.df[i,31] = condition_control$logFC[
                                (condition_control$Gene_name == input_vec[i]) &
                                    (grepl("Tinaja",condition_control$Class))]
                            output.df[i,35] = condition_control$PValue[
                                (condition_control$Gene_name == input_vec[i]) &
                                    (grepl("Tinaja",condition_control$Class))]
                        }else{
                            output.df[i,31] = NA
                            output.df[i,35] = NA
                        }
                    }else{
                        output.df[i,28] = NA
                        output.df[i,29] = NA
                        output.df[i,30] = NA
                        output.df[i,31] = NA
                        output.df[i,32] = NA
                        output.df[i,33] = NA
                        output.df[i,34] = NA
                        output.df[i,35] = NA
                    }
                    
                    geneName = T
                    geneID = F
                    # If element is present in all_IDs, consider vector a vector of IDs and 
                    # output all associated elements to output dataframe
                }else if(input_vec[i] %in% all.genes_IDs$all_IDs){
                    output.df[i,1] = all.genes_IDs$all_genes[all.genes_IDs$all_IDs == input_vec[i]]
                    output.df[i,2] = input_vec[i]
                    
                    # If the current gene is present in the position table, output position 
                    # table info. If not, output all NA
                    if(output.df[i,2] %in% position_table$V16){
                        output.df[i,3] = position_table$V1[position_table$V16 == output.df[i,2]]
                        output.df[i,4] = position_table$V4[position_table$V16 == output.df[i,2]]
                        output.df[i,5] = position_table$V5[position_table$V16 == output.df[i,2]]
                    }else{
                        output.df[i,3] = NA
                        output.df[i,4] = NA
                        output.df[i,5] = NA
                    }
                    
                    # Check if gene is present in GO term table. If so, output GO terms. If 
                    # not, output NA
                    if(output.df[i,2] %in% GeneToGO$Gene.names){
                        output.df[i,7] = GeneToGO$Gene.ontology.IDs[GeneToGO$Gene.names == output.df[i,2]]
                    }else{
                        output.df[i,7] = NA
                    }
                    
                    # Check if gene is present in stat table. If so, output info. If not, 
                    # output NAs
                    if(input_vec[i] %in% stat_table$Stable_Gene_ID){
                        output.df[i,6] = stat_table$Gene_Description[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,8] = stat_table$Pi_RioChoy[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,9] = stat_table$Pi_Pachon[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,10] = stat_table$Pi_Molino[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,11] = stat_table$Pi_Tinaja[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,12] = stat_table$Pi_Rascon[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,13] = stat_table$Dxy_RioChoy.Pachon[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,14] = stat_table$Dxy_RioChoy.Molino[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,15] = stat_table$Dxy_RioChoy.Tinaja[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,16] = stat_table$Dxy_Rascon.Pachon[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,17] = stat_table$Dxy_Rascon.Tinaja[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,18] = stat_table$Fst_RioChoy.Pachon[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,19] = stat_table$Fst_RioChoy.Molino[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,20] = stat_table$Fst_RioChoy.Tinaja[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,21] = stat_table$Fst_Pachon.Rascon[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,22] = stat_table$Fst_Rascon.Tinaja[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,23] = stat_table$TajimasD_RioChoy[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,24] = stat_table$TajimasD_Pachon[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,25] = stat_table$TajimasD_Molino[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,26] = stat_table$TajimasD_Tinaja[stat_table$Stable_Gene_ID == input_vec[i]]
                        output.df[i,27] = stat_table$TajimasD_Rascon[stat_table$Stable_Gene_ID == input_vec[i]]
                    }else{
                        output.df[i,6] = NA
                        output.df[i,8] = NA
                        output.df[i,9] = NA
                        output.df[i,10] = NA
                        output.df[i,11] = NA
                        output.df[i,12] = NA
                        output.df[i,13] = NA
                        output.df[i,14] = NA
                        output.df[i,15] = NA
                        output.df[i,16] = NA
                        output.df[i,17] = NA
                        output.df[i,18] = NA
                        output.df[i,19] = NA
                        output.df[i,20] = NA
                        output.df[i,21] = NA
                        output.df[i,22] = NA
                        output.df[i,23] = NA
                        output.df[i,24] = NA
                        output.df[i,25] = NA
                        output.df[i,26] = NA
                        output.df[i,27] = NA
                    }
                    
                    # Check if current gene is found in transcription data. If so, output 
                    # associated information. If not, output NAs
                    if(input_vec[i] %in% condition_control$ï..Gene_stable_ID){
                        if(length(condition_control$logFC[
                            (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                            (grepl("Choy",condition_control$Class))]) != 0){
                            output.df[i,28] = condition_control$logFC[
                                (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                                    (grepl("Choy",condition_control$Class))]
                            output.df[i,32] = condition_control$PValue[
                                (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                                    (grepl("Choy",condition_control$Class))]
                        }else{
                            output.df[i,28] = NA
                            output.df[i,32] = NA
                        }
                        if(length(condition_control$logFC[
                            (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                            (grepl("Pachon",condition_control$Class))]) != 0){
                            output.df[i,29] = condition_control$logFC[
                                (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                                    (grepl("Pachon",condition_control$Class))]
                            output.df[i,33] = condition_control$PValue[
                                (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                                    (grepl("Pachon",condition_control$Class))]
                        }else{
                            output.df[i,29] = NA
                            output.df[i,33] = NA
                        }
                        if(length(condition_control$logFC[
                            (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                            (grepl("Molino",condition_control$Class))]) != 0){
                            output.df[i,30] = condition_control$logFC[
                                (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                                    (grepl("Molino",condition_control$Class))]
                            output.df[i,34] = condition_control$PValue[
                                (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                                    (grepl("Molino",condition_control$Class))]
                        }else{
                            output.df[i,30] = NA
                            output.df[i,34] = NA
                        }
                        if(length(condition_control$logFC[
                            (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                            (grepl("Tinaja",condition_control$Class))]) != 0){
                            output.df[i,31] = condition_control$logFC[
                                (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                                    (grepl("Tinaja",condition_control$Class))]
                            output.df[i,35] = condition_control$PValue[
                                (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                                    (grepl("Tinaja",condition_control$Class))]
                        }else{
                            output.df[i,31] = NA
                            output.df[i,35] = NA
                        }
                    }else{
                        output.df[i,28] = NA
                        output.df[i,29] = NA
                        output.df[i,30] = NA
                        output.df[i,31] = NA
                        output.df[i,32] = NA
                        output.df[i,33] = NA
                        output.df[i,34] = NA
                        output.df[i,35] = NA
                    }
                    geneName = T
                    geneID = F
                    # If element is present in neither, output an error  
                }else if(!(input_vec[i] %in% all.genes_IDs$all_genes) & 
                         !(input_vec[i] %in% all.genes_IDs$all_IDs)){
                    return(paste(c("ERROR: No transcription or statistical data present for 
                       the input", input_vec[i], "."), collapse = " "))
                }
            }
            # If input is NOT a comma-separated string, check whether input is a ID, name,
            # or phrase
        }else if(!grepl(comma, input)){
            if(input %in% all.genes_IDs$all_genes){
                # Store input and associated values in appropriate objects and mark input
                # as a gene name
                output.df <- data.frame(matrix(ncol = 35, nrow = 1))
                output.df[1,1] = input
                output.df[1,2] = all.genes_IDs$all_IDs[all.genes_IDs$all_genes == input]
                
                # If the current gene is present in the position table, output position 
                # table info. If not, output all NA
                if(input %in% position_table$V16){
                    output.df[1,3] = position_table$V1[position_table$V16 == input]
                    output.df[1,4] = position_table$V4[position_table$V16 == input]
                    output.df[1,5] = position_table$V5[position_table$V16 == input]
                }else{
                    output.df[1,3] = NA
                    output.df[1,4] = NA
                    output.df[1,5] = NA
                }
                
                # Check if gene is present in GO term table. If so, output GO terms. If 
                # not, output NA
                if(input %in% GeneToGO$Gene.names){
                    output.df[1,7] = GeneToGO$Gene.ontology.IDs[GeneToGO$Gene.names == input]
                }else{
                    output.df[1,7] = NA
                }
                
                # Check if gene is present in stat table. If so, output info. If not, 
                # output NAs
                if(input %in% stat_table$Gene_Name){
                    output.df[1,6] = stat_table$Gene_Description[stat_table$Gene_Name == input]
                    output.df[1,8] = stat_table$Pi_RioChoy[stat_table$Gene_Name == input]
                    output.df[1,9] = stat_table$Pi_Pachon[stat_table$Gene_Name == input]
                    output.df[1,10] = stat_table$Pi_Molino[stat_table$Gene_Name == input]
                    output.df[1,11] = stat_table$Pi_Tinaja[stat_table$Gene_Name == input]
                    output.df[1,12] = stat_table$Pi_Rascon[stat_table$Gene_Name == input]
                    output.df[1,13] = stat_table$Dxy_RioChoy.Pachon[stat_table$Gene_Name == input]
                    output.df[1,14] = stat_table$Dxy_RioChoy.Molino[stat_table$Gene_Name == input]
                    output.df[1,15] = stat_table$Dxy_RioChoy.Tinaja[stat_table$Gene_Name == input]
                    output.df[1,16] = stat_table$Dxy_Rascon.Pachon[stat_table$Gene_Name == input]
                    output.df[1,17] = stat_table$Dxy_Rascon.Tinaja[stat_table$Gene_Name == input]
                    output.df[1,18] = stat_table$Fst_RioChoy.Pachon[stat_table$Gene_Name == input]
                    output.df[1,19] = stat_table$Fst_RioChoy.Molino[stat_table$Gene_Name == input]
                    output.df[1,20] = stat_table$Fst_RioChoy.Tinaja[stat_table$Gene_Name == input]
                    output.df[1,21] = stat_table$Fst_Pachon.Rascon[stat_table$Gene_Name == input]
                    output.df[1,22] = stat_table$Fst_Rascon.Tinaja[stat_table$Gene_Name == input]
                    output.df[1,23] = stat_table$TajimasD_RioChoy[stat_table$Gene_Name == input]
                    output.df[1,24] = stat_table$TajimasD_Pachon[stat_table$Gene_Name == input]
                    output.df[1,25] = stat_table$TajimasD_Molino[stat_table$Gene_Name == input]
                    output.df[1,26] = stat_table$TajimasD_Tinaja[stat_table$Gene_Name == input]
                    output.df[1,27] = stat_table$TajimasD_Rascon[stat_table$Gene_Name == input]
                }else{
                    output.df[1,6] = NA
                    output.df[1,8] = NA
                    output.df[1,9] = NA
                    output.df[1,10] = NA
                    output.df[1,11] = NA
                    output.df[1,12] = NA
                    output.df[1,13] = NA
                    output.df[1,14] = NA
                    output.df[1,15] = NA
                    output.df[1,16] = NA
                    output.df[1,17] = NA
                    output.df[1,18] = NA
                    output.df[1,19] = NA
                    output.df[1,20] = NA
                    output.df[1,21] = NA
                    output.df[1,22] = NA
                    output.df[1,23] = NA
                    output.df[1,24] = NA
                    output.df[1,25] = NA
                    output.df[1,26] = NA
                    output.df[1,27] = NA
                }
                
                # Check if current gene is found in transcription data. If so, output 
                # associated information. If not, output NAs
                if(input %in% condition_control$Gene_name){
                    if(length(condition_control$logFC[
                        (condition_control$Gene_name == input_vec[i]) &
                        (grepl("Choy",condition_control$Class))]) != 0){
                        output.df[1,28] = condition_control$logFC[
                            (condition_control$Gene_name == input_vec[i]) &
                                (grepl("Choy",condition_control$Class))]
                        output.df[1,32] = condition_control$PValue[
                            (condition_control$Gene_name == input_vec[i]) &
                                (grepl("Choy",condition_control$Class))]
                    }else{
                        output.df[1,28] = NA
                        output.df[1,32] = NA
                    }
                    if(length(condition_control$logFC[
                        (condition_control$Gene_name == input_vec[i]) &
                        (grepl("Pachon",condition_control$Class))]) != 0){
                        output.df[1,29] = condition_control$logFC[
                            (condition_control$Gene_name == input_vec[i]) &
                                (grepl("Pachon",condition_control$Class))]
                        output.df[1,33] = condition_control$PValue[
                            (condition_control$Gene_name == input_vec[i]) &
                                (grepl("Pachon",condition_control$Class))]
                    }else{
                        output.df[1,29] = NA
                        output.df[1,33] = NA
                    }
                    if(length(condition_control$logFC[
                        (condition_control$Gene_name == input_vec[i]) &
                        (grepl("Molino",condition_control$Class))]) != 0){
                        output.df[1,30] = condition_control$logFC[
                            (condition_control$Gene_name == input_vec[i]) &
                                (grepl("Molino",condition_control$Class))]
                        output.df[1,34] = condition_control$PValue[
                            (condition_control$Gene_name == input_vec[i]) &
                                (grepl("Molino",condition_control$Class))]
                    }else{
                        output.df[1,30] = NA
                        output.df[1,34] = NA
                    }
                    if(length(condition_control$logFC[
                        (condition_control$Gene_name == input_vec[i]) &
                        (grepl("Tinaja",condition_control$Class))]) != 0){
                        output.df[1,31] = condition_control$logFC[
                            (condition_control$Gene_name == input_vec[i]) &
                                (grepl("Tinaja",condition_control$Class))]
                        output.df[1,35] = condition_control$PValue[
                            (condition_control$Gene_name == input_vec[i]) &
                                (grepl("Tinaja",condition_control$Class))]
                    }else{
                        output.df[1,31] = NA
                        output.df[1,35] = NA
                    }
                }else{
                    output.df[1,28] = NA
                    output.df[1,29] = NA
                    output.df[1,30] = NA
                    output.df[1,31] = NA
                    output.df[1,32] = NA
                    output.df[1,33] = NA
                    output.df[1,34] = NA
                    output.df[1,35] = NA
                }
                
                geneName = T
                geneID = F
            }else if(input %in% all.genes_IDs$all_IDs){
                # Store input and associated values in appropriate objects and mark input
                # as a gene ID
                output.df <- data.frame(matrix(ncol = 35, nrow = 1))
                output.df[1,1] = all.genes_IDs$all_genes[all.genes_IDs$all_IDs == input]
                output.df[1,2] = input
                
                # If the current gene is present in the position table, output position 
                # table info. If not, output all NA
                if(output.df[1,2] %in% position_table$V16){
                    output.df[1,3] = position_table$V1[position_table$V16 == output.df[1,2]]
                    output.df[1,4] = position_table$V4[position_table$V16 == output.df[1,2]]
                    output.df[1,5] = position_table$V5[position_table$V16 == output.df[1,2]]
                }else{
                    output.df[1,3] = NA
                    output.df[1,4] = NA
                    output.df[1,5] = NA
                }
                
                # Check if gene is present in GO term table. If so, output GO terms. If 
                # not, output NA
                if(output.df[1,2] %in% GeneToGO$Gene.names){
                    output.df[1,7] = GeneToGO$Gene.ontology.IDs[GeneToGO$Gene.names == output.df[1,2]]
                }else{
                    output.df[1,7] = NA
                }
                
                # Check if gene is present in stat table. If so, output info. If not, 
                # output NAs
                if(input %in% stat_table$Stable_Gene_ID){
                    output.df[1,6] = stat_table$Gene_Description[stat_table$Stable_Gene_ID == input]
                    output.df[1,8] = stat_table$Pi_RioChoy[stat_table$Stable_Gene_ID == input]
                    output.df[1,9] = stat_table$Pi_Pachon[stat_table$Stable_Gene_ID == input]
                    output.df[1,10] = stat_table$Pi_Molino[stat_table$Stable_Gene_ID == input]
                    output.df[1,11] = stat_table$Pi_Tinaja[stat_table$Stable_Gene_ID == input]
                    output.df[1,12] = stat_table$Pi_Rascon[stat_table$Stable_Gene_ID == input]
                    output.df[1,13] = stat_table$Dxy_RioChoy.Pachon[stat_table$Stable_Gene_ID == input]
                    output.df[1,14] = stat_table$Dxy_RioChoy.Molino[stat_table$Stable_Gene_ID == input]
                    output.df[1,15] = stat_table$Dxy_RioChoy.Tinaja[stat_table$Stable_Gene_ID == input]
                    output.df[1,16] = stat_table$Dxy_Rascon.Pachon[stat_table$Stable_Gene_ID == input]
                    output.df[1,17] = stat_table$Dxy_Rascon.Tinaja[stat_table$Stable_Gene_ID == input]
                    output.df[1,18] = stat_table$Fst_RioChoy.Pachon[stat_table$Stable_Gene_ID == input]
                    output.df[1,19] = stat_table$Fst_RioChoy.Molino[stat_table$Stable_Gene_ID == input]
                    output.df[1,20] = stat_table$Fst_RioChoy.Tinaja[stat_table$Stable_Gene_ID == input]
                    output.df[1,21] = stat_table$Fst_Pachon.Rascon[stat_table$Stable_Gene_ID == input]
                    output.df[1,22] = stat_table$Fst_Rascon.Tinaja[stat_table$Stable_Gene_ID == input]
                    output.df[1,23] = stat_table$TajimasD_RioChoy[stat_table$Stable_Gene_ID == input]
                    output.df[1,24] = stat_table$TajimasD_Pachon[stat_table$Stable_Gene_ID == input]
                    output.df[1,25] = stat_table$TajimasD_Molino[stat_table$Stable_Gene_ID == input]
                    output.df[1,26] = stat_table$TajimasD_Tinaja[stat_table$Stable_Gene_ID == input]
                    output.df[1,27] = stat_table$TajimasD_Rascon[stat_table$Stable_Gene_ID == input]
                }else{
                    output.df[1,6] = NA
                    output.df[1,8] = NA
                    output.df[1,9] = NA
                    output.df[1,10] = NA
                    output.df[1,11] = NA
                    output.df[1,12] = NA
                    output.df[1,13] = NA
                    output.df[1,14] = NA
                    output.df[1,15] = NA
                    output.df[1,16] = NA
                    output.df[1,17] = NA
                    output.df[1,18] = NA
                    output.df[1,19] = NA
                    output.df[1,20] = NA
                    output.df[1,21] = NA
                    output.df[1,22] = NA
                    output.df[1,23] = NA
                    output.df[1,24] = NA
                    output.df[1,25] = NA
                    output.df[1,26] = NA
                    output.df[1,27] = NA
                }
                
                # Check if current gene is found in transcription data. If so, output 
                # associated information. If not, output NAs
                if(input %in% condition_control$ï..Gene_stable_ID){
                    if(length(condition_control$logFC[
                        (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                        (grepl("Choy",condition_control$Class))]) != 0){
                        output.df[1,28] = condition_control$logFC[
                            (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                                (grepl("Choy",condition_control$Class))]
                        output.df[1,32] = condition_control$PValue[
                            (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                                (grepl("Choy",condition_control$Class))]
                    }else{
                        output.df[1,28] = NA
                        output.df[1,32] = NA
                    }
                    if(length(condition_control$logFC[
                        (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                        (grepl("Pachon",condition_control$Class))]) != 0){
                        output.df[1,29] = condition_control$logFC[
                            (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                                (grepl("Pachon",condition_control$Class))]
                        output.df[1,33] = condition_control$PValue[
                            (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                                (grepl("Pachon",condition_control$Class))]
                    }else{
                        output.df[1,29] = NA
                        output.df[1,33] = NA
                    }
                    if(length(condition_control$logFC[
                        (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                        (grepl("Molino",condition_control$Class))]) != 0){
                        output.df[1,30] = condition_control$logFC[
                            (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                                (grepl("Molino",condition_control$Class))]
                        output.df[1,34] = condition_control$PValue[
                            (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                                (grepl("Molino",condition_control$Class))]
                    }else{
                        output.df[1,30] = NA
                        output.df[1,34] = NA
                    }
                    if(length(condition_control$logFC[
                        (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                        (grepl("Tinaja",condition_control$Class))]) != 0){
                        output.df[1,31] = condition_control$logFC[
                            (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                                (grepl("Tinaja",condition_control$Class))]
                        output.df[1,35] = condition_control$PValue[
                            (condition_control$ï..Gene_stable_ID == input_vec[i]) &
                                (grepl("Tinaja",condition_control$Class))]
                    }else{
                        output.df[1,31] = NA
                        output.df[1,35] = NA
                    }
                }else{
                    output.df[1,28] = NA
                    output.df[1,29] = NA
                    output.df[1,30] = NA
                    output.df[1,31] = NA
                    output.df[1,32] = NA
                    output.df[1,33] = NA
                    output.df[1,34] = NA
                    output.df[1,35] = NA
                }
                geneID = T
                geneName = F
            }else if(!(input %in% all.genes_IDs$all_genes) & 
                     !(input %in% all.genes_IDs$all_IDs)){
                geneName = F
                geneID = F
            }
        } 
        
        
        # If the input string contains NO commas but is not WHOLLY comprised of IDs or
        # names, consider the string a phrase
        if((!geneID) & (!geneName)){
            # Find all gene names whose description, GO terms, or gene name contains the
            # phrase-of-interest
            
            # Must first check if phrase is found in any column-of-interest prior to 
            # appending phrase to column of interest
            input_vec <- c()
            if((sum(grepl(input, stat_table$Gene_Description)) != 0) |
               (sum(grepl(input, stat_table$GO_Terms)) != 0)){
                input_vec <- append(input_vec, stat_table$Gene_Name[
                    (grepl(input, stat_table$Gene_Description)) |
                        (grepl(input, stat_table$GO_Terms))]
                )
            }
            if((sum(grepl(input, GeneToGO$Gene.ontology..biological.process.)) != 0) |
               (sum(grepl(input, GeneToGO$Gene.ontology..cellular.component.)) != 0) |
               (sum(grepl(input, GeneToGO$Gene.ontology..molecular.function.)) != 0)){
                input_vec <- append(input_vec, GeneToGO$Gene.names[
                    grepl(input, GeneToGO$Gene.ontology..biological.process.),
                    grepl(input, GeneToGO$Gene.ontology..cellular.component.),
                    grepl(input, GeneToGO$Gene.ontology..molecular.function.)
                ]
                )
            }
            if(sum(grepl(input, all.genes_IDs$all_genes)) != 0){
                input_vec <- append(input_vec, 
                                    all.genes_IDs$all_genes[grepl(input, 
                                                                  all.genes_IDs$all_genes)])
            }
            
            # Remove duplicate genes from input vector
            input_vec <- input_vec[!duplicated(input_vec)]
            
            # If no gene descriptions contain the phrase, return an error
            if(length(input_vec) == 0){
                return(paste(c("ERROR: No genes-of-interest can be described by the phrase",
                               input), collapse = " "))
            }
            # Initialize df in which to store output based on number of genes
            output.df <- data.frame(matrix(nrow = length(input_vec), ncol = 35))
            
            # For each stable ID, collect all values associated with the stable ID
            for(i in 1:length(input_vec)){
                output.df[i,1] = input_vec[i]
                output.df[i,2] = all.genes_IDs$all_IDs[all.genes_IDs$all_genes == input_vec[i]]
                
                # If the current gene is present in the position table, output position 
                # table info. If not, output all NA
                if(input_vec[i] %in% position_table$V16){
                    output.df[i,3] = position_table$V1[position_table$V16 == input_vec[i]]
                    output.df[i,4] = position_table$V4[position_table$V16 == input_vec[i]]
                    output.df[i,5] = position_table$V5[position_table$V16 == input_vec[i]]
                }else{
                    output.df[i,3] = NA
                    output.df[i,4] = NA
                    output.df[i,5] = NA
                }
                
                # Check if gene is present in GO term table. If so, output GO terms. If 
                # not, output NA
                if(input_vec[i] %in% GeneToGO$Gene.names){
                    output.df[i,7] = GeneToGO$Gene.ontology.IDs[GeneToGO$Gene.names == input_vec[i]]
                }else{
                    output.df[i,7] = NA
                }
                
                # Check if gene is present in stat table. If so, output info. If not, 
                # output NAs
                if(input_vec[i] %in% stat_table$Gene_Name){
                    output.df[i,6] = stat_table$Gene_Description[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,8] = stat_table$Pi_RioChoy[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,9] = stat_table$Pi_Pachon[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,10] = stat_table$Pi_Molino[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,11] = stat_table$Pi_Tinaja[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,12] = stat_table$Pi_Rascon[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,13] = stat_table$Dxy_RioChoy.Pachon[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,14] = stat_table$Dxy_RioChoy.Molino[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,15] = stat_table$Dxy_RioChoy.Tinaja[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,16] = stat_table$Dxy_Rascon.Pachon[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,17] = stat_table$Dxy_Rascon.Tinaja[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,18] = stat_table$Fst_RioChoy.Pachon[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,19] = stat_table$Fst_RioChoy.Molino[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,20] = stat_table$Fst_RioChoy.Tinaja[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,21] = stat_table$Fst_Pachon.Rascon[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,22] = stat_table$Fst_Rascon.Tinaja[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,23] = stat_table$TajimasD_RioChoy[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,24] = stat_table$TajimasD_Pachon[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,25] = stat_table$TajimasD_Molino[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,26] = stat_table$TajimasD_Tinaja[stat_table$Gene_Name == input_vec[i]]
                    output.df[i,27] = stat_table$TajimasD_Rascon[stat_table$Gene_Name == input_vec[i]]
                }else{
                    output.df[i,6] = NA
                    output.df[i,8] = NA
                    output.df[i,9] = NA
                    output.df[i,10] = NA
                    output.df[i,11] = NA
                    output.df[i,12] = NA
                    output.df[i,13] = NA
                    output.df[i,14] = NA
                    output.df[i,15] = NA
                    output.df[i,16] = NA
                    output.df[i,17] = NA
                    output.df[i,18] = NA
                    output.df[i,19] = NA
                    output.df[i,20] = NA
                    output.df[i,21] = NA
                    output.df[i,22] = NA
                    output.df[i,23] = NA
                    output.df[i,24] = NA
                    output.df[i,25] = NA
                    output.df[i,26] = NA
                    output.df[i,27] = NA
                }
                
                # Check if current gene is found in transcription data. If so, output 
                # associated information. If not, output NAs
                if(input_vec[i] %in% condition_control$Gene_name){
                    if(length(condition_control$logFC[
                        (condition_control$Gene_name == input_vec[i]) &
                        (grepl("Choy",condition_control$Class))]) != 0){
                        output.df[i,28] = condition_control$logFC[
                            (condition_control$Gene_name == input_vec[i]) &
                                (grepl("Choy",condition_control$Class))]
                        output.df[i,32] = condition_control$PValue[
                            (condition_control$Gene_name == input_vec[i]) &
                                (grepl("Choy",condition_control$Class))]
                    }else{
                        output.df[i,28] = NA
                        output.df[i,32] = NA
                    }
                    if(length(condition_control$logFC[
                        (condition_control$Gene_name == input_vec[i]) &
                        (grepl("Pachon",condition_control$Class))]) != 0){
                        output.df[i,29] = condition_control$logFC[
                            (condition_control$Gene_name == input_vec[i]) &
                                (grepl("Pachon",condition_control$Class))]
                        output.df[i,33] = condition_control$PValue[
                            (condition_control$Gene_name == input_vec[i]) &
                                (grepl("Pachon",condition_control$Class))]
                    }else{
                        output.df[i,29] = NA
                        output.df[i,33] = NA
                    }
                    if(length(condition_control$logFC[
                        (condition_control$Gene_name == input_vec[i]) &
                        (grepl("Molino",condition_control$Class))]) != 0){
                        output.df[i,30] = condition_control$logFC[
                            (condition_control$Gene_name == input_vec[i]) &
                                (grepl("Molino",condition_control$Class))]
                        output.df[i,34] = condition_control$PValue[
                            (condition_control$Gene_name == input_vec[i]) &
                                (grepl("Molino",condition_control$Class))]
                    }else{
                        output.df[i,30] = NA
                        output.df[i,34] = NA
                    }
                    if(length(condition_control$logFC[
                        (condition_control$Gene_name == input_vec[i]) &
                        (grepl("Tinaja",condition_control$Class))]) != 0){
                        output.df[i,31] = condition_control$logFC[
                            (condition_control$Gene_name == input_vec[i]) &
                                (grepl("Tinaja",condition_control$Class))]
                        output.df[i,35] = condition_control$PValue[
                            (condition_control$Gene_name == input_vec[i]) &
                                (grepl("Tinaja",condition_control$Class))]
                    }else{
                        output.df[i,31] = NA
                        output.df[i,35] = NA
                    }
                }else{
                    output.df[i,28] = NA
                    output.df[i,29] = NA
                    output.df[i,30] = NA
                    output.df[i,31] = NA
                    output.df[i,32] = NA
                    output.df[i,33] = NA
                    output.df[i,34] = NA
                    output.df[i,35] = NA
                }
            }
        }
        
        # Output all values obtained for gene(s) of interest
        names(output.df) <- c(
            "Gene Name",
            "Gene Stable ID",
            "Scaffold",
            "Start Position",
            "Stop Position",
            "Gene Description",
            "GO Term(s)",
            "Pi_RioChoy",
            "Pi_Pachon",
            "Pi_Molino",
            "Pi_Tinaja",
            "Pi_Rascon",
            "Dxy_RioChoy:Pachon",
            "Dxy_RioChoy:Molino",
            "Dxy_RioChoy:Tinaja",
            "Dxy_Rascon:Pachon",
            "Dxy_Rascon:Tinaja",
            "Fst_RioChoy:Pachon",
            "Fst_RioChoy:Molino",
            "Fst_RioChoy:Tinaja",
            "Fst_Pachon:Rascon",
            "Fst_Rascon:Tinaja",
            "TajimasD_RioChoy",
            "TajimasD_Pachon",
            "TajimasD_Molino",
            "TajimasD_Tinaja",
            "TajimasD_Rascon",
            "SD_log(FC)_RioChoy",
            "SD_log(FC)_Pachon",
            "SD_log(FC)_Molino",
            "SD_log(FC)_Tinaja",
            "p-value for SD_log(FC)_RioChoy",
            "p-value for SD_log(FC)_Pachon",
            "p-value for SD_log(FC)_Molino",
            "p-value for SD_log(FC)_Tinaja"
        )
        return(output.df)
    }
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
                GeneCentOutput()
        }else if(typeof(GeneCentOutput()) == "character"){
            data.frame()
        }
    })
    output$GeneCent_warnings <- renderText({
        if(typeof(GeneCentOutput()) == "list"){
            "No warnings or errors"
        }else if(typeof(GeneCentOutput()) == "character"){
            GeneCentOutput()
        }
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
