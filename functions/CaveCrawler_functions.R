########################### Load Required Libraries ############################

library(shinyWidgets)
library(shiny)
library(ggplot2)
library(plotly)
library(WVPlots)
library(stringr)
library(tibble)
library(gridExtra)
library(dplyr)

################################## Load Data ###################################
# Make dataframe of latitudes and longitudes for all morphs
# Citation:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3282648/pdf/1471-2148-12-9.pdf
Latit_Longit_unedited <- data.frame(
  Population = c("Pachón","Yerbaniz","Japonés","Arroyo","Tinaja","Curva","Toro",
                 "Chica","Molino","Caballo Moro","Subterráneo","Río Frío",
                 "Arroyo Sarco", "Chamal","Río Meco","Río Tantáon","Río Florído",
                 "RióTampaón","Nacimiento del Río Santa Clara",
                 "San Rafael Los Castros","Rio Subterráneo Valley"),
  Latitude = c(22.60,22.20,22.10,22.20,22.08,21.98,21.85,21.85,23.06,22.92,
               22.10,22.99,22.02,22.84,22.82,22.37,21.98,21.85,22.50,22.75,22.13
  ),
  Longitude = c(-99.05,-98.97,-98.95,-98.97,-98.95,-98.93,-98.93,-98.93,-99.16,
                -99.20,-99.18,-99.15,-99.32,-99.20,-99.31,-98.90,-98.77,-98.94,-98.9,
                -99.02,-99.17
  )
)
# Rascon and Rio Choy data obtained from google maps
Latit_Longit <- rbind(Latit_Longit_unedited, data.frame(
  Population = c("Rascón","Río Choy"),
  Latitude = c(21.9750,21.9998),
  Longitude = c(-99.2578,-98.7785)
))
# Edit data to include only the populations currently available on the website
Latit_Longit <- Latit_Longit[
  Latit_Longit$Population %in% 
    c("Chica","Molino","Tinaja","Pachón","Rascón","Río Choy"),]

# Get a map of the world
world_map_1 <- map_data("world")
# get map of USA, Mexico, Belize, and Guatemala
world_map <- subset(world_map_1, region %in% c("USA", "Mexico", "Belize","Guatemala"))

position_table <- read.csv("data/AmexPositionTable.csv", fill = TRUE)
# Remove duplicates so multiple rows with the same values are not outputted
position_table$Gene_Name <- tolower(position_table$Gene_Name)
position_table <- position_table[!duplicated(position_table$Gene_Name),]

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
chica_table$Publication_Name <- rep("Chica paper (Provisional)", 
                                    nrow(chica_table))
names(chica_table) <- names(outlier_table)
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

################################## Functions ###################################
# Gene Search Page: Input a single or comma-separated list of genes or gene IDs
# or a phrase associated with a gene-of-interest and output all available data
GeneCentered <- function(input, stat_table, GeneToGO, condition_control,
                         position_table){
  comma <- ", "
  
  # Dataframe in which to store inputs and associated values
  output.df <- data.frame(matrix(ncol = 38))
  
  # If input is a comma-separated string, parse string and separate elements,
  # then determine whether vector contains gene names or gene IDs
  if(grepl(comma, input)){
    input_vec <- str_split(input, pattern = comma)[[1]]
    
    # Next, iterate through each element of vector
    for(i in 1:length(input_vec)){
      # If element is present in all_genes, consider vector a vector of gene
      # names and output all associated elements to output dataframe
      if(input_vec[i] %in% c(position_table$Gene_Name, stat_table$Gene_Name,
                             condition_control$Gene_name)){
        # Find number of times this gene occurs in stat_table
        num_stat_copies <- sum(str_count(stat_table$Gene_Name, input_vec[i]))
        # If number of copies of gene in statistic data is 0, copy the gene at
        # least once so at least one copy is present for transcription data
        if(num_stat_copies == 0){
          num_stat_copies = 1
        }
        # Create temporary dataframe to be appended to final
        temp.df <- data.frame(matrix(nrow = num_stat_copies, ncol = 38))
        
        temp.df[,1] = rep(input_vec[i], num_stat_copies)
        temp.df[,2] = rep(all.genes_IDs$all_IDs[
          all.genes_IDs$all_genes == input_vec[i]], num_stat_copies)
        
        # If the current gene is present in the position table, output position
        # table info. If not, output all NA
        if(input_vec[i] %in% position_table$Gene_Name){
          temp.df[,3] = rep(position_table$Scaffold[
            position_table$Gene_Name == input_vec[i]], num_stat_copies)
          temp.df[,4] = rep(position_table$Start_Locus[
            position_table$Gene_Name == input_vec[i]], num_stat_copies)
          temp.df[,5] = rep(position_table$End_Locus[
            position_table$Gene_Name == input_vec[i]], num_stat_copies)
        }else{
          temp.df[,3:5] = rep(NA, num_stat_copies)
        }
        
        # Check if gene is present in GO term table. If so, output GO terms. If
        # not, output NA
        if(input_vec[i] %in% GeneToGO$Gene.names){
          temp.df[,7] = rep(GeneToGO$Gene.ontology.IDs[
            GeneToGO$Gene.names == input_vec[i]], num_stat_copies)
        }else{
          temp.df[,7] = rep(NA, num_stat_copies)
        }
        
        # Check if gene is present in stat table. If so, output info. If not,
        # output NAs
        if(input_vec[i] %in% stat_table$Gene_Name){
          temp.df[,6] = stat_table$Gene_Description[stat_table$Gene_Name == input_vec[i]]
          
          temp.df[,8] = stat_table$Pi_RioChoy[stat_table$Gene_Name == input_vec[i]]
          temp.df[,9] = stat_table$Pi_Pachon[stat_table$Gene_Name == input_vec[i]]
          temp.df[,10] = stat_table$Pi_Molino[stat_table$Gene_Name == input_vec[i]]
          temp.df[,11] = stat_table$Pi_Tinaja[stat_table$Gene_Name == input_vec[i]]
          temp.df[,12] = stat_table$Pi_Rascon[stat_table$Gene_Name == input_vec[i]]
          
          temp.df[,13] = stat_table$Dxy_RioChoy.Pachon[stat_table$Gene_Name == input_vec[i]]
          temp.df[,14] = stat_table$Dxy_RioChoy.Molino[stat_table$Gene_Name == input_vec[i]]
          temp.df[,15] = stat_table$Dxy_RioChoy.Tinaja[stat_table$Gene_Name == input_vec[i]]
          temp.df[,16] = stat_table$Dxy_Rascon.Pachon[stat_table$Gene_Name == input_vec[i]]
          temp.df[,17] = stat_table$Dxy_Rascon.Tinaja[stat_table$Gene_Name == input_vec[i]]
          temp.df[,18] = stat_table$Dxy_Chica1.Chica2[stat_table$Gene_Name == input_vec[i]]
          
          temp.df[,19] = stat_table$Fst_RioChoy.Pachon[stat_table$Gene_Name == input_vec[i]]
          temp.df[,20] = stat_table$Fst_RioChoy.Molino[stat_table$Gene_Name == input_vec[i]]
          temp.df[,21] = stat_table$Fst_RioChoy.Tinaja[stat_table$Gene_Name == input_vec[i]]
          temp.df[,22] = stat_table$Fst_Pachon.Rascon[stat_table$Gene_Name == input_vec[i]]
          temp.df[,23] = stat_table$Fst_Rascon.Tinaja[stat_table$Gene_Name == input_vec[i]]
          temp.df[,24] = stat_table$Fst_Chica1.Chica2[stat_table$Gene_Name == input_vec[i]]
          
          temp.df[,25] = stat_table$TajimasD_RioChoy[stat_table$Gene_Name == input_vec[i]]
          temp.df[,26] = stat_table$TajimasD_Pachon[stat_table$Gene_Name == input_vec[i]]
          temp.df[,27] = stat_table$TajimasD_Molino[stat_table$Gene_Name == input_vec[i]]
          temp.df[,28] = stat_table$TajimasD_Tinaja[stat_table$Gene_Name == input_vec[i]]
          temp.df[,29] = stat_table$TajimasD_Rascon[stat_table$Gene_Name == input_vec[i]]
          temp.df[,38] = stat_table$Publication_Name[stat_table$Gene_Name == input_vec[i]]
        }else{
          temp.df[,6] = rep(NA, num_stat_copies)
          temp.df[,8:29] = rep(NA, num_stat_copies)
          temp.df[,38] = rep(NA, num_stat_copies)
        }
        
        # Check if current gene is found in transcription data. If so, output
        # associated information. If not, output NAs
        if(input_vec[i] %in% condition_control$Gene_name){
          if(length(condition_control$logFC[
            (condition_control$Gene_name == input_vec[i]) &
            (grepl("Choy",condition_control$Class))]) != 0){
            temp.df[,30] = rep(condition_control$logFC[
              (condition_control$Gene_name == input_vec[i]) &
                (grepl("Choy",condition_control$Class))], num_stat_copies)
            temp.df[,34] = rep(condition_control$PValue[
              (condition_control$Gene_name == input_vec[i]) &
                (grepl("Choy",condition_control$Class))], num_stat_copies)
          }else{
            temp.df[,30] = rep(NA, num_stat_copies)
            temp.df[,34] = rep(NA, num_stat_copies)
          }
          if(length(condition_control$logFC[
            (condition_control$Gene_name == input_vec[i]) &
            (grepl("Pachon",condition_control$Class))]) != 0){
            temp.df[,31] = rep(condition_control$logFC[
              (condition_control$Gene_name == input_vec[i]) &
                (grepl("Pachon",condition_control$Class))], num_stat_copies)
            temp.df[,35] = rep(condition_control$PValue[
              (condition_control$Gene_name == input_vec[i]) &
                (grepl("Pachon",condition_control$Class))], num_stat_copies)
          }else{
            temp.df[,31] = rep(NA, num_stat_copies)
            temp.df[,35] = rep(NA, num_stat_copies)
          }
          if(length(condition_control$logFC[
            (condition_control$Gene_name == input_vec[i]) &
            (grepl("Molino",condition_control$Class))]) != 0){
            temp.df[,32] = rep(condition_control$logFC[
              (condition_control$Gene_name == input_vec[i]) &
                (grepl("Molino",condition_control$Class))], num_stat_copies)
            temp.df[,36] = rep(condition_control$PValue[
              (condition_control$Gene_name == input_vec[i]) &
                (grepl("Molino",condition_control$Class))], num_stat_copies)
          }else{
            temp.df[,32] = rep(NA, num_stat_copies)
            temp.df[,36] = rep(NA, num_stat_copies)
          }
          if(length(condition_control$logFC[
            (condition_control$Gene_name == input_vec[i]) &
            (grepl("Tinaja",condition_control$Class))]) != 0){
            temp.df[,33] = rep(condition_control$logFC[
              (condition_control$Gene_name == input_vec[i]) &
                (grepl("Tinaja",condition_control$Class))], num_stat_copies)
            temp.df[,37] = rep(condition_control$PValue[
              (condition_control$Gene_name == input_vec[i]) &
                (grepl("Tinaja",condition_control$Class))], num_stat_copies)
          }else{
            temp.df[,33] = rep(NA, num_stat_copies)
            temp.df[,37] = rep(NA, num_stat_copies)
          }
        }else{
          temp.df[,30:37] = rep(NA, num_stat_copies)
        }
        
        geneName = T
        geneID = F
        
        # If element is present in all_IDs, consider vector a vector of IDs and
        # output all associated elements to output dataframe
      }else if(input_vec[i] %in% c(position_table$Gene_ID, 
                                   stat_table$Stable_Gene_ID, 
                                   condition_control$Gene_stable_ID)){
        # Find number of times this gene ID occurs in the statistics data
        num_stat_copies <- sum(str_count(stat_table$Stable_Gene_ID, 
                                         input_vec[i]))
        # If number of copies of gene in statistic data is 0, copy the gene at
        # least once so at least one copy is present for transcription data
        if(num_stat_copies == 0){
          num_stat_copies = 1
        }
        # Create temporary dataframe to be appended to final
        temp.df <- data.frame(matrix(nrow = num_stat_copies, ncol = 38))
        
        temp.df[,1] = rep(all.genes_IDs$all_genes[
          all.genes_IDs$all_IDs == input_vec[i]], num_stat_copies)
        temp.df[,2] = rep(input_vec[i], num_stat_copies)
        
        # If the current gene is present in the position table, output position
        # table info. If not, output all NA
        if(temp.df[i,2] %in% position_table$Gene_Name){
          temp.df[,3] = rep(position_table$Scaffold[
            position_table$Gene_Name == temp.df[,2]], num_stat_copies)
          temp.df[,4] = rep(position_table$Start_Locus[
            position_table$Gene_Name == temp.df[,2]], num_stat_copies)
          temp.df[,5] = rep(position_table$End_Locus[
            position_table$Gene_Name == temp.df[,2]], num_stat_copies)
        }else{
          temp.df[,3:5] = rep(NA, num_stat_copies)
        }
        
        # Check if gene is present in GO term table. If so, output GO terms. If
        # not, output NA
        if(temp.df[i,1] %in% GeneToGO$Gene.names){
          temp.df[,7] = rep(GeneToGO$Gene.ontology.IDs[
            GeneToGO$Gene.names == temp.df[i,1]], num_stat_copies)
        }else{
          temp.df[,7] = rep(NA, num_stat_copies)
        }
        
        # Check if gene is present in stat table. If so, output info. If not,
        # output NAs
        if(input_vec[i] %in% stat_table$Stable_Gene_ID){
          temp.df[,6] = stat_table$Gene_Description[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,8] = stat_table$Pi_RioChoy[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,9] = stat_table$Pi_Pachon[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,10] = stat_table$Pi_Molino[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,11] = stat_table$Pi_Tinaja[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,12] = stat_table$Pi_Rascon[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,13] = stat_table$Dxy_RioChoy.Pachon[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,14] = stat_table$Dxy_RioChoy.Molino[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,15] = stat_table$Dxy_RioChoy.Tinaja[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,16] = stat_table$Dxy_Rascon.Pachon[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,17] = stat_table$Dxy_Rascon.Tinaja[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,18] = stat_table$Dxy_Chica1.Chica2[stat_table$Stable_Gene_ID == input_vec[i]]
          
          temp.df[,19] = stat_table$Fst_RioChoy.Pachon[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,20] = stat_table$Fst_RioChoy.Molino[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,21] = stat_table$Fst_RioChoy.Tinaja[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,22] = stat_table$Fst_Pachon.Rascon[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,23] = stat_table$Fst_Rascon.Tinaja[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,24] = stat_table$Fst_Chica1.Chica2[stat_table$Stable_Gene_ID == input_vec[i]]
          
          temp.df[,25] = stat_table$TajimasD_RioChoy[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,26] = stat_table$TajimasD_Pachon[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,27] = stat_table$TajimasD_Molino[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,28] = stat_table$TajimasD_Tinaja[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,29] = stat_table$TajimasD_Rascon[stat_table$Stable_Gene_ID == input_vec[i]]
          temp.df[,38] = stat_table$Publication_Name[stat_table$Stable_Gene_ID == input_vec[i]]
        }else{
          temp.df[,6] = rep(NA, num_stat_copies)
          temp.df[,8:29] = rep(NA, num_stat_copies)
          temp.df[,38] = rep(NA, num_stat_copies)
        }
        
        # Check if current gene is found in transcription data. If so, output
        # associated information. If not, output NAs
        if(input_vec[i] %in% condition_control$Gene_stable_ID){
          if(length(condition_control$logFC[
            (condition_control$Gene_stable_ID == input_vec[i]) &
            (grepl("Choy",condition_control$Class))]) != 0){
            temp.df[,30] = rep(condition_control$logFC[
              (condition_control$Gene_stable_ID == input_vec[i]) &
                (grepl("Choy",condition_control$Class))], num_stat_copies)
            temp.df[,34] = rep(condition_control$PValue[
              (condition_control$Gene_stable_ID == input_vec[i]) &
                (grepl("Choy",condition_control$Class))], num_stat_copies)
          }else{
            temp.df[,30] = rep(NA, num_stat_copies)
            temp.df[,34] = rep(NA, num_stat_copies)
          }
          if(length(condition_control$logFC[
            (condition_control$Gene_stable_ID == input_vec[i]) &
            (grepl("Pachon",condition_control$Class))]) != 0){
            temp.df[,31] = rep(condition_control$logFC[
              (condition_control$Gene_stable_ID == input_vec[i]) &
                (grepl("Pachon",condition_control$Class))], num_stat_copies)
            temp.df[,35] = rep(condition_control$PValue[
              (condition_control$Gene_stable_ID == input_vec[i]) &
                (grepl("Pachon",condition_control$Class))], num_stat_copies)
          }else{
            temp.df[,31] = rep(NA, num_stat_copies)
            temp.df[,35] = rep(NA, num_stat_copies)
          }
          if(length(condition_control$logFC[
            (condition_control$Gene_stable_ID == input_vec[i]) &
            (grepl("Molino",condition_control$Class))]) != 0){
            temp.df[,32] = rep(condition_control$logFC[
              (condition_control$Gene_stable_ID == input_vec[i]) &
                (grepl("Molino",condition_control$Class))], num_stat_copies)
            temp.df[,36] = rep(condition_control$PValue[
              (condition_control$Gene_stable_ID == input_vec[i]) &
                (grepl("Molino",condition_control$Class))], num_stat_copies)
          }else{
            temp.df[,32] = rep(NA, num_stat_copies)
            temp.df[,36] = rep(NA, num_stat_copies)
          }
          if(length(condition_control$logFC[
            (condition_control$Gene_stable_ID == input_vec[i]) &
            (grepl("Tinaja",condition_control$Class))]) != 0){
            temp.df[,33] = rep(condition_control$logFC[
              (condition_control$Gene_stable_ID == input_vec[i]) &
                (grepl("Tinaja",condition_control$Class))], num_stat_copies)
            temp.df[,37] = rep(condition_control$PValue[
              (condition_control$Gene_stable_ID == input_vec[i]) &
                (grepl("Tinaja",condition_control$Class))], num_stat_copies)
          }else{
            temp.df[,33] = rep(NA, num_stat_copies)
            temp.df[,37] = rep(NA, num_stat_copies)
          }
        }else{
          temp.df[,30:37] = NA
        }
        geneName = T
        geneID = F
        # If element is present in neither, output an error
      }else if(!(input_vec[i] %in% all.genes_IDs$all_genes) &
               !(input_vec[i] %in% all.genes_IDs$all_IDs)){
        return(paste(c("ERROR: No transcription or statistical data present for
                       the gene", input_vec[i], "."), collapse = " "))
      }
      # Bind temporary datframe to output.df
      output.df <- rbind(output.df, temp.df)
    }
    # If input is NOT a comma-separated string, check whether input is a ID, name,
    # or phrase
  }else if(!grepl(comma, input)){
    if(input %in% c(position_table$Gene_Name, stat_table$Gene_Name,
                    condition_control$Gene_name)){
      # Count how many times this gene occurs in the stat_table
      num_stat_copies <- sum(str_count(stat_table$Gene_Name, input))
      # If number of copies of gene in statistic data is 0, copy the gene at
      # least once so at least one copy is present for transcription data
      if(num_stat_copies == 0){
        num_stat_copies = 1
      }
      # Create temporary dataframe to be appended to final
      temp.df <- data.frame(matrix(nrow = num_stat_copies, ncol = 38))
      # Store input and associated values in appropriate objects and mark input
      # as a gene name
      temp.df[,1] = rep(input, num_stat_copies)
      temp.df[,2] = rep(all.genes_IDs$all_IDs[
        all.genes_IDs$all_genes == input], num_stat_copies)
      
      # If the current gene is present in the position table, output position
      # table info. If not, output all NA
      if(input %in% position_table$Gene_Name){
        temp.df[,3] = rep(position_table$Scaffold[
          position_table$Gene_Name == input], num_stat_copies)
        temp.df[,4] = rep(position_table$Start_Locus[
          position_table$Gene_Name == input], num_stat_copies)
        temp.df[,5] = rep(position_table$End_Locus[
          position_table$Gene_Name == input], num_stat_copies)
      }else{
        temp.df[,3:5] = rep(NA, num_stat_copies)
      }
      
      # Check if gene is present in GO term table. If so, output GO terms. If
      # not, output NA
      if(input %in% GeneToGO$Gene.names){
        temp.df[,7] = rep(GeneToGO$Gene.ontology.IDs[
          GeneToGO$Gene.names == input], num_stat_copies)
      }else{
        temp.df[,7] = rep(NA, num_stat_copies)
      }
      
      # Check if gene is present in stat table. If so, output info. If not,
      # output NAs
      if(input %in% stat_table$Gene_Name){
        temp.df[,6] = stat_table$Gene_Description[stat_table$Gene_Name == input]
        temp.df[,8] = stat_table$Pi_RioChoy[stat_table$Gene_Name == input]
        temp.df[,9] = stat_table$Pi_Pachon[stat_table$Gene_Name == input]
        temp.df[,10] = stat_table$Pi_Molino[stat_table$Gene_Name == input]
        temp.df[,11] = stat_table$Pi_Tinaja[stat_table$Gene_Name == input]
        temp.df[,12] = stat_table$Pi_Rascon[stat_table$Gene_Name == input]
        temp.df[,13] = stat_table$Dxy_RioChoy.Pachon[stat_table$Gene_Name == input]
        temp.df[,14] = stat_table$Dxy_RioChoy.Molino[stat_table$Gene_Name == input]
        temp.df[,15] = stat_table$Dxy_RioChoy.Tinaja[stat_table$Gene_Name == input]
        temp.df[,16] = stat_table$Dxy_Rascon.Pachon[stat_table$Gene_Name == input]
        temp.df[,17] = stat_table$Dxy_Rascon.Tinaja[stat_table$Gene_Name == input]
        temp.df[,18] = stat_table$Dxy_Chica1.Chica2[stat_table$Gene_Name == input]
        
        temp.df[,19] = stat_table$Fst_RioChoy.Pachon[stat_table$Gene_Name == input]
        temp.df[,20] = stat_table$Fst_RioChoy.Molino[stat_table$Gene_Name == input]
        temp.df[,21] = stat_table$Fst_RioChoy.Tinaja[stat_table$Gene_Name == input]
        temp.df[,22] = stat_table$Fst_Pachon.Rascon[stat_table$Gene_Name == input]
        temp.df[,23] = stat_table$Fst_Rascon.Tinaja[stat_table$Gene_Name == input]
        temp.df[,24] = stat_table$Fst_Chica1.Chica2[stat_table$Gene_Name == input]
        
        temp.df[,25] = stat_table$TajimasD_RioChoy[stat_table$Gene_Name == input]
        temp.df[,26] = stat_table$TajimasD_Pachon[stat_table$Gene_Name == input]
        temp.df[,27] = stat_table$TajimasD_Molino[stat_table$Gene_Name == input]
        temp.df[,28] = stat_table$TajimasD_Tinaja[stat_table$Gene_Name == input]
        temp.df[,29] = stat_table$TajimasD_Rascon[stat_table$Gene_Name == input]
        temp.df[,38] = stat_table$Publication_Name[stat_table$Gene_Name == input]
        
      }else{
        temp.df[,6] = rep(NA, num_stat_copies)
        temp.df[,8:29] = rep(NA, num_stat_copies)
        temp.df[,38] = rep(NA, num_stat_copies)
      }
      
      # Check if current gene is found in transcription data. If so, output
      # associated information. If not, output NAs
      if(input %in% condition_control$Gene_name){
        if(length(condition_control$logFC[
          (condition_control$Gene_name == input_vec[i]) &
          (grepl("Choy",condition_control$Class))]) != 0){
          temp.df[,30] = rep(condition_control$logFC[
            (condition_control$Gene_name == input_vec[i]) &
              (grepl("Choy",condition_control$Class))], num_stat_copies)
          temp.df[,34] = rep(condition_control$PValue[
            (condition_control$Gene_name == input_vec[i]) &
              (grepl("Choy",condition_control$Class))], num_stat_copies)
        }else{
          temp.df[,30] = rep(NA, num_stat_copies)
          temp.df[,34] = rep(NA, num_stat_copies)
        }
        if(length(condition_control$logFC[
          (condition_control$Gene_name == input_vec[i]) &
          (grepl("Pachon",condition_control$Class))]) != 0){
          temp.df[,31] = rep(condition_control$logFC[
            (condition_control$Gene_name == input_vec[i]) &
              (grepl("Pachon",condition_control$Class))], num_stat_copies)
          temp.df[,35] = rep(condition_control$PValue[
            (condition_control$Gene_name == input_vec[i]) &
              (grepl("Pachon",condition_control$Class))], num_stat_copies)
        }else{
          temp.df[,31] = rep(NA, num_stat_copies)
          temp.df[,35] = rep(NA, num_stat_copies)
        }
        if(length(condition_control$logFC[
          (condition_control$Gene_name == input_vec[i]) &
          (grepl("Molino",condition_control$Class))]) != 0){
          temp.df[,32] = rep(condition_control$logFC[
            (condition_control$Gene_name == input_vec[i]) &
              (grepl("Molino",condition_control$Class))], num_stat_copies)
          temp.df[,36] = rep(condition_control$PValue[
            (condition_control$Gene_name == input_vec[i]) &
              (grepl("Molino",condition_control$Class))], num_stat_copies)
        }else{
          temp.df[,32] = rep(NA, num_stat_copies)
          temp.df[,36] = rep(NA, num_stat_copies)
        }
        if(length(condition_control$logFC[
          (condition_control$Gene_name == input_vec[i]) &
          (grepl("Tinaja",condition_control$Class))]) != 0){
          temp.df[,33] = rep(condition_control$logFC[
            (condition_control$Gene_name == input_vec[i]) &
              (grepl("Tinaja",condition_control$Class))], num_stat_copies)
          temp.df[,37] = rep(condition_control$PValue[
            (condition_control$Gene_name == input_vec[i]) &
              (grepl("Tinaja",condition_control$Class))], num_stat_copies)
        }else{
          temp.df[,33] = rep(NA, num_stat_copies)
          temp.df[,37] = rep(NA, num_stat_copies)
        }
      }else{
        temp.df[,30:37] = rep(NA, num_stat_copies)
      }
      
      geneName = T
      geneID = F
      
      output.df <- rbind(output.df, temp.df)
      
    }else if(input %in% c(position_table$Gene_ID, stat_table$Stable_Gene_ID,
                          condition_control$Gene_stable_ID)){
      
      # Count how many times this ID occurs in the stat_table
      num_stat_copies <- sum(str_count(stat_table$Stable_Gene_ID, input))
      # If number of copies of gene in statistic data is 0, copy the gene at
      # least once so at least one copy is present for transcription data
      if(num_stat_copies == 0){
        num_stat_copies = 1
      }
      # Create temporary dataframe to later be bound to final df
      temp.df <- data.frame(matrix(nrow = num_stat_copies, ncol = 38))
      # Store input and associated values in appropriate objects and mark input
      # as a gene ID
      temp.df[,1] = rep(all.genes_IDs$all_genes[
        all.genes_IDs$all_IDs == input], num_stat_copies)
      temp.df[,2] = rep(input, num_stat_copies)
      
      # If the current gene is present in the position table, output position
      # table info. If not, output all NA
      if(temp.df[1,2] %in% position_table$Gene_Name){
        temp.df[,3] = rep(position_table$Scaffold[
          position_table$Gene_Name == temp.df[1,2]], num_stat_copies)
        temp.df[,4] = rep(position_table$Start_Locus[
          position_table$Gene_Name == temp.df[1,2]], num_stat_copies)
        temp.df[,5] = rep(position_table$End_Locus[
          position_table$Gene_Name == temp.df[1,2]], num_stat_copies)
      }else{
        temp.df[,3:5] = rep(NA, num_stat_copies)
      }
      
      # Check if gene is present in GO term table. If so, output GO terms. If
      # not, output NA
      if(temp.df[1,2] %in% GeneToGO$Gene.names){
        temp.df[,7] = rep(GeneToGO$Gene.ontology.IDs[
          GeneToGO$Gene.names == temp.df[1,2]], num_stat_copies)
      }else{
        temp.df[,7] = rep(NA, num_stat_copies)
      }
      
      # Check if gene is present in stat table. If so, output info. If not,
      # output NAs
      if(input %in% stat_table$Stable_Gene_ID){
        temp.df[,6] = stat_table$Gene_Description[stat_table$Stable_Gene_ID == input]
        temp.df[,8] = stat_table$Pi_RioChoy[stat_table$Stable_Gene_ID == input]
        temp.df[,9] = stat_table$Pi_Pachon[stat_table$Stable_Gene_ID == input]
        temp.df[,10] = stat_table$Pi_Molino[stat_table$Stable_Gene_ID == input]
        temp.df[,11] = stat_table$Pi_Tinaja[stat_table$Stable_Gene_ID == input]
        temp.df[,12] = stat_table$Pi_Rascon[stat_table$Stable_Gene_ID == input]
        temp.df[,13] = stat_table$Dxy_RioChoy.Pachon[stat_table$Stable_Gene_ID == input]
        temp.df[,14] = stat_table$Dxy_RioChoy.Molino[stat_table$Stable_Gene_ID == input]
        temp.df[,15] = stat_table$Dxy_RioChoy.Tinaja[stat_table$Stable_Gene_ID == input]
        temp.df[,16] = stat_table$Dxy_Rascon.Pachon[stat_table$Stable_Gene_ID == input]
        temp.df[,17] = stat_table$Dxy_Rascon.Tinaja[stat_table$Stable_Gene_ID == input]
        temp.df[,18] = stat_table$Dxy_Chica1.Chica2[stat_table$Stable_Gene_ID == input]
        
        temp.df[,19] = stat_table$Fst_RioChoy.Pachon[stat_table$Stable_Gene_ID == input]
        temp.df[,20] = stat_table$Fst_RioChoy.Molino[stat_table$Stable_Gene_ID == input]
        temp.df[,21] = stat_table$Fst_RioChoy.Tinaja[stat_table$Stable_Gene_ID == input]
        temp.df[,22] = stat_table$Fst_Pachon.Rascon[stat_table$Stable_Gene_ID == input]
        temp.df[,23] = stat_table$Fst_Rascon.Tinaja[stat_table$Stable_Gene_ID == input]
        temp.df[,24] = stat_table$Fst_Chica1.Chica2[stat_table$Stable_Gene_ID == input]
        
        temp.df[,25] = stat_table$TajimasD_RioChoy[stat_table$Stable_Gene_ID == input]
        temp.df[,26] = stat_table$TajimasD_Pachon[stat_table$Stable_Gene_ID == input]
        temp.df[,27] = stat_table$TajimasD_Molino[stat_table$Stable_Gene_ID == input]
        temp.df[,28] = stat_table$TajimasD_Tinaja[stat_table$Stable_Gene_ID == input]
        temp.df[,29] = stat_table$TajimasD_Rascon[stat_table$Stable_Gene_ID == input]
        temp.df[,38] = stat_table$Publication_Name[stat_table$Stable_Gene_ID == input]
      }else{
        temp.df[,6] = rep(NA, num_stat_copies)
        temp.df[,8:29] = rep(NA, num_stat_copies)
        temp.df[,38] = rep(NA, num_stat_copies)
      }
      
      # Check if current gene is found in transcription data. If so, output
      # associated information. If not, output NAs
      if(input %in% condition_control$Gene_stable_ID){
        if(length(condition_control$logFC[
          (condition_control$Gene_stable_ID == input_vec[i]) &
          (grepl("Choy",condition_control$Class))]) != 0){
          temp.df[,30] = rep(condition_control$logFC[
            (condition_control$Gene_stable_ID == input_vec[i]) &
              (grepl("Choy",condition_control$Class))], num_stat_copies)
          temp.df[,34] = rep(condition_control$PValue[
            (condition_control$Gene_stable_ID == input_vec[i]) &
              (grepl("Choy",condition_control$Class))], num_stat_copies)
        }else{
          temp.df[,30] = rep(NA, num_stat_copies)
          temp.df[,34] = rep(NA, num_stat_copies)
        }
        if(length(condition_control$logFC[
          (condition_control$Gene_stable_ID == input_vec[i]) &
          (grepl("Pachon",condition_control$Class))]) != 0){
          temp.df[,31] = rep(condition_control$logFC[
            (condition_control$Gene_stable_ID == input_vec[i]) &
              (grepl("Pachon",condition_control$Class))], num_stat_copies)
          temp.df[,35] = rep(condition_control$PValue[
            (condition_control$Gene_stable_ID == input_vec[i]) &
              (grepl("Pachon",condition_control$Class))], num_stat_copies)
        }else{
          temp.df[,31] = rep(NA, num_stat_copies)
          temp.df[,35] = rep(NA, num_stat_copies)
        }
        if(length(condition_control$logFC[
          (condition_control$Gene_stable_ID == input_vec[i]) &
          (grepl("Molino",condition_control$Class))]) != 0){
          temp.df[,32] = rep(condition_control$logFC[
            (condition_control$Gene_stable_ID == input_vec[i]) &
              (grepl("Molino",condition_control$Class))], num_stat_copies)
          temp.df[,36] = rep(condition_control$PValue[
            (condition_control$Gene_stable_ID == input_vec[i]) &
              (grepl("Molino",condition_control$Class))], num_stat_copies)
        }else{
          temp.df[,32] = rep(NA, num_stat_copies)
          temp.df[,36] = rep(NA, num_stat_copies)
        }
        if(length(condition_control$logFC[
          (condition_control$Gene_stable_ID == input_vec[i]) &
          (grepl("Tinaja",condition_control$Class))]) != 0){
          temp.df[,33] = rep(condition_control$logFC[
            (condition_control$Gene_stable_ID == input_vec[i]) &
              (grepl("Tinaja",condition_control$Class))], num_stat_copies)
          temp.df[,37] = rep(condition_control$PValue[
            (condition_control$Gene_stable_ID == input_vec[i]) &
              (grepl("Tinaja",condition_control$Class))], num_stat_copies)
        }else{
          temp.df[,33] = rep(NA, num_stat_copies)
          temp.df[,37] = rep(NA, num_stat_copies)
        }
      }else{
        temp.df[,30:37] = rep(NA, num_stat_copies)
      }
      geneID = T
      geneName = F
      
      output.df <- rbind(output.df, temp.df)
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
      input_vec <- append(input_vec,
                          c(GeneToGO$Gene.names[
                            grepl(input, GeneToGO$Gene.ontology..biological.process.)
                          ],
                          GeneToGO$Gene.names[
                            grepl(input, GeneToGO$Gene.ontology..cellular.component.)
                          ],
                          GeneToGO$Gene.names[
                            grepl(input, GeneToGO$Gene.ontology..molecular.function.)
                          ]
                          )
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
    
    # For each gene, collect all values associated with the gene
    for(i in 1:length(input_vec)){
      # Find number of copies of this gene in stat table
      num_stat_copies <- sum(str_count(stat_table$Gene_Name, input_vec[i]))
      # If number of copies of gene in statistic data is 0, copy the gene at
      # least once so at least one copy is present for transcription data
      if(num_stat_copies == 0){
        num_stat_copies = 1
      }
      # Create temporary dataframe to be appended to final
      temp.df <- data.frame(matrix(nrow = num_stat_copies, ncol = 38))
      
      temp.df[,1] = rep(input_vec[i], num_stat_copies)
      if(input_vec[i] %in% all.genes_IDs$all_genes){
        temp.df[,2] = rep(all.genes_IDs$all_IDs[
          all.genes_IDs$all_genes == input_vec[i]], num_stat_copies)
      }else{
        temp.df[,2] = rep(NA, num_stat_copies)
      }
      
      # If the current gene is present in the position table, output position
      # table info. If not, output all NA
      if(input_vec[i] %in% position_table$Gene_Name){
        temp.df[,3] = rep(position_table$Scaffold[
          position_table$Gene_Name == input_vec[i]], num_stat_copies)
        temp.df[,4] = rep(position_table$Start_Locus[
          position_table$Gene_Name == input_vec[i]], num_stat_copies)
        temp.df[,5] = rep(position_table$End_Locus[
          position_table$Gene_Name == input_vec[i]], num_stat_copies)
      }else{
        temp.df[,3:5] = rep(NA, num_stat_copies)
      }
      
      # Check if gene is present in GO term table. If so, output GO terms. If
      # not, output NA
      if(input_vec[i] %in% GeneToGO$Gene.names){
        temp.df[,7] = rep(GeneToGO$Gene.ontology.IDs[
          GeneToGO$Gene.names == input_vec[i]], num_stat_copies)
      }else{
        temp.df[,7] = rep(NA, num_stat_copies)
      }
      
      # Check if gene is present in stat table. If so, output info. If not,
      # output NAs
      if(input_vec[i] %in% stat_table$Gene_Name){
        temp.df[,6] = stat_table$Gene_Description[stat_table$Gene_Name == input_vec[i]]
        temp.df[,8] = stat_table$Pi_RioChoy[stat_table$Gene_Name == input_vec[i]]
        temp.df[,9] = stat_table$Pi_Pachon[stat_table$Gene_Name == input_vec[i]]
        temp.df[,10] = stat_table$Pi_Molino[stat_table$Gene_Name == input_vec[i]]
        temp.df[,11] = stat_table$Pi_Tinaja[stat_table$Gene_Name == input_vec[i]]
        temp.df[,12] = stat_table$Pi_Rascon[stat_table$Gene_Name == input_vec[i]]
        temp.df[,13] = stat_table$Dxy_RioChoy.Pachon[stat_table$Gene_Name == input_vec[i]]
        temp.df[,14] = stat_table$Dxy_RioChoy.Molino[stat_table$Gene_Name == input_vec[i]]
        temp.df[,15] = stat_table$Dxy_RioChoy.Tinaja[stat_table$Gene_Name == input_vec[i]]
        temp.df[,16] = stat_table$Dxy_Rascon.Pachon[stat_table$Gene_Name == input_vec[i]]
        temp.df[,17] = stat_table$Dxy_Rascon.Tinaja[stat_table$Gene_Name == input_vec[i]]
        temp.df[,18] = stat_table$Dxy_Chica1.Chica2[stat_table$Gene_Name == input_vec[i]]
        
        temp.df[,19] = stat_table$Fst_RioChoy.Pachon[stat_table$Gene_Name == input_vec[i]]
        temp.df[,20] = stat_table$Fst_RioChoy.Molino[stat_table$Gene_Name == input_vec[i]]
        temp.df[,21] = stat_table$Fst_RioChoy.Tinaja[stat_table$Gene_Name == input_vec[i]]
        temp.df[,22] = stat_table$Fst_Pachon.Rascon[stat_table$Gene_Name == input_vec[i]]
        temp.df[,23] = stat_table$Fst_Rascon.Tinaja[stat_table$Gene_Name == input_vec[i]]
        temp.df[,24] = stat_table$Fst_Chica1.Chica2[stat_table$Gene_Name == input_vec[i]]
        
        temp.df[,25] = stat_table$TajimasD_RioChoy[stat_table$Gene_Name == input_vec[i]]
        temp.df[,26] = stat_table$TajimasD_Pachon[stat_table$Gene_Name == input_vec[i]]
        temp.df[,27] = stat_table$TajimasD_Molino[stat_table$Gene_Name == input_vec[i]]
        temp.df[,28] = stat_table$TajimasD_Tinaja[stat_table$Gene_Name == input_vec[i]]
        temp.df[,29] = stat_table$TajimasD_Rascon[stat_table$Gene_Name == input_vec[i]]
        temp.df[,38] = stat_table$Publication_Name[stat_table$Gene_Name == input_vec[i]]
      }else{
        temp.df[,6] = rep(NA, num_stat_copies)
        temp.df[,8:29] = rep(NA, num_stat_copies)
        temp.df[,38] = rep(NA, num_stat_copies)
      }
      
      # Check if current gene is found in transcription data. If so, output
      # associated information. If not, output NAs
      if(input_vec[i] %in% condition_control$Gene_name){
        if(length(condition_control$logFC[
          (condition_control$Gene_name == input_vec[i]) &
          (grepl("Choy",condition_control$Class))]) != 0){
          temp.df[,30] = rep(condition_control$logFC[
            (condition_control$Gene_name == input_vec[i]) &
              (grepl("Choy",condition_control$Class))], num_stat_copies)
          temp.df[,34] = rep(condition_control$PValue[
            (condition_control$Gene_name == input_vec[i]) &
              (grepl("Choy",condition_control$Class))], num_stat_copies)
        }else{
          temp.df[,30] = rep(NA, num_stat_copies)
          temp.df[,34] = rep(NA, num_stat_copies)
        }
        if(length(condition_control$logFC[
          (condition_control$Gene_name == input_vec[i]) &
          (grepl("Pachon",condition_control$Class))]) != 0){
          temp.df[,31] = rep(condition_control$logFC[
            (condition_control$Gene_name == input_vec[i]) &
              (grepl("Pachon",condition_control$Class))], num_stat_copies)
          temp.df[,35] = rep(condition_control$PValue[
            (condition_control$Gene_name == input_vec[i]) &
              (grepl("Pachon",condition_control$Class))], num_stat_copies)
        }else{
          temp.df[,31] = rep(NA, num_stat_copies)
          temp.df[,35] = rep(NA, num_stat_copies)
        }
        if(length(condition_control$logFC[
          (condition_control$Gene_name == input_vec[i]) &
          (grepl("Molino",condition_control$Class))]) != 0){
          temp.df[,32] = rep(condition_control$logFC[
            (condition_control$Gene_name == input_vec[i]) &
              (grepl("Molino",condition_control$Class))], num_stat_copies)
          temp.df[,36] = rep(condition_control$PValue[
            (condition_control$Gene_name == input_vec[i]) &
              (grepl("Molino",condition_control$Class))], num_stat_copies)
        }else{
          temp.df[,32] = rep(NA, num_stat_copies)
          temp.df[,36] = rep(NA, num_stat_copies)
        }
        if(length(condition_control$logFC[
          (condition_control$Gene_name == input_vec[i]) &
          (grepl("Tinaja",condition_control$Class))]) != 0){
          temp.df[,33] = rep(condition_control$logFC[
            (condition_control$Gene_name == input_vec[i]) &
              (grepl("Tinaja",condition_control$Class))], num_stat_copies)
          temp.df[,37] = rep(condition_control$PValue[
            (condition_control$Gene_name == input_vec[i]) &
              (grepl("Tinaja",condition_control$Class))], num_stat_copies)
        }else{
          temp.df[,33] = rep(NA, num_stat_copies)
          temp.df[,37] = rep(NA, num_stat_copies)
        }
      }else{
        temp.df[,30:37] = rep(NA, num_stat_copies)
      }
      # Add current temporary dataframe to output dataframe
      output.df <- rbind(output.df, temp.df)
    }
  }
  # Remove first row, as it is empty
  output.df <- output.df[-1,]
  
  # Output all values obtained for gene(s) of interest
  names(output.df) <- c(
    "Gene Name",
    "Gene Stable ID",
    "Scaffold",
    "Start Position",
    "Stop Position",
    "Gene Description",
    "GO ID(s)",
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
    "Dxy_Chica1:Chica2",
    "Fst_RioChoy:Pachon",
    "Fst_RioChoy:Molino",
    "Fst_RioChoy:Tinaja",
    "Fst_Pachon:Rascon",
    "Fst_Rascon:Tinaja",
    "Fst_Chica1:Chica2",
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
    "p-value for SD_log(FC)_Tinaja",
    "Publication Name (Population Genetics Data)"
  )
  return(output.df)
}
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
        G_IDs <- append(G_IDs, morph1.rows$Gene_stable_ID[
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

StatDistTable <- function(in_type, UL, stat, thresh, stat_table, pops){
  # Create 2 vectors of the indices corresponding to each population or pop
  # pair's statistic values
  indices <- c()
  pop_strings <- c()
  
  # Create vector into which warnings will be stored for later output
  wrnings <- c("Notes:\n")
  
  # Check if statistic of interest makes comparisons between TWO populations
  if((stat == "Fst") | (stat == "Dxy")){
    # Tell the program that a two-population statistic was entered so the
    # program outputs the specified number of genes for EACH population
    stat_type = "Two Pop"
    # Find all population pairs
    two_pops <- combn(pops,2)
    # If so, iterate through each combination of pops and output the indices
    # and string version of the population pair
    for(pair in 1:ncol(two_pops)){
      # If pops is a matrix, read the strings, find the column corresponding
      # to the stat of interest for the populations in the vector, and set
      # "index" equal to the column housing this statistic
      val <- which(grepl(two_pops[1, pair], names(stat_table))
                   & grepl(two_pops[2, pair], names(stat_table))
                   & grepl(stat, names(stat_table)))
      # If statistic for populations-of-interest is not present, return error
      if(length(val) == 0){
        wrnings <- append(wrnings, (paste(c("Statistic",stat,
                                            "is not present for the populations",two_pops[1, pair],"and",
                                            two_pops[2, pair], "\n"),collapse = " ")))
        # If these populations-of-interest are the only pops which were
        # inputted, return a warning
        if(ncol(two_pops) == 1){
          null.df <- data.frame(matrix(nrow = 1, ncol = 7))
          names(null.df) <- c(
            "Rank",
            "Population(s)",
            stat,
            "Scaffold",
            "Gene Name",
            "GO Term(s)",
            "Publication Name"
          )
          return(list(paste(c("Statistic ",stat,
                              " is not present for the populations ",
                              two_pops[1, pair]," and ",
                              two_pops[2, pair]),collapse = ""), 
                      null.df))
        }
      }else{
        # Create a row to add to the indices dataframe
        temp_str <- paste(c(two_pops[1, pair],"-",two_pops[2, pair]), collapse = "")
        pop_strings <- append(pop_strings,temp_str)
        indices <- append(indices,val)
      }
      # If pops is NOT a matrix, return an error
    }
  }else if((stat == "TajimasD") | (stat == "Pi")){
    # Tell the program that a one-population statistic was entered so the
    # program outputs the specified number of genes across ALL populations
    stat_type = "One Pop"
    # Iterate through each individual population
    for(p in 1:length(pops)){
      # If pops is a string, set indx equal to the column housing the stat of
      # interest for this population
      val <- which(grepl(pops[p], names(stat_table))
                   & grepl(stat, names(stat_table)))
      # If statistic for populations-of-interest is not present, return a
      # warning
      if(length(val) == 0){
        wrnings <- append(wrnings, (paste(c("Statistic",stat,
                                            "is not present for the population",
                                            pops[p],"\n"),collapse = " ")))
      }else{
        pop_strings <- append(pop_strings,pops[p])
        indices <- append(indices,val)
      }
    }
  }else{
    return("ERROR: Invalid statistic name")
  }
  
  # Initialize vectors of gene names, populations, statistic values, and
  # publication names
  genes <- c()
  DF_pops <- c()
  stat_vals <- c()
  pub_names <- c()
  
  
  # Check if statistic value or gene count was entered
  # If value was entered...
  if(in_type == "Statistic Value"){
    # Check whether top or bottom proportion was requested
    # If higher proportion was requested, iterate through each index
    if(UL == "top"){
      for(i in 1:length(indices)){
        # First, remove all NA values for this index
        stat_table <- stat_table[!is.na(stat_table[,indices[i]]),]
        # For each index, collect all genes, scaffolds, populations, values,
        # and publication names whose stat values fall above the entered value
        genes <- append(genes,stat_table$Gene_Name[
          stat_table[,indices[i]] >= thresh])
        DF_pops <- append(DF_pops,rep(pop_strings[i],
                                      length(stat_table$Gene_Name[
                                        stat_table[,indices[i]] >= thresh])))
        stat_vals <- append(stat_vals,stat_table[
          stat_table[,indices[i]] >= thresh,indices[i]])
        pub_names <- append(pub_names,stat_table$Publication_Name[
          stat_table[,indices[i]] >= thresh])
      }
      # If lower tail was requested, iterate through each index
    }else if(UL == "bottom"){
      for(i in 1:length(indices)){
        # First, remove all NA values for this index
        stat_table <- stat_table[!is.na(stat_table[,indices[i]]),]
        # For each index, collect all genes, populations, and statistic values
        # whose values fall in lowest tail
        genes <- append(genes,stat_table$Gene_Name[
          stat_table[,indices[i]] <= thresh])
        DF_pops <- append(DF_pops,rep(pop_strings[i],
                                      length(stat_table$Gene_Name[
                                        stat_table[,indices[i]] <= thresh])))
        stat_vals <- append(stat_vals,stat_table[
          stat_table[,indices[i]] <= thresh,indices[i]])
        pub_names <- append(pub_names,stat_table$Publication_Name[
          stat_table[,indices[i]] <= thresh])
      }
    }
    # If count was entered...
  }else if(in_type == "Gene Count"){
    # Check whether top or bottom proportion was requested
    # If top count was requested, iterate through each index
    if(UL == "top"){
      # Check if statistic is a one or a two population statistic
      # If statistic is a two-population statistic, output the top N genes for ALL
      # possible population pairs
      if(stat_type == "Two Pop"){
        # For each index...
        for(i in 1:length(indices)){
          # First, remove all NA values for this index
          stat_table <- stat_table[!is.na(stat_table[,indices[i]]),]
          # Collect positions of top N genes for the current index
          top_genes <- order(
            stat_table[,indices[i]], decreasing = T)[1:thresh]
          # Collect the genes with the highest values, as well as the associated
          # populations and values
          genes <- append(genes,stat_table$Gene_Name[top_genes])
          DF_pops <- append(DF_pops,rep(pop_strings[i],length(top_genes)))
          stat_vals <- append(stat_vals,stat_table[top_genes,indices[i]])
          pub_names <- append(pub_names,stat_table$Publication_Name[
            top_genes])
        }
        # If statistic is a one-population statistic, output collect the N genes
        # with the HIGHEST stat value, regardless of pop
      }else if(stat_type == "One Pop"){
        # Collect stat values for ALL indices into a 3 vectors: row
        # in one column, population name in another,and stat value in the other
        all_stats <- c()
        all_pops <- c()
        all_genes <- c()
        all_pubs <- c()
        for(i in 1:length(indices)){
          # First, remove all NA values for this index
          stat_table <- stat_table[!is.na(stat_table[,indices[i]]),]
          all_stats <- append(all_stats, stat_table[,indices[i]])
          all_pops <- append(all_pops,rep(pop_strings[i],
                                          length(stat_table[,indices[i]])))
          all_genes <- append(all_genes, stat_table$Gene_Name)
          all_pubs <- append(all_genes, stat_table$Publication_Name)
        }
        # Organize vectors into a dataframe
        temp_df <- data.frame(
          all_stats,
          all_pops,
          all_genes,
          all_pubs
        )
        # Retrieve the parallel indices for the N highest genes
        par_indices <- order(temp_df$all_stats, decreasing = T)[1:thresh]
        # Retrieve the gene names, population names, stat values, and 
        # publication names
        genes <- append(genes,temp_df$all_genes[par_indices])
        DF_pops <- append(DF_pops,temp_df$all_pops[par_indices])
        stat_vals <- append(stat_vals,temp_df$all_stats[par_indices])
        pub_names <- append(pub_names,temp_df$all_pubs[par_indices])
      }
      # If lower proportion was requested, iterate through each index
    }else if(UL == "bottom"){
      # Check if statistic is a one or a two population statistic
      # If statistic is a two-population statistic, output the bottom N genes for ALL
      # possible population pairs
      if(stat_type == "Two Pop"){
        # For each index...
        for(i in 1:length(indices)){
          # First, remove all NA values for this index
          stat_table <- stat_table[!is.na(stat_table[,indices[i]]),]
          # Collect positions of bottom N genes for the current index
          bottom_genes <- order(stat_table[,indices[i]], decreasing = F)[1:thresh]
          # Collect the genes with the highest values, as well as the 
          # associated populations, values, and publications
          genes <- append(genes,stat_table$Gene_Name[bottom_genes])
          DF_pops <- append(DF_pops,rep(pop_strings[i],length(bottom_genes)))
          stat_vals <- append(stat_vals,stat_table[bottom_genes,indices[i]])
          pub_names <- append(pub_names,stat_table$Publication_Name[bottom_genes])
        }
        # If statistic is a one-population statistic, output collect the N genes
        # with the HIGHEST stat value, regardless of pop
      }else if(stat_type == "One Pop"){
        # Collect stat values for ALL indices into a 3 vectors: row
        # in one column, population name in another,and stat value in the other
        all_stats <- c()
        all_pops <- c()
        all_genes <- c()
        all_pubs <- c()
        for(i in 1:length(indices)){
          # First, remove all NA values for this index
          stat_table <- stat_table[!is.na(stat_table[,indices[i]]),]
          all_stats <- append(all_stats, stat_table[,indices[i]])
          all_pops <- append(all_pops,rep(pop_strings[i],
                                          length(stat_table[,indices[i]])))
          all_genes <- append(all_genes, stat_table$Gene_Name)
          all_pubs <- append(all_genes, stat_table$Publication_Name)
        }
        # Organize vectors into a dataframe
        temp_df <- data.frame(
          all_stats,
          all_pops,
          all_genes,
          all_pubs
        )
        # Retreat the parallel indeices for the N highest genes
        par_indices <- order(temp_df$all_stats, decreasing = F)[1:thresh]
        # Retrieve the gene names, population names
        genes <- append(genes,temp_df$all_genes[par_indices])
        DF_pops <- append(DF_pops,temp_df$all_pops[par_indices])
        stat_vals <- append(stat_vals,temp_df$all_stats[par_indices])
        pub_names <- append(pub_names,temp_df$all_pubs[par_indices])
      }
    }
  }
  # If no population pairs were found, output an error
  if(length(indices) == 0){
    null.df <- data.frame(matrix(nrow = 1, ncol = 7))
    names(null.df) <- c(
      "Rank",
      "Population(s)",
      stat,
      "Scaffold",
      "Gene Name",
      "GO Term(s)",
      "Publication Name"
    )
    return(list(paste(c("Statistic",stat,"not present for the selected population(s)"),
                      collapse = " "), null.df))
  }
  
  # Initialize a vector of GO terms and a vector of scaffolds
  scaffs <- character(length = length(genes))
  DF_GOs <- character(length = length(genes))
  
  # For each gene in the gene vector, output the associated scaffolds to the
  # vector of scaffolds and output the associated GO terms to the vector of GO
  # terms
  for(g in 1:length(genes)){
    # If gene is present in positions table and GO table, obtain scaffold and GO
    # terms
    if((genes[g] %in% position_table$Gene_Name) & (genes[g] %in% GeneToGO$Gene.names)){
      scaffs[g] = position_table$Scaffold[position_table$Gene_Name == genes[g]]
      DF_GOs[g] = GeneToGO$Gene.ontology.IDs[GeneToGO$Gene.names == genes[g]]
      # If gene is present in position table but NOT GO table, output NA for GO
      # term but output real scaffold
    }else if((genes[g] %in% position_table$Gene_Name) & !(genes[g] %in% GeneToGO$Gene.names)){
      scaffs[g] = position_table$Scaffold[position_table$Gene_Name == genes[g]]
      DF_GOs[g] = "Not applicable"
      # If gene is present in GO table but NOT position table, output NA for scaffold
      # term but output real GO
    }else if(!(genes[g] %in% position_table$Gene_Name) & (genes[g] %in% GeneToGO$Gene.names)){
      scaffs[g] = "Not applicable"
      DF_GOs[g] = GeneToGO$Gene.ontology.IDs[GeneToGO$Gene.names == genes[g]]
      # If gene is present in neither GO nor position tables, output NA for scaff
      # and GO
    }else if(!(genes[g] %in% position_table$Gene_Name) & !(genes[g] %in% GeneToGO$Gene.names)){
      scaffs[g] = "Not applicable"
      DF_GOs[g] = "Not applicable"
    }
  }
  # Create a dataframe of scaffolds, gene names, GO terms, statistic types, and
  # stat values
  prelim_df <- data.frame(
    DF_pops,
    stat_vals,
    scaffs,
    genes,
    DF_GOs,
    pub_names
  )
  # Sort into a final dataframe
  final_df <- data.frame()
  
  
  # Sort df based on statistic values
  # If "top" was checked, sort dataframe in DESCENDING order from top to bottom
  # but retain population groups
  if(UL == "top"){
    # If statistic is a two-population statistic, sort WITHIN populations PRIOR
    # to sorting by statistic
    if(stat_type == "Two Pop"){
      # Create vector in which to store number of genes falling above threshhold
      # for each population pair
      ranks <- c()
      for(i in 1:length(pop_strings)){
        rows <- prelim_df[prelim_df$DF_pops == pop_strings[i],]
        new_rows <- rows[order(rows$stat_vals, decreasing = T),]
        # Collect number of genes found above threshold for THIS population
        ranks <- append(ranks, 1:nrow(new_rows))
        final_df <- rbind(final_df, new_rows)
      }
      # Add column which ranks within each population pair
      final_df <- add_column(final_df, ranks,.before = T)
      # If statistic is a one-population statistic, sort ONLY by statistic value
    }else if(stat_type == "One Pop"){
      final_df <- prelim_df[order(stat_vals, decreasing = T),]
      final_df <- add_column(final_df, 1:nrow(prelim_df),
                             .before = T)
    }
    # If "bottom" was checked, sort dataframe in ASCENDING order from top to bottom
    # but retain population groups
  }else if(UL == "bottom"){
    if(stat_type == "Two Pop"){
      # Create vector in which to store number of genes falling above threshhold
      # for each population pair
      ranks <- c()
      for(i in 1:length(pop_strings)){
        rows <- prelim_df[prelim_df$DF_pops == pop_strings[i],]
        new_rows <- rows[order(rows$stat_vals, decreasing = F),]
        # Collect number of genes found above threshhold for THIS population
        ranks <- append(ranks, 1:nrow(new_rows))
        final_df <- rbind(final_df, new_rows)
      }
      # Add column which ranks within each population pair
      final_df <- add_column(final_df, ranks,.before = T)
    }else if(stat_type == "One Pop"){
      final_df <- prelim_df[order(stat_vals, decreasing = F),]
      # Add rank column
      final_df <- add_column(final_df, 1:nrow(prelim_df),
                             .before = T)
    }
  }
  
  # Rename df stat type and stat value columns with specific statistic's name
  names(final_df) <- c(
    "Rank",
    "Population(s)",
    stat,
    "Scaffold",
    "Gene Name",
    "GO Term(s)",
    "Publication Name"
  )
  
  # Output df and warnings
  if(is.na(final_df[1,2])){
    final_df <- final_df[-1,]
  }
  return(list(wrnings, final_df))
}
StatDistPlot <- function(stat, UL, thresh, stat_table, pops){
  # EDIT: Check if threshold is within appropriate range for statistic-of-interest. If
  # not, return an error.
  
  # Initialize vector into which notes will be stored
  wrnings <- c("Notes: ")
  # Create an error plot to be returned if an error occurs
  error_plot <- plot(x = c(.95, .955,  .975,  1, 1.025, 1.05,   1.07,  1.075, 0.99, 1.035), 
                     y = c( 1, 1.5, 1.75, 2, 2,    1.75,  1.5,   1,    3.5,    3.5), pch = 16, axes = F,
                     ylab = "",
                     xlab = "Whoops! Something went wrong. See errors")
  # Initialize vector of indices
  indices <- c()
  
  # Check if statistic of interest makes comparisons between TWO populations
  if((stat == "Fst") | (stat == "Dxy")){
    # Tell the program that a two-population statistic was entered so the
    # program outputs the specified number of genes for EACH population
    stat_type = "Two Pop"
    
    if(length(pops) == 1){
      return(list("ERROR: Only one population supplied for a two-population statistic",
                  error_plot))
    }
    # Find all population pairs
    two_pops <- combn(pops,2)
    # If so, iterate through each combination of pops and output the indices
    # and string version of the population pair
    for(pair in 1:ncol(two_pops)){
      # Read the strings, find the column corresponding
      # to the stat of interest for the populations in the vector, and set
      # "indx" equal to the column housing this statistic
      indx <- which(grepl(two_pops[1, pair], names(stat_table))
                    & grepl(two_pops[2, pair], names(stat_table))
                    & grepl(stat, names(stat_table)))
      
      # If statistic for populations-of-interest is not present, return error
      if(length(indx) == 0){
        wrnings <- append(wrnings, paste(c("Statistic",stat,
                                           "is not present for the populations",two_pops[1, pair],"and",
                                           two_pops[2, pair]),collapse = " "))
        # If these populations-of-interest are the only pops which were
        # inputted, return a warning
        if(ncol(two_pops) == 1){
          return(list(paste(c("Statistic ",stat,
                              " is not present for the populations ",
                              two_pops[1, pair]," and ",
                              two_pops[2, pair]),collapse = ""), 
                      error_plot))
        }
      }else{
        indices <- append(indices, indx)
        next
      }
    }
    # Check if statistic of interest makes comparisons between ONE population
  }else if((stat == "TajimasD") | (stat == "Pi")){
    # Tell the program that a one-population statistic was entered so the
    # program outputs the specified number of genes across ALL populations
    stat_type = "One Pop"
    # Iterate through each individual population
    for(p in 1:length(pops)){
      # If pops is a string, set indx equal to the column housing the stat of
      # interest for this population
      indx <- append(indices, which(grepl(pops[p], names(stat_table))
                                    & grepl(stat, names(stat_table))))
      # If statistic for populations-of-interest is not present, return note
      if(length(indx) == 0){
        wrnings <- append(wrnings, paste(c("Statistic",stat,
                                           "is not present for the population",pops[p]),
                                         collapse = " "))
        next
      }else{
        indices <- append(indices, indx)
        next
      }
    }
  }
  # If no population pairs were found, output an error
  if(length(indices) == 0){
    return(list(paste(c("Statistic",stat,"not present for the selected population(s)"),
                      collapse = " "), error_plot))
  }
  
  # Create dataframe with values of stat-of-interest for every single pop-of
  # -interest
  temp = c()
  for(i in 1:length(indices)){
    temp <- append(temp, stat_table[,indices[i]])
  }
  # Remove missing vlaues so plotting function will work
  temp <- temp[!is.na(temp)]
  
  # Coerce statistics into a one column dataframe so you can use 
  # ShadedDensity()
  filt_table <- as.data.frame(temp)
  
  # Rename df statistics column so you can reference the statistic values i
  # in the plotting function
  names(filt_table)[1] <- "Statistic_Values"
  
  # Pick a direction to shade in based on whether "top" or "bottom" was 
  # specified
  if(UL == "top"){
    t = "right"
  }else if(UL == "bottom"){
    t = "left"
  }
  
  # Create a string to be used in the title of the density plot
  pop_string = ""
  for(p in 1:length(pops)){
    if(p == 1){
      pop_string <- pops[p]
    }else if((p == 2) & (p == length(pops))){
      pop_string <- paste(c(pop_string, " and ", pops[p]), collapse = "")
    }else if((p != length(pops)) & (p > 2)){
      pop_string <- paste(c(pop_string, ", ", pops[p]), collapse = "")
    }else if((p == length(pops)) & (p != 2)){
      pop_string <- paste(c(pop_string, pops[p]), collapse = ", and ")
    }
  }
  
  # Create density plot. Title must be different depending on whether you 
  # have one or two populations
  if(stat_type == "Two Pop"){
    plot_title <- paste(c("Pairwise ", stat, " values for ",pop_string), 
                        collapse = "")
  }else{
    plot_title <- paste(c(stat, " values for ", pop_string),collapse = "")
  }
  
  output_plot <- ShadedDensity(frame = filt_table,
                               xvar = "Statistic_Values",
                               threshold = thresh,
                               title = plot_title,
                               tail = t,
                               shading = "maroon")
  return(list(wrnings, output_plot))
}

StatByChrTable <- function(GOTerm, GeneToGo, GoIDToNames, UpperLower, 
                           stat_vec, position_table, stat_table, pops){
  # Initialize vector for all GO terms of interest
  GOs <- c()
  # Initialize vector in which to store warnings
  wrnings <- c("Notes: ")
  
  # Check if user inputted a word/phrase or a GO ID. If the user inputted a 
  # word/phrase, find the name(s) which contains that word/phrase and set the GO
  # vector equal to the corresponding GO IDs.
  if(!(GOTerm %in% GoIDToNames$GO.ID)){
    if(length(GoIDToNames$GO.ID[which(grepl(GOTerm, GoIDToNames$GO.Term, 
                                            ignore.case = T))]) != 0){
      GOs <- GoIDToNames$GO.ID[which(grepl(GOTerm, GoIDToNames$GO.Term, 
                                           ignore.case = T))]
      # If the input is not found in the GO terms but does have a comma,
      # output both as GO IDs and check whether they are valid GO IDs
    }else if((length(GoIDToNames$GO.ID[which(grepl(GOTerm, GoIDToNames$GO.Term, 
                                                   ignore.case = T))]) == 0)
             & (grepl(", ", GOTerm, ignore.case = T))){
      GOs <- str_split(GOTerm, ", ")[[1]]
      # If the input is not found in the GO terms and does not have a comma,
      # return an error
    }else if((length(GoIDToNames$GO.ID[which(grepl(GOTerm, GoIDToNames$GO.Term, 
                                                   ignore.case = T))]) == 0)
             & !(grepl(", ", GOTerm, ignore.case = T))){
      null.df <- data.frame(matrix(nrow = 1, ncol = 9))
      names(null.df) <- c("Gene",
                          "Scaffold",
                          "Start_Position",
                          "End_Position",
                          "GO_IDs",
                          "Statistic_Type",
                          "Population",
                          "Statistic_Value",
                          "Publication_Name")
      return(list("ERROR: Input is neither a valid GO ID nor a valid phrase", 
                  null.df))
    }
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
    found = F
    # If the GO term appears in the data frame of names AND the corresponding gene
    # occurs in the statistics vector, add the GO term and gene name
    for(i in 1:length(GeneToGO$Gene.ontology.IDs)){
      if(grepl(GOs[g], GeneToGO$Gene.ontology.IDs[i])){
        found = T
      }else{
        next
      }
    }
    if(found){
      gene_vec <- append(gene_vec, GeneToGO$Gene.names[grepl(GOs[g], 
                                                             GeneToGO$Gene.ontology.IDs)])
      found_GOs <- append(found_GOs, rep(GOs[g], 
                                         length(GeneToGO$Gene.names[grepl(GOs[g], 
                                                                          GeneToGO$Gene.ontology.IDs)])))
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
  if(("Fst" %in% stat_vec) | ("Dxy" %in% stat_vec)){
    all_pops <- combn(pops,2) 
    two_pop <- T
  }else{
    two_pop <- F
  }
  
  # Create a temporary dataframe for all possible combinations of statistic
  # types, columns, and the populations to which they correspond 
  stat_pop_combos <- data.frame(matrix(ncol = 3))
  colnames(stat_pop_combos) <- c("Stat_Type", "Pop", "Col")
  
  for(s in 1:length(stat_vec)){
    # Check if statistic of interest makes comparisons between TWO populations
    if((stat_vec[s] == "Fst") | (stat_vec[s] == "Dxy")){
      # Check if all_pops was created
      if(two_pop){
        # If so, iterate through each combination and output stat value for that
        # combination
        for(pair in 1:ncol(all_pops)){
          
          # If pops is a matrix, read the strings, find the column corresponding 
          # to the stat of interest for the populations in the vector, and set 
          # "index" equal to the column housing this statistic
          val <- which(grepl(all_pops[1, pair], names(stat_table))
                       & grepl(all_pops[2, pair], names(stat_table))
                       & grepl(stat_vec[s], names(stat_table)))
          # If statistic for populations-of-interest is not present, add a
          # warning
          if(length(val) == 0){
            wrnings <- append(wrnings, paste(c("Statistic ",stat_vec[s],
                                               " is not present for the populations ",
                                               all_pops[1, pair]," and ",
                                               all_pops[2, pair], " | "),collapse = ""))
            # If these populations-of-interest are the only pops which were
            # inputted, return a warning
            if(ncol(all_pops) == 1){
              null.df <- data.frame(matrix(nrow = 1, ncol = 9))
              names(null.df) <- c("Gene",
                                  "Scaffold",
                                  "Start_Position",
                                  "End_Position",
                                  "GO_IDs",
                                  "Statistic_Type",
                                  "Population",
                                  "Statistic_Value",
                                  "Publication_Name")
              return(list(paste(c("Statistic ",stat_vec[s],
                                  " is not present for the populations ",
                                  all_pops[1, pair]," and ",
                                  all_pops[2, pair]),collapse = "")), 
                     null.df)
            }
          }else{
            # Create a row to add to the indices data frame
            temp_str <- paste(c(all_pops[1, pair],"-",all_pops[2, pair]), 
                              collapse = "")
            temp_vec <- c(stat_vec[s],temp_str,val)
            stat_pop_combos <- rbind(stat_pop_combos, temp_vec)
          }
          # If pops is NOT a matrix, return an error
        }
      }else{
        null.df <- data.frame(matrix(nrow = 1, ncol = 9))
        names(null.df) <- c("Gene",
                            "Scaffold",
                            "Start_Position",
                            "End_Position",
                            "GO_IDs",
                            "Statistic_Type",
                            "Population",
                            "Statistic_Value",
                            "Publication_Name")
        return(list(paste(c("ERROR: Only one population, ", pops, 
                            ", supplied for the two-population statistic ", 
                            stat_vec[s], "."),
                          collapse = ""), null.df))
        
        # Check if statistic of interest makes comparisons between ONE population
      }
    }else if((stat_vec[s] == "TajimasD") | (stat_vec[s] == "Pi")){
      # Iterate through each individual population
      for(p in 1:length(pops)){
        # If pops is a string, set indx equal to the column housing the stat of
        # interest for this population
        val <- which(grepl(pops[p], names(stat_table))
                     & grepl(stat_vec[s], names(stat_table)))
        # If statistic for populations-of-interest is not present, return a 
        # warning 
        if(length(val) == 0){
          wrnings <- append(wrnings, paste(c("Statistic ",stat_vec[s],
                                             " is not present for the population ",
                                             pops[p], " | "), collapse = ""))
        }else{
          # Create a row to add to the indices dataframe
          temp_vec <- c(stat_vec[s],pops[p],val)
          stat_pop_combos <- rbind(stat_pop_combos, temp_vec)
        }
      }
    }
  }
  # If NONE of the populations-of-interest had values for the statistics-of-
  # interest, output an error
  if(nrow(stat_pop_combos) == 1){
    null.df <- data.frame(matrix(nrow = 1, ncol = 9))
    names(null.df) <- c("Gene",
                        "Scaffold",
                        "Start_Position",
                        "End_Position",
                        "GO_IDs",
                        "Statistic_Type",
                        "Population",
                        "Statistic_Value",
                        "Publication_Name")
    return(list("ERROR: None of the input statistics are present for any of the input populations", 
                null.df)) 
  }
  # If dataframe was valid, remove first row (as this is all NAs from)
  # initialization
  stat_pop_combos <- stat_pop_combos[-1,]
  
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
  # 9. Publication from which the statistic value was derived
  
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
  Pub_Name <- c()
  
  for(i in 1:nrow(stat_pop_combos)){
    for(g in 1:length(geneGOs$Gene)){
      # Check if current gene is in stat AND position table
      # If so...
      if((geneGOs$Gene[g] %in% stat_table$Gene_Name) & 
         (geneGOs$Gene[g] %in% position_table$Gene_Name)){
        # Find the number of copies of the current gene in the statistic 
        # table
        copies <- length(stat_table[stat_table$Gene_Name == geneGOs$Gene[g], 
                                    as.numeric(stat_pop_combos$Col[i])])
        # Output the current gene as many times as there are copies
        Gene <- append(Gene, rep(geneGOs$Gene[g], copies))
        # Output scaffold of current gene as many times as there are copies
        Scaffold <- append(Scaffold, 
                           rep(position_table$Scaffold[position_table$Gene_Name == geneGOs$Gene[g]], 
                               copies))
        # Output starting position of the current gene as many times as 
        # there are copies of the gene
        Start_Position <- append(Start_Position, 
                                 rep(position_table$Start_Locus[position_table$Gene_Name == geneGOs$Gene[g]],
                                     copies))
        # Output the ending position of the current gene as many times as 
        # there are copies of the gene
        End_Position <- append(End_Position, 
                               rep(position_table$End_Locus[position_table$Gene_Name == geneGOs$Gene[g]],
                                   copies))
        # Output ALL GO terms associated with the current gene as many times
        # as there are copies of the gene
        GO_IDs <- append(GO_IDs,
                         rep(
                           paste(geneGOs$GO_ID[geneGOs$Gene == geneGOs$Gene[g]], collapse = "; "),
                           copies))
        # Output the current statistic type as many times as there are copies
        # of the gene
        Statistic_Type <- append(Statistic_Type, 
                                 rep(stat_pop_combos$Stat_Type[i], 
                                     copies))
        # Output the population(s) which correspond to the current column,
        # and output as many times as there are copies of the gene.
        Population <- append(Population, rep(stat_pop_combos$Pop[i], 
                                             copies))
        # Output the statistic value. Will automatically output value for
        # each copy of the gene.
        Statistic_Value <- append(Statistic_Value, 
                                  stat_table[stat_table$Gene_Name == geneGOs$Gene[g], 
                                             as.numeric(stat_pop_combos$Col[i])])
        # Output the publication from which the statistic was obtained
        Pub_Name <- append(Pub_Name, stat_table$Publication_Name[which(
          stat_table$Gene_Name == geneGOs$Gene[g]
        )])
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
                          Statistic_Value,
                          Pub_Name
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
                                            & (Full_Table$Statistic_Type == unique_stats[s])],
        Pub_Name = Full_Table$Pub_Name[(Full_Table$Scaffold == unique_scaffs[c])
                                       & (Full_Table$Statistic_Type == unique_stats[s])]
      )
      
      # Plot the statistic values for all genes which occur on the scaffold
      # of interest, colored by population or population pair
      plist[[p]] <- ggplot(plot_dat, aes(x=Locus, y=Value, color=Populations)) + 
        ylab(stat_vec[s]) +
        xlab("Locus") +
        ggtitle(paste(c("Scaffold:", unique_scaffs[c]), collapse = " ")) +
        geom_point(aes(shape = Pub_Name)) +
        theme_bw()
      #Store the plot in a vector
    }
  }
  # Add names to plot list
  names(plist) <- names_vec
  return(plist)
}
# In this function, the user inputs GO ID(s) or a phrase and the function 
# outputs the class, lower-level GO IDs, and GO term associated with each
# relevant GO ID
GOInfo <- function(GO_input, GO_classes, GOIDToNames, UpperLower, all.GO_IDs){
  # Initialize a vector in which to store GO IDs
  GO_ID_vec <- c()
  
  # Initialize an empty dataframe to output in the event of an error
  error.df <- data.frame(matrix(nrow = 1, ncol = 4))
  names(error.df) <- c("GO ID",
                       "GO Term",
                       "Namespace",
                       "All Nested GO IDs")
  
  # Initialize a vector in which to output warnings/errors
  wrnings <- c("Notes: ")
  
  # Find whether input is comma-separated list of GO IDs, single GO ID, or
  # phrase and fill GO_ID_vec with appropriate GO IDs based on answer
  if((grepl(", ", GO_input))){
    # If input is a comma-separated list of GO IDs, parse string and add
    # each GO ID to the vector of GO IDs
    temp_GO_ID_vec <- str_split(string = GO_input, pattern = ", ")[[1]]
    # Check to ensure that at least one of the inputted IDs is a real GO ID
    for(i in 1:length(temp_GO_ID_vec)){
      if(temp_GO_ID_vec[i] %in% all.GO_IDs){
        GO_ID_vec <- append(GO_ID_vec, temp_GO_ID_vec[i])
        next
      }else{
        wrnings <- append(wrnings, paste(c("Input ",temp_GO_ID_vec[i], 
                                           " is not a GO ID."), 
                                         collapse = ""))
        next
      }
    }
    # If NONE of the inputs are real GO IDs, output an error 
    if(is.null(GO_ID_vec)){
      return(list("ERROR: None of the inputs in the comma-separated list are GO IDs",
                  error.df))
    }
    
  }else if(!(grepl(", ", GO_input)) & (GO_input %in% all.GO_IDs)){
    # If input is a single GO ID, add it to the vector of GO IDs
    GO_ID_vec <- GO_input
  }else if(!(grepl(", ", GO_input)) & !(GO_input %in% all.GO_IDs)){
    # If input is a phrase, find all GO terms associated with that phrase
    if(sum(grepl(GO_input, GoIDToNames$GO.Term, ignore.case = T)) != 0){
      GO_ID_vec <- GoIDToNames$GO.ID[grepl(GO_input, GoIDToNames$GO.Term, 
                                           ignore.case = T)]
      # If NO GO terms are associated with the phrase, output an error
    }else{
      return(list(paste(c("ERROR: Input ", GO_input,
                          " is neither a GO ID nor a phrase associated with any recorded GO IDs"),
                        collapse = ""), error.df))
    }
    
  }
  
  # Create a dataframe with the GO terms, namespaces, and lower-level GO IDs
  # associated with all input GO IDs
  GO_df <- data.frame(matrix(nrow = length(GO_ID_vec), ncol = 4))
  names(GO_df) <- c("GO ID",
                    "GO Term",
                    "Namespace",
                    "All Nested GO IDs")
  for(GO in 1:length(GO_ID_vec)){
    # Output GO ID
    GO_df$`GO ID`[GO] <- GO_ID_vec[GO]
    # Output GO term, if present
    if(GO_ID_vec[GO] %in% GoIDToNames$GO.ID){
      GO_df$`GO Term`[GO] <- GoIDToNames$GO.Term[GoIDToNames$GO.ID == GO_ID_vec[GO]]
    }else{
      GO_df$`GO Term`[GO] <- "Not applicable"
    }
    # Output class of GO ID
    if(GO_ID_vec[GO] %in% GO_classes$GO_ID){
      GO_df$Namespace[GO] <- GO_classes$Namespace[GO_classes$GO_ID == GO_ID_vec[GO]]
    }else{
      GO_df$Namespace[GO] <- "Not applicable"
    }
    # Output all lower-level GO IDs associated with current GO ID
    # First, add the current GO ID to a vector
    lower_GOs <- GO_ID_vec[GO]
    for(g in 1:length(lower_GOs)){
      # Add all "Lower" GO IDs which occur on a row where the current vector entry
      # is an "Upper" to the vector of GO IDs, then move to the next GO ID
      if(lower_GOs[g] %in% UpperLower$Upper){
        lower_GOs <- append(lower_GOs, 
                            UpperLower$Lower[UpperLower$Upper == lower_GOs[g]])
        # If the current GO ID does NOT occur anywhere in the "Upper" column,
        # skip it
      }else{
        next
      }
    }
    # Remove the first GO ID from the vector of lower GO IDs, as this is the
    # original GO ID
    lower_GOs <- lower_GOs[-1]
    if(length(lower_GOs) != 0){
      GO_df$`All Nested GO IDs`[GO] <- paste(lower_GOs, collapse = "; ")
    }else{
      GO_df$`All Nested GO IDs`[GO] <- "No lower-level GO IDs"
    }
    
  }
  
  # Return data frame and warnings, if applicable
  return(list(wrnings, GO_df))
}

# This function finds the minimum and maximum values for a population-specific
# statistic in a statistical table
MinMax <- function(mm_pops, mm_stat, stat_table){
  # Create matrix to store all columns which house the morphs-of-interest
  cols.w.morphs <- data.frame(matrix(nrow = length(mm_pops), 
                                     ncol = ncol(stat_table)))
  colnames(cols.w.morphs) <- colnames(stat_table)
  
  # Find all columns of stat_table which house at least one morph-of-interest
  for(m in 1:length(mm_pops)){
    cols.w.morphs[m,] <- grepl(mm_pops[m], names(stat_table), 
                               ignore.case = T)
  }
  cols.of.interest <- c()
  for(i in 1:ncol(cols.w.morphs)){
    if(sum(cols.w.morphs[,i]) > 0){
      cols.of.interest <- append(cols.of.interest, names(cols.w.morphs)[i])
    }
  }
  # Of the columns of stat_table which house at least one morph-of-interest,
  # exclude the columns which have morphs which are NOT morphs of interest
  bad.morphs <- c("Pachon", "Rascon", "RioChoy", "Molino", "Tinaja")
  for(i in 1:length(mm_pops)){
    bad.morphs <- bad.morphs[bad.morphs != mm_pops[i]]
  }
  for(b in 1:length(bad.morphs)){
    if(sum(grepl(bad.morphs[b], cols.of.interest,ignore.case = T)) != 0){
      cols.of.interest <- cols.of.interest[!(grepl(bad.morphs[b], 
                                                   cols.of.interest, 
                                                   ignore.case = T))]
    }
  }
  # Find columns-of-interest which house stat-of-interest
  cols.of.interest <- cols.of.interest[grepl(mm_stat, 
                                             cols.of.interest, 
                                             ignore.case = T)]
  # Combine columns-of-interest into a single vector of numbers
  num_pool <- c()
  for(i in 1:length(cols.of.interest)){
    num_pool <- append(num_pool, 
                       stat_table[,names(stat_table) == cols.of.interest[i]])
  }
  
  # Output the minimum and maximum values in that vector
  minimum <- num_pool[which.min(num_pool)]
  maximum <- num_pool[which.max(num_pool)]
  return(list(minimum,maximum))
}


