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
library(maps)
library(ggrepel)

################################## Load Data ###################################
# Dataframe of latitudes and longitudes for all morphs
# Citation:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3282648/pdf/1471-2148-12-9.pdf
# Rascon and Rio Choy data obtained from google maps
Latit_Longit <- read.csv("data/PopulationLocations.csv")

# Get a map of the world
world_map_1 <- map_data("world")
# get map of USA, Mexico, Belize, and Guatemala
world_map <- subset(world_map_1, region %in% c("USA", "Mexico", "Belize","Guatemala"))


# Read in gene position data
position_table <- read.csv("data/AmexPositionTable.csv", fill = TRUE)
# Remove duplicates so multiple rows with the same values are not outputted
position_table <- position_table[!duplicated(position_table$Gene_ID),]

condition_control <- read.csv("data/Transcription_Morph-Control_Remapped2022.csv")
condition_control$Publication[condition_control$Publication == "McGaugh_et_al_2020"] <- rep("5", sum(condition_control$Publication == "McGaugh_et_al_2020"))

morph1.morph2 <- read.csv("data/Transcription_Morph-Morph.csv")

morph1.morph2$Publication[morph1.morph2$Publication == "Mack_et_al_2020"] <- rep("4", sum(morph1.morph2$Publication == "Mack_et_al_2020"))

morph1.morph2$Publication[morph1.morph2$Publication == "McGaugh_et_al_2020"] <- rep("5", sum(morph1.morph2$Publication == "McGaugh_et_al_2020"))

GeneToGO <- read.csv("data/AMexGOTerms.csv", fill = T)

MasterGO <- read.csv("data/MasterGO.csv", fill = T)

UpperLower <- read.table("data/GOTermAssociations.txt", fill = T, sep = "\t", header = T)

# Read in population genetics table and change publication names
stat_table <- read.csv("data/PopulationGenetics.csv")
stat_table$Publication[stat_table$Publication == "Herman_et_al_2018"] <- "1"
stat_table$Publication[stat_table$Publication == "Moran_et_al_2022"] <- "2"

################################## Functions ###################################
# Gene Search Page: Input a single or comma-separated list of genes, gene IDs,
# GO terms, or a phrase associated with a gene-of-interest and output all 
# available data into distinct tables based on data type
GeneSearch <- function(input, posBool, transcBool, popgenBool, GOBool, 
                       position_table, morph1.morph2, condition_control,
                       stat_table, GeneToGO){
  # Initialize warnings
  warn <- c("Notes: ")
  
  # Output warnings if data was not requested so the user does not get confused
  # by the lack of output
  if(posBool == F){
    warn <- append(warn, "Position table is blank because position information was not requested.")
  }
  if(transcBool == F){
    warn <- append(warn, "Transcription table is blank because transcription information was not requested.")
  }
  if(popgenBool == F){
    warn <- append(warn, "Population genetics table is blank because Population genetics information was not requested.")
  }
  if(GOBool == F){
    warn <- append(warn, "GO table is blank because GO information was not requested.")
  }
  # Define a comma
  comma <- ", "
  
  
  # Define a vector in which to store all inputs
  input_vec <- c()
  
  # Check if comma occurs in input string
  if(grepl(comma, input)){
    # If so, split the input string and add each individual input to the input
    # vector
    input_vec <- str_split(input, pattern = comma)[[1]]
  }else{
    # If not, add the whole string to the input vector
    input_vec <- input
  }
  
  # Define a vector in which to store gene IDs
  geneIDs <- c()
  
  # Iterate through each entry in vector of inputs to find the IDs corresponding
  # to these inputs
  for(i in 1:length(input_vec)){
    # Check if current entry is a gene ID
    if(grepl("ENSAMXG", input_vec[i])){
      # If current entry IS a gene ID, simply append to vector of IDs
      geneIDs <- append(geneIDs, input_vec[i])
    }else{
      # If current entry is NOT a gene ID, convert to lowercase and grep all
      # lower-case gene name and gene info columns for this string.
      searchTerm <- tolower(input_vec[i])
      tempIDs <- c()
      
      # grep every gene name in position table (use grep instead of %in% for 
      # partial gene name matching)
      for(g in 1:nrow(position_table)){
        if(grepl(searchTerm, tolower(position_table$Gene_Name[g]))){
          tempIDs <- append(tempIDs, position_table$Gene_ID[g])
        }
      }
      
      # statistic table
      for(s in 1:nrow(stat_table)){
        # If current entry is part of the gene name, description, OR GO terms,
        # add the ID to the vector of IDs
        if((grepl(searchTerm, tolower(stat_table$Gene_Name[s]))) |
           (grepl(searchTerm, tolower(stat_table$Gene_Description[s]))) |
           (grepl(searchTerm, tolower(stat_table$GO_Terms[s])))){
          tempIDs <- append(tempIDs, stat_table$Stable_Gene_ID[s])
        }
      }
      
      # transcription tables
      for(t in 1:nrow(morph1.morph2)){
        # If current entry is part of the gene name or study-specific details,
        # add the ID to the vector of IDs
        if((grepl(searchTerm, tolower(morph1.morph2$Gene_name[t]))) |
           (grepl(searchTerm, tolower(morph1.morph2$study_specific_gene_details[t])))){
          tempIDs <- append(tempIDs, morph1.morph2$Gene_stable_ID[t])
        }
      }  
      
      for(t in 1:nrow(condition_control)){
        # If current entry is part of the gene name, gene description, or Ensembl
        # description, add the ID to the vector of IDs
        if((grepl(searchTerm, tolower(condition_control$Gene_name[t]))) |
           (grepl(searchTerm, tolower(condition_control$Gene_description[t]))) |
           (grepl(searchTerm, tolower(condition_control$Ensembl_Family_Description[t])))){
          tempIDs <- append(tempIDs, condition_control$Gene_stable_ID[t])
        }
      }
      
      # GO table
      for(go in 1:nrow(GeneToGO)){
        # If current entry is found in the gene name or any of the GO info, add
        # the ID corresponding to this information to the vector of gene IDs
        if((grepl(searchTerm, tolower(GeneToGO$Gene.names[go]))) |
           (grepl(searchTerm, tolower(GeneToGO$Gene.ontology..biological.process.[go]))) |
           (grepl(searchTerm, tolower(GeneToGO$Gene.ontology..cellular.component.[go]))) |
           (grepl(searchTerm, tolower(GeneToGO$Gene.ontology..molecular.function.[go]))) |
           (grepl(searchTerm, tolower(GeneToGO$Gene.ontology.IDs[go])))){
          tempIDs <- append(tempIDs, GeneToGO$Ensembl_GeneID[go])
        }
      }
      
      # If no info is found for the current ID, append a warning and skip
      if(length(tempIDs) == 0){
        warn <- append(warn, paste(c("No genes found corresponding to the input",
                                     input_vec[i]), collapse = " "))
      }else{
        # Once ID(s) corresponding to this entry are found, remove duplicate/NA IDs,
        # append the IDs to the gene ID vector, then move to the next entry
        geneIDs <- append(geneIDs, tempIDs[!duplicated(tempIDs) & 
                                             !grepl("NA", tempIDs) &
                                             !is.na(tempIDs)])
        next
      }
    }
  }
  
  # Initialize final dataframes in which to store data corresponding to ALL gene
  # IDs
  finalPos <- data.frame(matrix(nrow = length(geneIDs), ncol = 6))
  finalGO <- data.frame(matrix(nrow = length(geneIDs), ncol = 7))
  # Do not give row names for transcription or popgen data because there could
  # be multiple rows per gene
  finalTransc <- data.frame(matrix(ncol = 11))
  finalPopgen <- data.frame(matrix(ncol = 7))
  
  names(finalPos) <- c("Gene ID",
                       "Gene name",
                       "Scaffold",
                       "Start Locus",
                       "End Locus",
                       "Publication")
  names(finalTransc) <- c("Gene ID",
                          "Gene name",
                          "Gene description",
                          "Study-specific information",
                          "Comparison",
                          "logFC",
                          "p-value",
                          "Condition",
                          "Age at Sampling",
                          "Tissue",
                          "Publication")
  names(finalPopgen) <- c("Gene ID",
                          "Gene name",
                          "Gene description",
                          "Statistic Type",
                          "Population(s)",
                          "Statistic Value",
                          "Publication")
  names(finalGO) <- c("Gene ID",
                      "Gene name",
                      "Gene Ontology IDs",
                      "Biological Process",
                      "Cellular Component",
                      "Molecular Function",
                      "Publication")
  # If no info was found for ANY of the entries, output an ERROR
  if(length(geneIDs) == 0){
    warn <- paste(c("Error: No genes can be described by any of these inputs:", 
                    input), collapse = "\n")
    return(list(finalPos, finalTransc, finalPopgen, tempGO, warn))
  }
  # Iterate through each gene ID in the vector of gene IDs
  for(i in 1:length(geneIDs)){
    # If position data was requested...
    if(posBool){
      # Check if this gene ID is present in the position table
      if(geneIDs[i] %in% position_table$Gene_ID){
        # If so, output the gene ID's info to the temp dataframe and bind this
        # temp dataframe
        finalPos$`Gene ID`[i] <- position_table$Gene_ID[
          position_table$Gene_ID == geneIDs[i]]
        finalPos$`Gene name`[i] <- position_table$Gene_Name[
          position_table$Gene_ID == geneIDs[i]]
        finalPos$Scaffold[i] <- position_table$Scaffold[
          position_table$Gene_ID == geneIDs[i]]
        finalPos$`Start Locus`[i] <- position_table$Start_Locus[
          position_table$Gene_ID == geneIDs[i]]
        finalPos$`End Locus`[i] <- position_table$End_Locus[
          position_table$Gene_ID == geneIDs[i]]
        finalPos$Publication[i] <- "6"
        
      }else{
        # If not, output all NAs at the corresponding row in the dataframe so the 
        # length of the dataframe stays consistent. Remove these rows later
        finalPos[i,] <- rep(NA, ncol(finalPos))
        
        # Output a warning saying that position data is not present for this 
        # ID
        warn <- append(warn, paste(c("Positional data not present for gene ID ", 
                                     geneIDs[i]), collapse = ""))
        
      }
    }
    
    # If transcription data was requested...
    if(transcBool){
      # Start with the assumption that there is no transcription data present
      # for this gene
      transcAbsent <- T
      if(geneIDs[i] %in% morph1.morph2$Gene_stable_ID){
        # If there is transcription data for morph-morph comparisons, add this
        # data to a temporary dataframe and bind the temp dataframe to the final
        transcAbsent <- F
        tempTransc <- data.frame(matrix(
          nrow = sum(morph1.morph2$Gene_stable_ID == geneIDs[i]), ncol = 11))
        names(tempTransc) <- c("Gene ID","Gene name","Gene description",
                               "Study-specific information","Comparison",
                               "logFC","p-value","Condition","Age at Sampling",
                               "Tissue","Publication")
        # Since there may be multiple rows corresponding to a single gene ID,
        # output as many gene IDs as there are copies
        tempTransc$`Gene ID` <- morph1.morph2$Gene_stable_ID[
          morph1.morph2$Gene_stable_ID == geneIDs[i]
        ]
        tempTransc$`Gene name` <- morph1.morph2$Gene_name[
          morph1.morph2$Gene_stable_ID == geneIDs[i]
        ]
        tempTransc$`Gene description` <- "Not available for this study"
        tempTransc$`Study-specific information` <- morph1.morph2$study_specific_gene_details[
          morph1.morph2$Gene_stable_ID == geneIDs[i]
        ]
        tempTransc$Comparison <- morph1.morph2$Comparison[
          morph1.morph2$Gene_stable_ID == geneIDs[i]
        ]
        tempTransc$logFC <- morph1.morph2$logFC[
          morph1.morph2$Gene_stable_ID == geneIDs[i]
        ]
        tempTransc$`p-value` <- morph1.morph2$PValue[
          morph1.morph2$Gene_stable_ID == geneIDs[i]
        ]
        tempTransc$Condition <- morph1.morph2$Condition[
          morph1.morph2$Gene_stable_ID == geneIDs[i]
        ]
        tempTransc$`Age at Sampling` <- morph1.morph2$Age_at_Sampling[
          morph1.morph2$Gene_stable_ID == geneIDs[i]
        ]
        tempTransc$Tissue <- morph1.morph2$Tissue[
          morph1.morph2$Gene_stable_ID == geneIDs[i]
        ]
        tempTransc$Publication <- morph1.morph2$Publication[
          morph1.morph2$Gene_stable_ID == geneIDs[i]
        ]
        
        finalTransc <- rbind(finalTransc, tempTransc)
        
      }
      
      if(geneIDs[i] %in% condition_control$Gene_stable_ID){
        # If there is transcription data for condition-control comparisons, add 
        # this data to a temporary dataframe and bind the temp dataframe to the 
        # final
        transcAbsent <- F
        tempTransc <- data.frame(matrix(
          nrow = sum(condition_control$Gene_stable_ID == geneIDs[i]), ncol = 11))
        names(tempTransc) <- c("Gene ID","Gene name","Gene description",
                               "Study-specific information","Comparison",
                               "logFC","p-value","Condition","Age at Sampling",
                               "Tissue","Publication")
        # Since there may be multiple rows corresponding to a single gene ID,
        # output as many gene IDs as there are copies
        tempTransc$`Gene ID` <- condition_control$Gene_stable_ID[
          condition_control$Gene_stable_ID == geneIDs[i]
        ]
        tempTransc$`Gene name` <- condition_control$Gene_name[
          condition_control$Gene_stable_ID == geneIDs[i]
        ]
        tempTransc$`Gene description` <- condition_control$Gene_description[
          condition_control$Gene_stable_ID == geneIDs[i]
        ]
        tempTransc$`Study-specific information` <- "Not available for this study"
        tempTransc$Comparison <- condition_control$Comparison[
          condition_control$Gene_stable_ID == geneIDs[i]
        ]
        tempTransc$logFC <- condition_control$logFC[
          condition_control$Gene_stable_ID == geneIDs[i]
        ]
        tempTransc$`p-value` <- condition_control$PValue[
          condition_control$Gene_stable_ID == geneIDs[i]
        ]
        tempTransc$Condition <- condition_control$Condition[
          condition_control$Gene_stable_ID == geneIDs[i]
        ]
        tempTransc$`Age at Sampling` <- condition_control$Age_at_Sampling[
          condition_control$Gene_stable_ID == geneIDs[i]
        ]
        tempTransc$Tissue <- condition_control$Tissue[
          condition_control$Gene_stable_ID == geneIDs[i]
        ]
        tempTransc$Publication <- condition_control$Publication[
          condition_control$Gene_stable_ID == geneIDs[i]
        ]
        
        finalTransc <- rbind(finalTransc, tempTransc)
        
      }
      # If there is no transcription data for this gene in either class, output
      # a warning
      if(transcAbsent){
        warn <- append(warn, paste(c("Transcriptional data not present for gene ID ", 
                                     geneIDs[i]), collapse = ""))
      }
    }
    
    # If popgen data was requested...
    if(popgenBool){
      if(geneIDs[i] %in% stat_table$Stable_Gene_ID){
        # Create a subset of the stat_table dataframe
        subsetPopgen <- subset(stat_table, subset = Stable_Gene_ID == geneIDs[i])
        
        # For each row in the subset dataframe, search each STATISTIC column and 
        # dissect the column information into a row of the temporary dataframe
        for(r in 1:nrow(subsetPopgen)){
          for(c in 1:ncol(stat_table)){
            # First, check that this is a statistic column
            if(typeof(stat_table[,c]) == "double"){
              # Ensure the current statistic value is not simply NA
              if(!is.na(subsetPopgen[r,c])){
                # Create a temporary dataframe with just a single row
                tempPopgen <- data.frame(matrix(nrow = 1, ncol = 7))
                names(tempPopgen) <- c("Gene ID","Gene name","Gene description",
                                       "Statistic Type","Population(s)","Statistic Value",
                                       "Publication")
                tempPopgen$`Gene ID`[1] <- geneIDs[i]
                tempPopgen$`Gene name`[1] <- subsetPopgen$Gene_Name[r]
                tempPopgen$`Gene description`[1] <- subsetPopgen$Gene_Description[r]
                names_pops <- str_split(names(subsetPopgen)[c], "_")[[1]]
                tempPopgen$`Statistic Type`[1] <- names_pops[1]
                tempPopgen$`Population(s)`[1] <- names_pops[2]
                tempPopgen$`Statistic Value`[1] <- subsetPopgen[r,c]
                
                tempPopgen$Publication[1] <- subsetPopgen$Publication[r]
                
                finalPopgen <- rbind(finalPopgen, tempPopgen)
              }
            }
          }
        }
        
        # Bind the temporary dataframe to the final
      }else{
        # If the gene ID is not present in the statistical data, output a warning
        warn <- append(warn, paste(c("Population genetics data not present for gene ID ", 
                                     geneIDs[i]), collapse = ""))
      }
    }
    
    # If GO data was requested...
    if(GOBool){
      # Check if this gene ID is present in the position table
      if(geneIDs[i] %in% GeneToGO$Ensembl_GeneID){
        # If so, output the gene ID's GO info to the final dataframe
        finalGO$`Gene ID`[i] <- geneIDs[i]
        finalGO$`Gene name`[i] <- GeneToGO$Gene.names[
          GeneToGO$Ensembl_GeneID == geneIDs[i]]
        finalGO$`Gene Ontology IDs`[i] <- GeneToGO$Gene.ontology.IDs[
          GeneToGO$Ensembl_GeneID == geneIDs[i]]
        finalGO$`Biological Process`[i] <- GeneToGO$Gene.ontology..biological.process.[
          GeneToGO$Ensembl_GeneID == geneIDs[i]]
        finalGO$`Cellular Component`[i] <- GeneToGO$Gene.ontology..cellular.component.[
          GeneToGO$Ensembl_GeneID == geneIDs[i]]
        finalGO$`Molecular Function`[i] <- GeneToGO$Gene.ontology..molecular.function.[
          GeneToGO$Ensembl_GeneID == geneIDs[i]]
        finalGO$Publication[i] <- "7"
        
      }else{
        # If not, output all NAs at the corresponding row in the dataframe so the 
        # length of the dataframe stays consistent. Remove these rows later
        finalGO[i,] <- rep(NA, ncol(finalGO))
        
        # Output a warning saying that GO data is not present for this 
        # ID
        warn <- append(warn, paste(c("GO data not present for gene ID ", 
                                     geneIDs[i]), collapse = ""))
        
      }
    }
  }
  
  # Once you have each of the dataframes, remove all NA rows from each dataframe
  # which is supposed to have and output, then output all 4 dataframes AND 
  # warnings as a list
  if(posBool){
    finalPos <- finalPos[!is.na(finalPos$`Gene ID`),]
  }
  if(transcBool){
    finalTransc <- finalTransc[!is.na(finalTransc$`Gene ID`),]
  }
  if(popgenBool){
    finalPopgen <- finalPopgen[!is.na(finalPopgen$`Gene ID`),]
  }
  if(GOBool){
    finalGO <- finalGO[!is.na(finalGO$`Gene ID`),]
  }
  
  return(list(finalPos,finalTransc,finalPopgen,finalGO,warn))
}

TranscTable <- function(morph1, morph2, condition, direction, tr.stat, tr.thresh, GOTable){
  wrnings <- c("Errors: ")
  # If condition is NOT "Between morph"...
  if(condition != "Between morph"){
    # Use transcription data of morph-control comparisons
    in_table <- condition_control
    # Set morph2 to control
    morph2 <- "Control"
    # Set grep pattern to morph1
    search_cond <- morph1
  # If searching between morphs, set search condition to morphs-of-comparison
  }else if(condition == "Between morph"){
    # Use transcription data for between-morph comparisons
    in_table <- morph1.morph2
    # Set grep pattern to comparison
    search_cond <- paste(c(morph1, morph2), collapse = "-")
  }
  # Store comparison for output. Note: "morph2" is control if "Between morph"
  # was not entered
  comp <- paste(c(morph1, morph2), collapse = "-")

  # If genes above threshold were requested...
  if(direction == "Above"){
    # Find all rows-of-interest (ROIs) for morph(s)-of-interest where logFC is
    # above specified value, and condition matches the input specification
    if(tr.stat == "logFC"){
      ROIs <- in_table[(grepl(search_cond, in_table$Comparison)
                        & (in_table$logFC > tr.thresh) &
                          (in_table$Condition == condition)), ]
      # Sort candidate rows with highest logFC values on top
      ROIs <- ROIs[order(ROIs$logFC, decreasing = T),]
    # Find all rows-of-interest (ROIs) for morph(s)-of-interest where logFC is
    # below specified value, and condition matches the input specification
    }else{
      ROIs <- in_table[(grepl(search_cond, in_table$Comparison)
                        & (in_table$PValue > tr.thresh) &
                          (in_table$Condition == condition)), ]
      # Sort candidate rows with highest p values on top
      ROIs <- ROIs[order(ROIs$PValue, decreasing = T),]
    }
    
  # If genes whose value is below the provided stat were requested...
  }else if(direction == "Below"){
    # Find all rows-of-interest (ROIs) for morph(s)-of-interest where logFC is
    # BELOW specified value, and condition matches the input specification
    if(tr.stat == "logFC"){
      ROIs <- in_table[(grepl(search_cond, in_table$Comparison)
                        & (in_table$logFC < tr.thresh) &
                          (in_table$Condition == condition)), ]
      # Sort candidate rows with SMALLEST logFC values on top
      ROIs <- ROIs[order(ROIs$logFC, decreasing = F),]
    }else{
      ROIs <- in_table[(grepl(search_cond, in_table$Comparison)
                        & (in_table$PValue < tr.thresh) &
                          (in_table$Condition == condition)), ]
      ROIs <- ROIs[order(ROIs$PValue, decreasing = F),]
    }
  }
  # Check if any genes were found for specified conditions
  if(nrow(ROIs) == 0){
    wrnings <- append(wrnings, "No genes found matching given parameters.")
    output.df <- as.data.frame(matrix(rep(NA,10), ncol = 10))
    names(output.df) <- c("Gene Name","Gene Stable ID","GO Term(s)","Comparison",
      "logFC","p-value","Age at Sampling","Tissue","Ensembl Family Description",
      "Publication"
    )
    return(list(output.df, wrnings))
  }
  # Obtain GO terms for ROIs
  GOTerms <- character(length = nrow(ROIs))
  for(i in 1:length(GOTerms)){
    if(length(grep(ROIs$Gene_stable_ID[i], GOTable$Ensembl)) != 0){
      GOTerms[i] = paste(GOTable$Gene.ontology.IDs[GOTable$Ensembl == ROIs$Gene_stable_ID[i]], collapse = " ")
    }else{
      GOTerms[i] = NA
    }
  }
  # Output gene names, gene stable IDs, GO terms, morph-of-comparison, logFC,
  # p-value, Ensembl family information (or study specific details, if Between
  # morph), and publication name to a dataframe
  if(condition == "Between morph"){
    special_info <- ROIs$study_specific_gene_details
    special_name <- "Study-Specific Details"
  }else{
    special_info <- ROIs$Ensembl_Family_Description
    special_name <- "Ensembl Family Description"
  }
  output.df <- data.frame(
    tolower(ROIs$Gene_name),
    ROIs$Gene_stable_ID,
    GOTerms,
    rep(comp, nrow(ROIs)),
    ROIs$logFC,
    ROIs$PValue,
    ROIs$Age_at_Sampling,
    ROIs$Tissue,
    special_info,
    ROIs$Publication
  )
  names(output.df) <- c(
    "Gene Name",
    "Gene Stable ID",
    "GO Term(s)",
    "Comparison",
    "logFC",
    "p-value",
    "Age at Sampling",
    "Tissue",
    special_name,
    "Publication"
  )
  return(list(output.df,wrnings))
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
            "Publication"
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
        temp_stat_table <- stat_table[!is.na(stat_table[,indices[i]]),]
        # For each index, collect all genes, scaffolds, populations, values,
        # and publication names whose stat values fall above the entered value
        genes <- append(genes,temp_stat_table$Gene_Name[
          temp_stat_table[,indices[i]] >= thresh])
        DF_pops <- append(DF_pops,rep(pop_strings[i],
                                      length(temp_stat_table$Gene_Name[
                                        temp_stat_table[,indices[i]] >= thresh])))
        stat_vals <- append(stat_vals,temp_stat_table[
          temp_stat_table[,indices[i]] >= thresh,indices[i]])
        pub_names <- append(pub_names,temp_stat_table$Publication[
          temp_stat_table[,indices[i]] >= thresh])
      }
      # If lower tail was requested, iterate through each index
    }else if(UL == "bottom"){
      for(i in 1:length(indices)){
        # First, remove all NA values for this index
        temp_stat_table <- stat_table[!is.na(stat_table[,indices[i]]),]
        # For each index, collect all genes, populations, and statistic values
        # whose values fall in lowest tail
        genes <- append(genes,temp_stat_table$Gene_Name[
          temp_stat_table[,indices[i]] <= thresh])
        DF_pops <- append(DF_pops,rep(pop_strings[i],
                                      length(temp_stat_table$Gene_Name[
                                        temp_stat_table[,indices[i]] <= thresh])))
        stat_vals <- append(stat_vals,temp_stat_table[
          temp_stat_table[,indices[i]] <= thresh,indices[i]])
        pub_names <- append(pub_names,temp_stat_table$Publication[
          temp_stat_table[,indices[i]] <= thresh])
      }
    }
    # If count was entered...
  }else if(in_type == "Gene Count"){
    # Check whether top or bottom proportion was requested
    # If top count was requested, iterate through each index
    if(UL == "top"){
      # For each index...
      for(i in 1:length(indices)){
          # First, remove all NA values for this index
          temp_stat_table <- stat_table[!is.na(stat_table[,indices[i]]),]
          # Collect positions of top N genes for the current index
          top_genes <- order(
            temp_stat_table[,indices[i]], decreasing = T)[1:thresh]
          # Collect the genes with the highest values, as well as the associated
          # populations and values
          genes <- append(genes,temp_stat_table$Gene_Name[top_genes])
          DF_pops <- append(DF_pops,rep(pop_strings[i],length(top_genes)))
          stat_vals <- append(stat_vals,temp_stat_table[top_genes,indices[i]])
          pub_names <- append(pub_names,temp_stat_table$Publication[
            top_genes])
        }
      # If lower proportion was requested, iterate through each index
    }else if(UL == "bottom"){
      # For each index...
      for(i in 1:length(indices)){
        # First, remove all NA values for this index
        temp_stat_table <- stat_table[!is.na(stat_table[,indices[i]]),]
        # Collect positions of bottom N genes for the current index
        bottom_genes <- order(temp_stat_table[,indices[i]], decreasing = F)[1:thresh]
        # Collect the genes with the highest values, as well as the
        # associated populations, values, and publications
        genes <- append(genes,temp_stat_table$Gene_Name[bottom_genes])
        DF_pops <- append(DF_pops,rep(pop_strings[i],length(bottom_genes)))
        stat_vals <- append(stat_vals,temp_stat_table[bottom_genes,indices[i]])
        pub_names <- append(pub_names,temp_stat_table$Publication[bottom_genes])
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
      "Publication"
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
      DF_GOs[g] = GeneToGO$Gene.ontology.IDs[GeneToGO$Gene.names == genes[g]][1]
      # If gene is present in position table but NOT GO table, output NA for GO
      # term but output real scaffold
    }else if((genes[g] %in% position_table$Gene_Name) & !(genes[g] %in% GeneToGO$Gene.names)){
      scaffs[g] = position_table$Scaffold[position_table$Gene_Name == genes[g]]
      DF_GOs[g] = "Not applicable"
      # If gene is present in GO table but NOT position table, output NA for scaffold
      # term but output real GO
    }else if(!(genes[g] %in% position_table$Gene_Name) & (genes[g] %in% GeneToGO$Gene.names)){
      scaffs[g] = "Not applicable"
      DF_GOs[g] = GeneToGO$Gene.ontology.IDs[GeneToGO$Gene.names == genes[g]][1]
      # If gene is present in neither GO nor position tables, output NA for scaff
      # and GO
    }else if(!(genes[g] %in% position_table$Gene_Name) & !(genes[g] %in% GeneToGO$Gene.names)){
      scaffs[g] = "Not applicable"
      DF_GOs[g] = "Not applicable"
    }
  }
  # Create a dataframe of scaffolds, gene names, GO terms, statistic types, and
  # stat values.
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
    "Publication"
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
  error_df <- data.frame(
    x = c(.95, .955,  .975,  1, 1.025, 1.05,   1.07,  1.075, 0.99, 1.035),
    y = c( 1, 1.5, 1.75, 2, 2,    1.75,  1.5,   1,    3.5,    3.5)
  )
  error_plot <- ggplot(data = error_df, aes(x = x, y = y)) +
    geom_point() +
    ylab("") +
    xlab("Whoops! Something went wrong. See errors")
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

StatByChrTable <- function(GOTerm, GeneToGo, MasterGO, UpperLower,
                           stat_vec, position_table, stat_table, pops){
  # Initialize vector for all GO terms of interest
  GOs <- c()
  # Initialize vector in which to store warnings
  wrnings <- c("Notes: ")

  # Check if user inputted a word/phrase or a GO ID. If the user inputted a
  # word/phrase, find the name(s) which contains that word/phrase and set the GO
  # vector equal to the corresponding GO IDs.
  if(!(GOTerm %in% MasterGO$GO_ID)){
    if(length(MasterGO$GO_ID[which(grepl(GOTerm, MasterGO$GO_Term,
                                            ignore.case = T))]) != 0){
      GOs <- MasterGO$GO_ID[which(grepl(GOTerm, MasterGO$GO_Term,
                                           ignore.case = T))]
      # If the input is not found in the GO terms but does have a comma,
      # output both as GO IDs and check whether they are valid GO IDs
    }else if((length(MasterGO$GO_ID[which(grepl(GOTerm, MasterGO$GO_Term,
                                                   ignore.case = T))]) == 0)
             & (grepl(", ", GOTerm, ignore.case = T))){
      GOs <- str_split(GOTerm, ", ")[[1]]
      # If the input is not found in the GO terms and does not have a comma,
      # return an error
    }else if((length(MasterGO$GO_ID[which(grepl(GOTerm, MasterGO$GO_Term,
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
  }else if(GOTerm %in% MasterGO$GO_ID){
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
  # Find all gene IDs associated with the current GO IDs and add to vector of gene IDs
  GeneID_vec <- c()
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
      GeneID_vec <- append(GeneID_vec, GeneToGO$Ensembl_GeneID[grepl(GOs[g],
                                                             GeneToGO$Gene.ontology.IDs)])
      found_GOs <- append(found_GOs, rep(GOs[g],
                                         length(GeneToGO$Ensembl_GeneID[grepl(GOs[g],
                                                                          GeneToGO$Gene.ontology.IDs)])))
      # If the GO term does NOT appear in the dataframe of names, skip it
    }else{
      next
    }
  }
  # If no gene IDs were found corresponding to the current GO ID, return an error
  if(is.null(GeneID_vec)){
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
    return(list("No genes were found corresponding to any of the input IDs or phrase",
           null.df))
  }else{
    geneGOs <- data.frame(
      GeneID = GeneID_vec,
      GO_ID = found_GOs
    )
    geneGOs <- geneGOs[!duplicated(geneGOs), ]
  }

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
                                  all_pops[2, pair]),collapse = ""), null.df))
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
  GeneID <- c()
  Scaffold <- c()
  Start_Position <- c()
  End_Position <- c()
  GO_IDs <- c()
  Statistic_Type <- c()
  Population <- c()
  Statistic_Value <- c()
  Publication_Number <- c()

  for(i in 1:nrow(stat_pop_combos)){
    for(g in 1:length(geneGOs$Gene)){
      # Check if current gene ID is in stat AND position table
      # If so...
      if((geneGOs$GeneID[g] %in% stat_table$Stable_Gene_ID) &
         (geneGOs$GeneID[g] %in% position_table$Gene_ID)){
        # Find the number of copies of the current gene ID in the statistic
        # table
        copies <- length(stat_table[stat_table$Stable_Gene_ID == geneGOs$GeneID[g],
                                    as.numeric(stat_pop_combos$Col[i])])
        # Output the current gene as many times as there are copies
        GeneID <- append(GeneID, rep(geneGOs$GeneID[g], copies))
        # Output scaffold of current gene as many times as there are copies
        Scaffold <- append(Scaffold,
                           rep(position_table$Scaffold[position_table$Gene_ID == geneGOs$GeneID[g]],
                               copies))
        # Output starting position of the current gene as many times as
        # there are copies of the gene
        Start_Position <- append(Start_Position,
                                 rep(position_table$Start_Locus[position_table$Gene_ID == geneGOs$GeneID[g]],
                                     copies))
        # Output the ending position of the current gene as many times as
        # there are copies of the gene
        End_Position <- append(End_Position,
                               rep(position_table$End_Locus[position_table$Gene_ID == geneGOs$GeneID[g]],
                                   copies))
        # Output ALL GO terms associated with the current gene ID as many times
        # as there are copies of the gene
        GO_IDs <- append(GO_IDs,
                         rep(
                           paste(geneGOs$GO_ID[geneGOs$GeneID == geneGOs$GeneID[g]], collapse = "; "),
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
                                  stat_table[stat_table$Stable_Gene_ID == geneGOs$GeneID[g],
                                             as.numeric(stat_pop_combos$Col[i])])
        # Output the publication from which the statistic was obtained
        Publication_Number <- append(Publication_Number,
                           stat_table$Publication[stat_table$Stable_Gene_ID == geneGOs$GeneID[g]])
        # If not, skip the gene
      }else{
        next
      }
    }
  }
  output_df <- data.frame(GeneID,
                          Scaffold,
                          Start_Position,
                          End_Position,
                          GO_IDs,
                          Statistic_Type,
                          Population,
                          Statistic_Value,
                          Publication_Number
  )
  # Remove all values where the stat-of-interest is NA
  output_df <- output_df[!is.na(output_df$Statistic_Value), ]
  return(list(wrnings, output_df))
}
StatByChrGraph <- function(Full_Table, stat_vec){
  # Check if table is completely empty. If so, yield an error plot
  if(sum(is.na(Full_Table)) == 9){
    plist <- vector(mode = "list", length = 1)
    plist[[1]] <- plot(x = c(.95, .955,  .975,  1, 1.025, 1.05,   1.07,  1.075, 0.99, 1.035),
                       y = c( 1, 1.5, 1.75, 2, 2,    1.75,  1.5,   1,    3.5,    3.5), pch = 16, axes = F,
                       ylab = "",
                       xlab = "Whoops! Something went wrong. See errors")
    return(plist)
  }
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
        Publication_Number = Full_Table$Publication_Number[(Full_Table$Scaffold == unique_scaffs[c])
                                       & (Full_Table$Statistic_Type == unique_stats[s])]
      )

      # Plot the statistic values for all genes which occur on the scaffold
      # of interest, colored by population or population pair
      plist[[p]] <- ggplot(plot_dat, aes(x=Locus, y=Value, color=Populations)) +
        ylab(stat_vec[s]) +
        xlab("Locus") +
        ggtitle(paste(c("Scaffold:", unique_scaffs[c]), collapse = " ")) +
        geom_point(aes(shape = Publication_Number)) +
        guides(shape=guide_legend(title="Publication")) +
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
GOInfo <- function(GO_input, MasterGO, UpperLower){
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
      if(temp_GO_ID_vec[i] %in% MasterGO$GO_ID){
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

  }else if(!(grepl(", ", GO_input)) & (GO_input %in% MasterGO$GO_ID)){
    # If input is a single GO ID, add it to the vector of GO IDs
    GO_ID_vec <- GO_input
  }else if(!(grepl(", ", GO_input)) & !(GO_input %in% MasterGO$GO_ID)){
    # If input is a phrase, find all GO terms associated with that phrase
    if(sum(grepl(GO_input, MasterGO$GO_Term, ignore.case = T)) != 0){
      GO_ID_vec <- MasterGO$GO_ID[grepl(GO_input, MasterGO$GO_Term,
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
    if(GO_ID_vec[GO] %in% MasterGO$GO_ID){
      GO_df$`GO Term`[GO] <- MasterGO$GO_Term[MasterGO$GO_ID == GO_ID_vec[GO]]
    }else{
      GO_df$`GO Term`[GO] <- "Not applicable"
    }
    # Output class of GO ID
    if(GO_ID_vec[GO] %in% MasterGO$GO_ID){
      GO_df$Namespace[GO] <- MasterGO$Namespace[MasterGO$GO_ID == GO_ID_vec[GO]]
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
MinMax <- function(mm_pops, mm_stat, in_table){
  # Create matrix to store all columns which house the morphs-of-interest
  cols.w.morphs <- data.frame(matrix(nrow = length(mm_pops),
                                     ncol = ncol(in_table)))
  colnames(cols.w.morphs) <- colnames(in_table)

  # Find all columns of in_table which house at least one morph-of-interest
  for(m in 1:length(mm_pops)){
    cols.w.morphs[m,] <- grepl(mm_pops[m], names(in_table),
                               ignore.case = T)
  }
  cols.of.interest <- c()
  for(i in 1:ncol(cols.w.morphs)){
    if(sum(cols.w.morphs[,i]) > 0){
      cols.of.interest <- append(cols.of.interest, names(cols.w.morphs)[i])
    }
  }
  # Of the columns of in_table which house at least one morph-of-interest,
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
                       in_table[,names(in_table) == cols.of.interest[i]])
  }

  # Output the minimum and maximum values in that vector
  minimum <- num_pool[which.min(num_pool)]
  maximum <- num_pool[which.max(num_pool)]
  return(list(minimum,maximum))
}


