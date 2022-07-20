# Annabel Perry
# July 2022
# When it was created in Summer 2022, this script was originally intended as a 
# fill-in-the-blank style exercise for Zelun Liu & Daniel Nguyen to learn how to
# code the QTL module function in CaveCrawler.

################################### Preparation ################################
# TODO: Throughout this script, there are comments marked with "TODO". These 
#       comments describe portions of code which need to be added (aka 
#       "pseudocode"). After the code described in a "TODO" comment has been 
#       successfully added and troubleshooted, please erase the "TODO" part of 
#       the comment, but don't erase the whole comment! Remember that you will 
#       have to describe this code to a NEW undergrad at some point (maybe 1-3 
#       years from now!). So, rephrase the comment in a manner which will help 
#       describe the code to an unfamiliar user. If a comment does NOT start 
#       with "TODO", there is no need to mess with it.

# TODO P: As you code the sub-modules AND the app, add any additional libraries
# here.
library(shiny)
library(shinyWidgets)
library(ggplot2)

# The QTL function takes the following inputs: 
chr_table <- read.csv("data/ChrTable.csv", fill = TRUE)
position_table <- read.csv("data/PositionTable.csv", fill = TRUE)
QTL_table <- read.csv("data/QTL.csv", fill = TRUE)

# TODO P: Read in the appropriate files for chr_table, position_table, and QTL_table

######################## QTL Function Input Descriptions #######################
# The next inputs are specific to sub-modules and will ONLY be assigned non-NA
# values if the user has activated the corresponding sub-module through the UI.

# Genomic Range Sub-Module (The portion of code activated when GR.bool == T):
#   GR.bool = Logical indicating whether the user would like to see markers/genes
#             occurring within a bp range on a scaffold.
#   GR.scaff = Character. The scaffold to be searched.
#   GR.start = Numeric. Starting bp to be searched on the scaffold.
#   GR.end = Numeric. Ending bp to be searched on the scaffold.

# Marker Range Sub-Module (The portion of code activated when Range.bool == T):
#   MR.bool = Logical indicating whether the user would like to see all genes
#             and markers within range of a central marker.
#   MR.search_term = Character. Name of the central marker.
#   MR.bp = Numeric. Number of bp up- and downstream of marker.

# Trait-to-Marker Sub-Module (The portion of code activated when TM.bool == T):
#   TM.bool = Logical indicating whether the user would like to see all markers
#             associated with a given trait.
#   TM.QT = Character. The quantitative trait of interest

######################### QTL Function Code Start ##############################
QTL <- function(chr_table, position_table, QTL_table, 
                GR.bool, GR.chr, GR.start, GR.end, 
                MR.bool, MR.search_term, MR.bp, 
                TM.bool, TM.QT){
  # Vector into which warnings will be appended
  QTL.wrnings <- c("Notes: ")
  # Dataframe for QTL marker(s) matching the user's search parameters
  QTL_Marker_Data <- data.frame(matrix(ncol = 9))
  colnames(QTL_Marker_Data) <- c("Marker (Peak)", "Scaffold", 
                                "Start Position on Scaffold", 
                                "End Position on Scaffold", "LOD (Peak)",
                                "Quantitative Trait", 
                                "Percent Variance Explained",
                                "Study-Specific Information", "Publication(s)")
  # Dataframe for QTL marker(s) matching the user's search parameters
  QTL_Gene_Data <- data.frame(matrix(ncol = 5))
  colnames(QTL_Gene_Data) <- c("Gene", "Scaffold", "Start Locus", "End Locus",
                               "Publication")
  # Object into which plots will be stored (This is currently a placeholder)
  chr_plot <- ggplot()
  
########################### Marker Range Pseudocode ###########################
  # TODO MR: Write & troubleshoot code which performs ALL the logic described in
  #          the following several lines of comments:
  #
  #          (Note: the number of lines indicated in the pseudocode does NOT
  #          necessarily match the number of lines you'll need for the actual
  #          code. I used as many lines as needed to communicate the logic. You
  #          may need more or fewer lines depending on how you choose to code)
  #
  # If MR.bool == T...
  #   Check whether MR.search_term occurs within the "Marker (Peak)" column of 
  #   QTL_table. If so...
  #     Use QTL_table to find the chromosome which the marker in MR.search_term 
  #     occurs on and the position of MR.search_term on this chromosome
  #       Set GR.chr equal to this chromosome
  #       Set GR.start equal to this position MINUS MR.bp
  #       Set GR.end equal to this position PLUS MR.bp
  #     Set GR.bool == T
  #     Permit code to fall through to the if-statement which checks whether
  #     GR.bool == T "
  #   If MR.search_term does NOT occur within the "Marker (Peak)" column of 
  #   QTL_table, check if it occurs in "Gene Names" column of position_table. 
  #   If so...
  #     Use position_table to find the chromosome which the gene in MR.search_term 
  #     occurs on and the position of MR.search_term on this chromosome
  #       Set GR.chr equal to this chromosome
  #       Set GR.start equal to this position MINUS MR.bp
  #       Set GR.end equal to this position PLUS MR.bp
  #     Set GR.bool == T
  #     Permit code to fall through to the if-statement which checks whether
  #     GR.bool == T "
  #   If MR.search_term does NOT occur within QTL_table OR within "Gene Names", 
  #   return a descriptive but concise (in other words, user-friendly) error
  
  # if second button is activated
  
  if (MR.bool) {
    for (i in 1:nrow(QTL_table)) {
      if (MR.search_term == QTL_table$Marker[i]) {
        GR.chr <- as.numeric(QTL_table$Chromosome[i])
        GR.start <- as.numeric(QTL_table$Start_Position_on_Chromosome[i] - MR.bp)
        GR.end <- as.numeric(QTL_table$End_Position_on_Chromosome[i] + MR.bp)
        GR.bool <- T
      }
    }
    if (GR.bool == F) {
      # TODO Edits for Zelun: There are two logical issues here
      # 1. Function will fail as soon as the script encounters a marker not found
      #   in the CSV. Think about what we discussed on Sunday
      # 2. User doesn't need a warning for EVERY instance where a marker in the
      #    CSV doesn't match the search parameters
      QTL.wrnings <- append(QTL.wrnings, paste("No marker is not found in data"))
      return(list(chr_plot,QTL_Marker_Data,QTL_Gene_Data,QTL.wrnings))
    }
  }
  
########################### Genomic Range Pseudocode ###########################
  # TODO GR: Write & troubleshoot code which performs ALL the logic described in
  #          the following several lines of comments:
  #
  #          (Note: the number of lines indicated in the pseudocode does NOT
  #          necessarily match the number of lines you'll need for the actual
  #          code. I used as many lines as needed to communicate the logic. You
  #          may need more or fewer lines depending on how you choose to code)
  #
  # If GR.bool == T...
  #   Search QTL_table for all markers where...
  #     1. Chromosome name can be coerced to numeric AND Chromosome == GR.chr
  #     2. Start Locus can be coerced to numeric AND Start Locus >= GR.start
  #     3. End Locus can be coerced to numeric AND End Locus <= GR.end
  #
  #     Add each marker matching these parameters to QTL_Marker_Data
  #   Search position_table for all genes where...
  #     1. Scaffold can be coerced to numeric AND Scaffold == GR.scaffold
	#	    2. Start Locus can be coerced to numeric AND Start Locus >= GR.start
	#	    3. End Locus can be coerced to numeric AND End Locus <= GR.end
	#	  ... and add the Gene Name, Scaffold, Start Locus, and End Locus for the 
  #   corresponding gene in the position_table to the corresponding columns of 
  #   QTL_Gene_Data
  #
	#   After all genes and markers have been added to QTL_Gene_Data and to 
  #   QTL_Marker_Data, check for duplicate rows in either table
  #
	#   Output chr_plot (currently empty) + QTL_Marker_Data + QTL_Gene_Data + 
  #   warnings as list
  #
  # ... close the "if" statement here.
    
  if (GR.bool) {
    # for each row in QTL table
    GR_markers <- c() 
    for (i in 1:nrow(QTL_table)) {
      if((grepl(GR.chr,QTL_table$Chromosome[i])) &
         (QTL_table$Start_Position_on_Chromosome[i] >= GR.start) &
         (QTL_table$End_Position_on_Chromosome[i] <= GR.end)){
        GR_markers <- append(GR_markers, QTL_table$Marker[i])
      }
    }
    QTL_Marker_Data <- data.frame(matrix(nrow = length(GR_markers),ncol = 9))
    colnames(QTL_Marker_Data) <- c("Marker (Peak)", "Scaffold", 
                                   "Start Position on Scaffold", 
                                   "End Position on Scaffold", "LOD (Peak)",
                                   "Quantitative Trait", 
                                   "Percent Variance Explained",
                                   "Study-Specific Information", "Publication(s)")
    j <- 1
    for(i in 1:nrow(QTL_table)){
      # You are resetting 'j' to 1 each time you enter a new row of the
      # QTL_table. This is wrong. Why?
      # Zelun: We want to see if after traversing through the loop, no marker is found,
      # so we want to set j <- 1 before entering the loop
      #if chromosome at that row is same as user input
      if ((grepl(GR.chr,QTL_table$Chromosome[i])) &
          (QTL_table$Start_Position_on_Chromosome[i] >= as.numeric(GR.start)) &
          (QTL_table$End_Position_on_Chromosome[i] <= as.numeric(GR.end))) {
        # TODO Edits for Zelun
        # As this is currently written, asigning a value to row 'j' of 
        # QTL_Marker_Data will fail if j > 1.
        #  To Troubleshoor:
        #   1. Set j = 2
        #   2. Run lines 59-65 to define QTL_Marker_Data
        #   3. Run this code: QTL_Marker_Data$`Marker (Peak)`[j] <- "enter anything here"
        #   You'll get an error. Troubleshoot it.
        #output this marker's info to the marker data frame, with j keeping track of order
        QTL_Marker_Data$`Marker (Peak)`[j] <- QTL_table$Marker[i]
        QTL_Marker_Data$`Scaffold`[j] <- QTL_table$Chromosome[i]
        QTL_Marker_Data$`Start Position on Scaffold`[j] <- QTL_table$Start_Position_on_Chromosome[i]
        QTL_Marker_Data$`End Position on Scaffold`[j] <- QTL_table$End_Position_on_Chromosome[i]
        QTL_Marker_Data$`LOD (Peak)`[j] <- QTL_table$LOD[i]
        QTL_Marker_Data$`Quantitative Trait`[j] <- QTL_table$Quantitative_Trait[i]
        QTL_Marker_Data$`Percent Variance Explained`[j] <- QTL_table$Percent_Variance_Explained[i]
        QTL_Marker_Data$`Study-Specific Information`[j] <- QTL_table$Study_Specific_Information[i]
        QTL_Marker_Data$`Publication(s)`[j] <- QTL_table$Publication[i]
        j <- j + 1
      }
    }
    if (j == 1) {
      QTL.wrnings <- append(QTL.wrnings, paste(c("No marker is not found in selected range: ", "from ", GR.start, " to ",GR.end, " on chromosome ", GR.chr), collapse = ""))
    }
    # Now for each item in position table
    GR_genes <- c() 
    for (i in 1:nrow(position_table)) {
      if((grepl(GR.chr,position_table$Scaffold[i])) &
         (position_table$Start_Locus[i] >= GR.start) &
         (position_table$End_Locus[i] <= GR.end)){
        GR_genes <- append(GR_genes, position_table$Gene_Name[i])
      }
    }
    QTL_Gene_Data <- data.frame(matrix(nrow = length(GR_genes),ncol = 5))
    colnames(QTL_Gene_Data) <- c("Gene", "Scaffold", "Start Locus", "End Locus",
                                 "Publication")
    j <- 1
    for(i in 1:nrow(position_table)){
      # if gene position is in range
      if ((grepl(GR.chr,position_table$Scaffold[i])) &
          (position_table$Start_Locus[i] >= as.numeric(GR.start)) &
          (position_table$End_Locus[i] <= as.numeric(GR.end))) {
        # add row in QTL_Gene_Data
        # TODO Edits for Zelun: See edit on line 173
        QTL_Gene_Data$`Gene`[j] <- position_table$Gene_ID[i]
        QTL_Gene_Data$`Scaffold`[j] <- position_table$Scaffold[i]
        QTL_Gene_Data$`Start Locus`[j] <- position_table$Start_Locus[i]
        QTL_Gene_Data$`End Locus`[j] <- position_table$End_Locus[i]
        # TODO Edits for Zelun: 
        # How are Publications formatted in all the other tables on CaveCrawler?
        QTL_Gene_Data$`Publication`[j] <- as.integer(6)
        j <- j + 1
      }
    }
    if (j == 1) {
      QTL.wrnings <- append(QTL.wrnings, paste(c("No gene is not found in selected range: ", "from ", GR.start, " to ",GR.end, " on chromosome ", GR.chr), collapse = ""))
    }
    #   Output chr_plot (currently empty) + QTL_Marker_Data + QTL_Gene_Data + 
    #   warnings as list
    return(list(chr_plot,QTL_Marker_Data,QTL_Gene_Data,QTL.wrnings))
  }
  
########################### Trait-to-Marker Pseudocode #########################
  # TODO TM: Write & troubleshoot code which performs ALL the logic described in
  #          the following several lines of comments:
  #
  #          (Note: the number of lines indicated in the pseudocode does NOT
  #          necessarily match the number of lines you'll need for the actual
  #          code. I used as many lines as needed to communicate the logic. You
  #          may need more or fewer lines depending on how you choose to code)
  #
  # If QT.bool == T...
  #   Search QTL_table for all markers which are correlated with the 
  #   quantitative trait stored in TM.QT
  #     Store these markers in QTL_Marker_Data
  #   After all markers have been added to QTL_Marker_Data, check for any 
  #   duplicate markers in this table
  #   Output QTL_Marker_Data + warnings as list
  
  # TODO Edits for Daniel
  # You only want to collect markers corresponding to a quantitative trait
  # if the Trait-to-Marker sub-module has been selected, right? So... which if
  # statement needs to be passed before running any of this code?
  
  #Define a vector in which to store markers
  TM_markers <- c()
  
  #for each trait
  for(i in 1:length(TM.QT)){
    #for each row in qtl_table, of trait matches, store marker. W/O case sensitivity.
    for(j in 1:nrow(QTL_table)){
      # TODO Edits for Daniel
      # This isn't an edit but I just wanted to say nicely done!! You did a 
      # great job anticipating that sometimes the case of the inputted quantitative
      # triat might not match the case of the quantitative trait in the table.
      # You're good at thinking these things through - you're gonna go far :)
      if(grepl(tolower(TM.QT[i]), tolower(QTL_table$Quantitative_Trait[j]))){
        TM_markers <- append(TM_markers, QTL_table$Marker[j])
      }
    }
  }
  
  
  
  #check for duplicates.
  TM_markers <- TM_markers[!duplicated(TM_markers) & !grepl("NA", TM_markers) &
                             !is.na(TM_markers)]
  
  #final marker data table
  for(i in 1:length(TM_markers)){
    #if TM was requested...
    # TODO Edits for Daniel
    # Putting this if statement here is not the most elegant solution becase
    # you end up evaluating the bool for EVERY marker. See the edit on Line 159
    if(TM.bool){
      #check if marker is present in position table
      # TODO Edits for Daniel
      # Where did you get the markers in "TM_markers" from? Given this, is this
      # if statement necessary..?
      if(TM_markers[i] %in% QTL_table$Marker){
        
        temp_TM <- data.frame(matrix(
          nrow = length(TM_markers), ncol = 9))
        names(temp_TM) <- c("Marker (Peak)", "Scaffold",
                            "Start Position on Scaffold",
                            "End Position on Scaffold", "LOD (Peak)",
                            "Quantitative Trait",
                            "Percent Variance Explained",
                            "Study-Specific Information", "Publication(s)")
        
        
        #if so, output the marker's info to the marker data frame
        # TODO Edits for Daniel
        # This is not incorrect, but it is redundant. Where did you get the 
        # markers stored in 'TM_markers' from?
        temp_TM$`Marker (Peak)`[i] <- QTL_table$Marker[
          QTL_table$Marker == TM_markers[i]]
        temp_TM$`Scaffold`[i] <- QTL_table$Chromosome[
          QTL_table$Marker == TM_markers[i]]
        temp_TM$`Start Position on Scaffold`[i] <- QTL_table$Start_Position_on_Chromosome[
          QTL_table$Marker == TM_markers[i]]
        temp_TM$`End Position on Scaffold`[i] <- QTL_table$End_Position_on_Chromosome[
          QTL_table$Marker == TM_markers[i]]
        temp_TM$`LOD (Peak)`[i] <- QTL_table$LOD[
          QTL_table$Marker == TM_markers[i]]
        temp_TM$`Quantitative Trait`[i] <- QTL_table$Quantitative_Trait[
          QTL_table$Marker == TM_markers[i]]
        temp_TM$`Percent Variance Explained`[i] <- QTL_table$Percent_Variance_Explained[
          QTL_table$Marker == TM_markers[i]]
        temp_TM$`Study-Specific Information`[i] <- QTL_table$Study_Specific_Information[
          QTL_table$Marker == TM_markers[i]]
        # TODO Edits for Daniel
        # How are Publications formatted in the rest of the tables on this website?
        # You'll want the formatting to be consistent with the rest of the tables
        temp_TM$`Publication(s)`[i] <- QTL_table$Publication[
          QTL_table$Marker == TM_markers[i]]
        
        QTL_Marker_Data <- rbind(QTL_Marker_Data, temp_TM)
        
      }else{
        #if marker is NOT present in position table...
        #output NAs at corresponding row, length data frame stays consistent.
        #remove these rows later
        QTL_Marker_Data[i,] <- rep(NA, ncol(QTL_Marker_Data))
        
        #output a warning saying that marker data is not present for this trait
        QTL.wrnings <- append(QTL.wrnings, paste(c("Marker datta not present for trait", TM_markers[i], collapse = "")))
      }
    }
  }
  
  # remove all NA rows from each data frame which is supposed to have output
  # then output data frame and warnings as a list
  
  # TODO Edits for Daniel
  # Remember the edit on Line 159? Given that... where would be the most 
  # elegant place to move this if statement?
  if(TM.bool){
    QTL_Marker_Data <- QTL_Marker_Data[!is.na(QTL_Marker_Data$`Marker (Peak)`),]
  }
  
  #return marker table and warnings as a list
  return(list(QTL_Marker_Data, QTL.wrnings))
  
}