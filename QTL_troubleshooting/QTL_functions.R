# Annabel Perry
# July 2022
# When it was created in Summer 2022, this script was originally intended as a 
# fill-in-the-blank style exercise for Zelun Liu & Daniel Nguyen to learn how to
# code the QTL module function in CaveCrawler.

library(shiny)
library(shinyWidgets)
library(ggplot2)
#TODO copy-and-paste this library
library(chromoMap)

# The QTL function takes the following inputs: 
chr_table <- read.csv("data/ChrTable.csv", fill = TRUE)
position_table <- read.csv("data/PositionTable.csv", fill = TRUE)
#TODO Add this publication-changing line beneath position_table reading
position_table$Publication[position_table$Publication == "Warren_et_al_2021"] <- 
  rep("6", sum(position_table$Publication == "Warren_et_al_2021"))

QTL_table <- read.csv("data/QTL.csv", fill = TRUE)
QTL_table$Publication[QTL_table$Publication == "Warren_et_al_2021"] <- rep("6", sum(QTL_table$Publication == "Warren_et_al_2021"))


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
QTL <- function(chr_table, position_table, QTL_table, 
                GR.bool, GR.chr, GR.start, GR.end, 
                MR.bool, MR.search_term, MR.bp, 
                TM.bool, TM.QT){
  # Vector into which warnings will be appended
  QTL.wrnings <- c("Notes: ")
  # Dataframe for QTL marker(s) matching the user's search parameters
  QTL_Marker_Data <- data.frame(matrix(nrow = 0, ncol = 10))
  colnames(QTL_Marker_Data) <- c("Marker (Peak)", "Scaffold", 
                                "Start Position on Scaffold", 
                                "End Position on Scaffold", 
                                "Cross",
                                "LOD (Peak)",
                                "Quantitative Trait", 
                                "Percent Variance Explained",
                                "Study-Specific Information", "Publication(s)")
  # Dataframe for QTL marker(s) matching the user's search parameters
  QTL_Gene_Data <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(QTL_Gene_Data) <- c("Gene", "Scaffold", "Start Locus", "End Locus",
                               "Publication")
  # Object into which plots will be stored (This is currently a placeholder)
  chr_plot <- ggplot()
  # Defines a dataframe to store the marker and gene locations to be plotted
  plot_annot <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(plot_annot) <- c("Element Name", "Chromosome Name", "Element Start",
                            "Element End")
  
  if (MR.bool) {
    for (i in 1:nrow(QTL_table)) {
      if (MR.search_term == QTL_table$Marker[i]) {
        GR.chr <- as.numeric(QTL_table$Chromosome[i])
        GR.start <- as.numeric(QTL_table$Start_Position_on_Chromosome[i]) - MR.bp
        GR.end <- as.numeric(QTL_table$End_Position_on_Chromosome[i]) + MR.bp
        GR.bool <- T
      }
    }
    for (i in 1:nrow(position_table)) {
      if (is.na(position_table$Gene_Name[i])) {
        next
      }
      if (MR.search_term == position_table$Gene_Name[i]) {
        GR.chr <- as.numeric(position_table$Scaffold[i])
        GR.start <- as.numeric(position_table$Start_Locus[i]) - MR.bp
        GR.end <- as.numeric(position_table$End_Locus[i]) + MR.bp
        GR.bool <- T
      }
    }
    if (GR.bool == F) {
      QTL.wrnings <- append(QTL.wrnings, paste("No marker is not found in data"))
      return(list(data.frame(),data.frame(),chr_plot,QTL.wrnings))
    }
  }

  if (GR.bool) {
    # For the marker table, retrieve the numbers corresponding to rows which 
    # match the input parameters
    marker_rows <- which((QTL_table$Chromosome == GR.chr) &
                       (QTL_table$Start_Position_on_Chromosome >= GR.start) &
                       (QTL_table$End_Position_on_Chromosome <= GR.end))
    # If no markers in the QTL table match the input parameters, output a 
    # warning and proceed to check for genes
    if(!length(marker_rows)){
      QTL.wrnings <- append(QTL.wrnings, 
                            paste(c("No markers are found in selected range: ", 
                                    "from ", GR.start, " to ",GR.end, 
                                    " on chromosome ", GR.chr), collapse = ""))
    # If you DID find rows matching the input parameters, store their data
    }else{
      # First, add as many rows to the dataframe as there are markers 
      # corresponding to the input parameters. Save the current column names,
      # because when you add the NA rows these will get messed up
      OriginalColNames <- colnames(QTL_Marker_Data)
      QTL_Marker_Data <- rbind(QTL_Marker_Data, 
                               data.frame(matrix(nrow = length(marker_rows),
                                                 ncol = ncol(QTL_Marker_Data),
                                                 dimnames = list(1:length(marker_rows),
                                                                 OriginalColNames))))
      colnames(QTL_Marker_Data) <- OriginalColNames
      # Next, fill in each row with the data from QTL_table corresponding to 
      # that row
      QTL_Marker_Data$`Marker (Peak)` <- QTL_table$Marker[marker_rows]
      QTL_Marker_Data$Scaffold <- QTL_table$Chromosome[marker_rows]
      QTL_Marker_Data$`Start Position on Scaffold` <- 
        QTL_table$Start_Position_on_Chromosome[marker_rows]
      QTL_Marker_Data$`End Position on Scaffold` <- 
        QTL_table$End_Position_on_Chromosome[marker_rows]
      QTL_Marker_Data$Cross <- QTL_table$Cross[marker_rows]
      QTL_Marker_Data$`LOD (Peak)` <- QTL_table$LOD[marker_rows]
      QTL_Marker_Data$`Quantitative Trait` <- QTL_table$Quantitative_Trait[
        marker_rows]
      QTL_Marker_Data$`Percent Variance Explained` <- 
        QTL_table$Percent_Variance_Explained[marker_rows]
      QTL_Marker_Data$`Study-Specific Information` <- 
        QTL_table$Study_Specific_Information[marker_rows]
      QTL_Marker_Data$`Publication(s)` <- QTL_table$Publication[marker_rows]
      
      # Add the name and positions of all markers to the annotation dataframe
      # First, create a temporary dataframe in which to store the markers
      temp_df <- data.frame(matrix(nrow = length(marker_rows),
                                   ncol = ncol(plot_annot)))
      colnames(temp_df) <- colnames(plot_annot)
      # Input all relevant data into the temporary dataframe
      temp_df$`Element Name` <- QTL_table$Marker[marker_rows]
      temp_df$`Chromosome Name` <- QTL_table$Chromosome[marker_rows]
      temp_df$`Element Start` <- 
        QTL_table$Start_Position_on_Chromosome[marker_rows]
      temp_df$`Element End` <- 
        QTL_table$End_Position_on_Chromosome[marker_rows]
      # Append the temporary dataframe to the plot_annot dataframe
      plot_annot <- rbind(plot_annot, temp_df)
    }
    
    # For the position table, retrieve the numbers corresponding to rows which 
    # match the input parameters
    gene_rows <- which((position_table$Scaffold == GR.chr) &
                       (position_table$Start_Locus >= GR.start) &
                       (position_table$End_Locus <= GR.end))
    # If no genes in the position table match the input parameters, output a 
    # warning
    if(!length(gene_rows)){
      QTL.wrnings <- append(QTL.wrnings, 
                            paste(c("No genes are found in selected range: ", 
                                    "from ", GR.start, " to ",GR.end, 
                                    " on chromosome ", GR.chr), collapse = ""))
      # If you DID find rows matching the input parameters, store their data
    }else{
      # First, add as many rows to the dataframe as there are genes 
      # corresponding to the input parameters. Save the current column names,
      # because when you add the NA rows these will get messed up
      OriginalColNames <- colnames(QTL_Gene_Data)
      QTL_Gene_Data <- rbind(QTL_Gene_Data,
                             data.frame(matrix(
                               nrow = length(gene_rows),
                               ncol = ncol(QTL_Gene_Data),
                               dimnames = list(1:length(gene_rows),
                                               OriginalColNames))))
      colnames(QTL_Gene_Data) <- OriginalColNames
      
      QTL_Gene_Data$Gene <- position_table$Gene_Name[gene_rows]
      QTL_Gene_Data$Scaffold <- position_table$Scaffold[gene_rows]
      QTL_Gene_Data$`Start Locus` <- position_table$Start_Locus[gene_rows]
      QTL_Gene_Data$`End Locus` <- position_table$End_Locus[gene_rows]
      QTL_Gene_Data$Publication <- position_table$Publication[gene_rows]
      
      # Add the name and positions of all genes to the annotation dataframe
      # First, create a temporary dataframe with as many rows as there are new
      # genes
      temp_df <- data.frame(matrix(nrow = length(gene_rows),
                                   ncol = ncol(plot_annot)))
      colnames(temp_df) <- colnames(plot_annot)
      # Next, add all relevant data on the genes to the temporary df
      temp_df$`Element Name` <- position_table$Gene_Name[gene_rows]
      temp_df$`Chromosome Name` <- position_table$Scaffold[gene_rows]
      temp_df$`Element Start` <- position_table$Start_Locus[gene_rows]
      temp_df$`Element End` <- position_table$End_Locus[gene_rows]
      # Now, bind this temporary dataframe to the final plotting dataframe
      plot_annot <- rbind(plot_annot, temp_df)
    }
    
    # If the dataframe containing markers to be plotted is NOT empty, plot the
    # markers
    if(nrow(plot_annot)){
      chr_plot <- chromoMap(ch.files = list(chr_table[chr_table$Chromosome == 
                                                        GR.chr,]),
                            data.files = list(plot_annot))
    }
    
    # Return tables, plot, and warnings
    return(list(QTL_Marker_Data, QTL_Gene_Data, chr_plot, QTL.wrnings))
  }
  
  if(TM.bool){
    # Define a vector in which to store markers
    TM_markers <- c()
    
    #for each trait
    for(i in 1:length(TM.QT)){
      #for each row in qtl_table, of trait matches, store marker. W/O case sensitivity.
      for(j in 1:nrow(QTL_table)){
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
      #check if marker is present in QTL table
      if(TM_markers[i] %in% QTL_table$Marker){
          
          temp_TM <- data.frame(matrix(
            nrow = length(TM_markers), ncol = 10))
          names(temp_TM) <- c("Marker (Peak)", "Scaffold",
                              "Start Position on Scaffold",
                              "End Position on Scaffold", "Cross", 
                              "LOD (Peak)",
                              "Quantitative Trait",
                              "Percent Variance Explained",
                              "Study-Specific Information", "Publication(s)")
          
          # if so, output the marker's info to the marker data frame
          temp_TM$`Marker (Peak)`[i] <- TM_markers[i]
          temp_TM$`Scaffold`[i] <- QTL_table$Chromosome[
            QTL_table$Marker == TM_markers[i]]
          temp_TM$`Start Position on Scaffold`[i] <- QTL_table$Start_Position_on_Chromosome[
            QTL_table$Marker == TM_markers[i]]
          temp_TM$`End Position on Scaffold`[i] <- QTL_table$End_Position_on_Chromosome[
            QTL_table$Marker == TM_markers[i]]
          temp_TM$Cross <- QTL_table$Cross[
            QTL_table$Marker == TM_markers[i]]
          temp_TM$`LOD (Peak)`[i] <- QTL_table$LOD[
            QTL_table$Marker == TM_markers[i]]
          temp_TM$`Quantitative Trait`[i] <- QTL_table$Quantitative_Trait[
            QTL_table$Marker == TM_markers[i]]
          temp_TM$`Percent Variance Explained`[i] <- QTL_table$Percent_Variance_Explained[
            QTL_table$Marker == TM_markers[i]]
          temp_TM$`Study-Specific Information`[i] <- QTL_table$Study_Specific_Information[
            QTL_table$Marker == TM_markers[i]]
          temp_TM$`Publication(s)`[i] <- QTL_table$Publication[
            QTL_table$Marker == TM_markers[i]]
          
          QTL_Marker_Data <- rbind(QTL_Marker_Data, temp_TM)
          
        }else{
          #if marker is NOT present in position table...
          #output NAs at corresponding row, length data frame stays consistent.
          #remove these rows later
          QTL_Marker_Data[i,] <- rep(NA, ncol(QTL_Marker_Data))
          
          #output a warning saying that marker data is not present for this trait
          QTL.wrnings <- append(QTL.wrnings, paste(c("Marker data not present for trait", TM_markers[i], collapse = "")))
        }
    }
    
    # remove all NA rows from each data frame which is supposed to have output
    QTL_Marker_Data <- QTL_Marker_Data[!is.na(QTL_Marker_Data$`Marker (Peak)`),]

    # return marker table and warnings as a list
    return(list(QTL_Marker_Data, data.frame(), chr_plot, QTL.wrnings))
  }
}
