# This script contains the plotting code for the MR and GR submodules. This code
# will need to be copy-pasted into the final QTL module
library(shiny)
# TODO load library into Z&D combined code
library(chromoMap)

setwd("C:/Users/knigh/Documents/GitHub/CaveCrawler/cavecrawler")

chr_table <- read.csv(file = "data/ChrTable.csv", fill = F)
position_table <- read.csv(file = "data/PositionTable.csv", fill = F)
QTL_table <- read.csv(file = "data/QTL.csv", fill = F)


QTL <- function(chr_table, position_table, QTL_table, 
                GR.bool = F, GR.chr = NA, GR.start = NA, GR.end = NA, 
                MR.bool = F, MR.search_term = NA, MR.bp = NA, 
                TM.bool = F, TM.QT = NA){
  # Vector into which warnings will be appended
  QTL.wrnings <- c()
  # Dataframe for QTL marker(s) matching the user's search parameters
  QTL_Marker_Data <- data.frame()
  colnames(QTL_Marker_Data) <- c("Marker (Peak)", "Scaffold", 
                                "Start Position on Scaffold", 
                                "End Position on Scaffold", "LOD (Peak)",
                                "Quantitative Trait", 
                                "Percent Variance Explained",
                                "Study-Specific Information", "Publication(s)")
  # Dataframe for QTL marker(s) matching the user's search parameters
  QTL_Gene_Data <- data.frame()
  colnames(QTL_Gene_Data) <- c("Gene", "Scaffold", "Start Locus", "End Locus",
                               "Publication")
  # TODO Remove the chr_plot placeholder 
  # TODO Copy-and-paste
  # Defines a dataframe to store the marker and gene locations to be plotted
  plot_annot <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(plot_annot) <- c("Element Name", "Chromosome Name", "Element Start",
                            "Element End")
  
  # TODO As Zelun's code fills the dataframes, store all marker and gene names
  #      + locations in plot_annot as well
  
  # TODO Copy-and-paste this plotting function into the end of the GR code
  chr_plot <- chromoMap(ch.files = list(chr_table[chr_table$Chromosome == GR.chr,]),
                        data.files = list(plot_annot))
  # TODO Output plot object
}