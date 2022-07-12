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

# The QTL function takes the following inputs: 
# chr_table = Dataframe describing the length of each chromosome in the Mexican
#             tetra genome assembly
# position_table = Dataframe describing location of every gene in Mexican tetra
#                  genome
# QTL_table = Dataframe describing each trait-associated marker

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
}