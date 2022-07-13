# 07-06-2022
# Annabel Perry
# This script organizes data from Warren et al 2021 into a QTL data file for
# CaveCrawler
library(stringr)
setwd("~/Summer 2022/CaveCrawler Maintenance/QTL")

# Read in "Supplemental Data 1" from Warren et al 2021
Warren2021 <- read.csv("41467_2021_21733_MOESM4_ESM.csv", skip = 1)

# Remove all markers for which there exist NO QTL for ANY studies as well as the
# columns which are wholly NA
Warren2021 <- Warren2021[!is.na(Warren2021$Protas_et_al.2008) |
                           !is.na(Warren2021$Protas_et_al.2007) |
                           !is.na(Warren2021$Yoshizawa.et.al.2012) |
                           !is.na(Warren2021$Yoshizawa.et.al.2015) |
                           !is.na(Warren2021$Kowalko.et.al.2013.PNAS) |
                           !is.na(Warren2021$Kowalko.et.al.2013.CurrentBiology) |
                           !is.na(Warren2021$O.Quin.et.al.2012), -c(24:28)]
# Create a table describing the crosses and study-specific information for each 
# study
Key <- data.frame(
  Publication = colnames(Warren2021)[14:20],
  Cross = c(rep("Pachon x Surface F2", 2), 
            rep("Pachon x Texas Surface F2 and F3",2), NA, "Tinaja x Surface F2", 
            "Pachon x Texas Surface F2"),
  Study_Specific_Information = c(
    rep(NA, 2),
    rep("In this study, quantitative traits were measured in both F2 and F3 fish", 2),
    rep(NA, 3))
  )

# Create matrix into which QTL data should be read
QTL <- data.frame(matrix(nrow = 1, ncol = 11))
colnames(QTL) <- c("Marker",	"Chromosome",	"Linkage_Group",	"Cross",	
                   "Start_Position_on_Chromosome",	"End_Position_on_Chromosome",	
                   "LOD",	"Quantitative_Trait",	"Percent_Variance_Explained",	
                   "Study_Specific_Information", "Publication")

comma <- ","
period <- '\\.'

# Iterate through each entry in Warren2021
for(r in 1:nrow(Warren2021)){
  # Count how many publications have information for this marker AND how many
  # QTL there are per publication and add the appropriate number of rows for 
  # this marker
  for(study in 14:20){
    # Check if the current marker has information for the current study
    if(!is.na(Warren2021[r,study])){
      # Check if current marker has only one quantitative trait
      if(!grepl(comma, Warren2021[r,study]) & !grepl(period, Warren2021[r,study])){
        # Create new row in QTL matrix with...
        QTL <- rbind(QTL, c(
          # The current marker
          Warren2021$Marker.ID[r], 
          # The Chromosome for the current marker
          Warren2021$Scaffold.ID[r],
          # The linkage group for the current marker
          Warren2021$LG[r], 
          # The fish populations used in the F2 cross for the current study
          Key$Cross[Key$Publication == colnames(Warren2021)[study]],
          # The start and end positions of the marker on the scaffold
          Warren2021$Scaffold.Alignmnet.Start[r],
          Warren2021$Scaffold.Alignmnet.End[r],
          # Filler NAs for LOD and PVE scores
          NA,
          # Entire quantitative trait
          Warren2021[r,study],
          NA,
          # Study-specific information for the current publication
          Key$Study_Specific_Information[Key$Publication == colnames(Warren2021)[study]],
          # The name of the current publication
          colnames(Warren2021)[study]))
        # If current marker has multiple quantitative traits...
        }else{
          # split the quantitative traits into elements of a vector
          if(grepl(comma, Warren2021[r,study])){
            QTs <- str_split(Warren2021[r,study], comma)[[1]]
          }else if(grepl(period, Warren2021[r,study])){
            QTs <- str_split(Warren2021[r,study], period)[[1]]
          }
          # add a new row for each trait
          new_rows <- data.frame(
            # The current marker
            rep(Warren2021$Marker.ID[r], length(QTs)),
            # The scaffold for the current marker
            rep(Warren2021$Scaffold.ID[r], length(QTs)),
            # The linkage group for the current marker
            rep(Warren2021$LG[r], length(QTs)), 
            # The fish populations used in the F2 cross for the current study
            rep(Key$Cross[Key$Publication == colnames(Warren2021)[study]], length(QTs)),
            # The start and end positions of the marker on the scaffold
            rep(Warren2021$Scaffold.Alignmnet.Start[r], length(QTs)),
            rep(Warren2021$Scaffold.Alignmnet.End[r], length(QTs)),
            # Filler NAs for the informartion to be added later
            rep(NA, length(QTs)),
            QTs,
            rep(NA, length(QTs)),
            # Study-specific information for the current publication
            rep(Key$Study_Specific_Information[Key$Publication == colnames(Warren2021)[study]], length(QTs)),
            # The name of the current publication
            rep(colnames(Warren2021)[study], length(QTs))
          )
          colnames(new_rows) <- colnames(QTL)
          QTL <- rbind(QTL, new_rows)
      }
    }
  }
}

# Find all possible quantitative traits after removing the NA rows
QTL <- QTL[!is.na(QTL),]
QTL <- QTL[-(1327:nrow(QTL)),]
AllQTs <- levels(as.factor(QTL$Quantitative_Trait))
# Create a key for all quantitative traits
QTKey <- data.frame(
  Current = AllQTs,
  Replacement = c(NA,"Anal Fin Rays","Relative Body Condition","Dark Preference",
                  "Dentary Bone Length", "Eye Size", "Feeding Angle", 
                  "Eye Lens Size (Standardized against eye size)",
                  rep("Eye Lens Size (Standardized against standard length)",2),
                  "Maxillary Bone Length", "Maxillary Tooth Count",
                  "Melanophore Count (Anal fin)","Melanophore",
                  "Melanophore Count (Below dorsal fin)",
                  "Melanophore Count (Above and forward of eye)",
                  "Melanophore Count (Lateral flank)", "Pupil Diameter",
                  rep("Ribs", 2), "Superficial Neuromasts in Eye Orbit",
                  "Width of Sub-Orbital Bones", rep("Taste Buds", 2),
                  "Vibration Attraction Behavior", "Activity","Anal Fin Rays",
                  rep("Chemical Sense", 4), "Relative Body Condition",
                  "Dentary Bone Length", rep("Eye Size", 3), "Feeding Angle",
                  rep("Fin Placement", 3),"Retina (Ganglion Cell Layer)",
                  "Head Depth", "Retina (Inner Nuclear Layer)", "Jaw Angle",
                  "Standard Length",
                  "Eye Lens Size (Standardized against eye size)",
                  "Eye Lens Size (Standardized against standard length)",
                  "Maxillary Bone Length", "Maxillary Tooth Count",
                  "Melanophore Count (Anal fin)",rep("Melanophore", 4),
                  "Melanophore Count (Below dorsal fin)",
                  "Melanophore Count (Above and forward of eye)",
                  "Melanophore Count (Lateral flank)",
                  "Retina (Outer Nuclear Layer)", 
                  "Retina (Outer Plexiform Layer)", rep("Peduncle Depth", 2),
                  "Peduncle Length", "Ribs", "Schooling", "Chemical Sense",
                  "Superficial Neuromasts in Eye Orbit",
                  rep("Width of Sub-Orbital Bones", 2), rep("Taste Buds", 2),
                  rep("Tooth Count", 5), rep("Vibration Attraction Behavior", 2),
                  rep("Weight Loss", 2))
)

for(r in 1:nrow(QTL)){
  # For each quantitative trait, replace it with the replacement in the key
  if(length(QTL$Quantitative_Trait[r]) == length(QTKey$Replacement[QTKey$Current == QTL$Quantitative_Trait[r]])){
    QTL$Quantitative_Trait[r] <- QTKey$Replacement[QTKey$Current == QTL$Quantitative_Trait[r]]
  }else{
    print(r)
    break
  }
  # Remove the "chromosome" prefix from all chromosome numbers
  QTL$Chromosome[r] <- substr(QTL$Chromosome[r], 
                            start = nchar("Chromosome ") + 1, 
                            stop = nchar(QTL$Chromosome[r]))
}

# Write the information to QTL.csv
write.csv(QTL, "QTL_remapped_noLOD.csv", row.names = F)

# Manually add the...
# 1. Cross (for Kowalko PNAS entries)
# 2. LOD
# 3. PVE
# ... for each study