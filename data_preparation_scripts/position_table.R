# Annabel Perry
# June 18 2022
# The purpose of this script is to organize an Ensembl genome into the format
# needed to serve as a position table for CaveCrawler (or an app based on 
# CaveCrawler code)

setwd("~/Summer 2022/CaveCrawler Maintenance/Population Genetics")

# Open GTF
gtf <- read.table("Astyanax_mexicanus.Astyanax_mexicanus-2.0.106.gtf/Astyanax_mexicanus.Astyanax_mexicanus-2.0.106.gtf", fill = T)

# Slice file down to JUST scaffold, gene names, gene IDs, and start + end loci
# Note: You may have to change which columns to keep based on the format of your
# GTF
gtf <- gtf[,c(1,4:5,10,16)]
colnames(gtf) <- c("Scaffold","Start_Locus","End_Locus","Gene_ID","Gene_Name")

# Remove duplicate gene names as well as transcript IDs (the data in CaveCrawler 
# currently does not have the precision to tell the user which specific 
# transcript is yielding a given signal)
gtf <- gtf[(!grepl("ENSAMXT", gtf$Gene_Name) & !duplicated(gtf$Gene_Name)), ]

# Save file as a position table to be drawn upon in the "functions.R" file
write.csv(gtf, "PositionTable.csv", row.names = F)
