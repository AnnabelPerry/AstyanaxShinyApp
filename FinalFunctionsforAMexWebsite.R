# Function which takes GOTerms, populations of interest, and statistics of interest
# and outputs a table
StatByChrTable <- function(GOTerm, GeneToGo, GoIDToNames, UpperLower, stat_vec, position_table, stat_table, all_pops){
  # Initialize vector for all GO terms of interest
  GOs <- c()
  # Check if user inputted a word/phrase or a GO ID. If the user inputted a 
  # word/phrase, find the name(s) which contains that word/phrase and set the GO
  # vector equal to the corresponding GO IDs.
  if(!(GOTerm %in% GoIDToNames$GO.ID)){
    GOs <- GoIDToNames$GO.ID[which(grepl(GOTerm, GoIDToNames$GO.Term, 
                                         ignore.case = T))]
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
    # If the GO term appears in the data frame of names AND the corresponding gene
    # occurs in the statistics vector, add the GO term and gene name
    if(GOs[g] %in% GeneToGO$Gene.ontology.IDs){
      gene_vec <- append(gene_vec, GeneToGO$Gene.names[grepl(GOs[g], GeneToGO$Gene.ontology.IDs)])
      found_GOs <- append(found_GOs, rep(GOs[g], length(gene_vec)))
      # If the GO term does NOT appear in the dataframe of names, skip it
    }else{
      next
    }
  }
  geneGOs <- data.frame(
    Gene = gene_vec,
    GO_ID = found_GOs
  )
  names(geneGOs) <- c("Gene", "GO_ID")
  # Extract all possible combinations of populations from populations of interest
  if(is.vector(all_pops)){
    pops <- combn(all_pops,2) 
  }
  
  # Create vectors in which to store values for later dataframe
  Statistic_Type_prelim <- c()
  Population_prelim <- c()
  Statistic_Value_prelim <- numeric()
  
  for(s in 1:length(stat_vec)){
    # Check if statistic of interest makes comparisons between TWO populations
    if((stat_vec[s] == "Fst") | (stat_vec[s] == "Dxy")){
      # Check if pops is a matrix
      if(is.matrix(pops)){
        # If so, iterate through each combination and output stat value for that
        # combination
        for(pair in 1:ncol(pops)){
          
          # If pops is a matrix, read the strings, find the column corresponding 
          # to the stat of interest for the populations in the vector, and set 
          # "index" equal to the column housing this statistic
          val <- which(grepl(pops[1, pair], names(stat_table))
                       & grepl(pops[2, pair], names(stat_table))
                       & grepl(stat_vec[s], names(stat_table)))
          # If statistic for populations-of-interest is not present, return error
          if(length(val) == 0){
            print(paste(c("Note: Statistic",stat_vec[s],
                          "is not present for the populations",pops[1, pair],"and",
                          pops[2, pair]),collapse = " "))
          }else{
            # Create a row to add to the indices dataframe
            temp_str <- paste(c(pops[1, pair],"-",pops[2, pair]), collapse = "")
            Statistic_Type_prelim <- append(Statistic_Type_prelim,stat_vec[s])
            Population_prelim <- append(Population_prelim,temp_str)
            Statistic_Value_prelim <- append(Statistic_Value_prelim,val)
          }
          # If pops is NOT a matrix, return an error
        }
      }else if(!is.matrix(pops)){
        print(paste(c("Warning: Only one population, ", pops, 
                      " supplied for the two-population statistic ", stat_vec[s], 
                      ". 'pops' likely has only one entry."), collapse = ""))
        # Check if statistic of interest makes comparisons between ONE population
      }
    }else if((stat_vec[s] == "TajimasD") | (stat_vec[s] == "Pi")){
      # Iterate through each individual population
      for(p in 1:length(all_pops)){
        # If pops is a string, set indx equal to the column housing the stat of
        # interest for this population
        val <- which(grepl(all_pops[p], names(stat_table))
                     & grepl(stat_vec[s], names(stat_table)))
        # If statistic for populations-of-interest is not present, return a 
        # warning 
        if(length(val) == 0){
          print(paste(c("Note: Statistic",stat_vec[s],
                        "is not present for the population",all_pops[p]),
                      collapse = " "))
        }else{
          # Create a row to add to the indices dataframe
          Statistic_Type_prelim <- append(Statistic_Type_prelim,stat_vec[s])
          Population_prelim <- append(Population_prelim,all_pops[p])
          Statistic_Value_prelim <- append(Statistic_Value_prelim,val)
        }
      }
    }else{
      return(paste(c("ERROR: Invalid statistic name: ", stat_vec[s]), 
                   collapse = ""))
    }
  }
  
  # For each GO Term in the vector of GO Terms-of-interest, find...
  # 1. All gene names which occur in the positions table AND in the statistics
  # table
  # 2. The scaffolds of those genes
  # 3. The starting positions of those genes
  # 4. The ending positions of those genes
  # 5. The GO terms associated with those genes
  # 6. Each of the 4 statistic types
  # 7. Each of the populations/population combinations
  # 8. Each of the statistical values
  
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
  
  for(s in 1:length(Statistic_Type_prelim)){
    for(g in 1:length(geneGOs$Gene)){
      # Check if current gene is in stat AND position table
      # If so...
      if((geneGOs$Gene[g] %in% stat_table$Gene_Name) & 
         (geneGOs$Gene[g] %in% position_table$V16)){
        # Output the current gene
        Gene <- append(Gene, geneGOs$Gene[g])
        # Output scaffold of current gene
        Scaffold <- append(Scaffold, 
                           position_table$V1[position_table$V16 == geneGOs$Gene[g]])
        # Output starting position of the current gene
        Start_Position <- append(Start_Position, 
                                 position_table$V4[position_table$V16 == geneGOs$Gene[g]])
        # Output the ending position of the current gene
        End_Position <- append(End_Position, 
                               position_table$V5[position_table$V16 == geneGOs$Gene[g]])
        # Output ALL GO terms associated with the current gene
        GO_IDs <- append(GO_IDs,
                         paste(geneGOs$GO_ID[geneGOs$Gene == geneGOs$Gene[g]], collapse = "; "))
        # Output the current statistic type
        Statistic_Type <- append(Statistic_Type, Statistic_Type_prelim[s])
        # Output the population(s)
        Population <- append(Population, Population_prelim[s])
        # Output the statistic value
        Statistic_Value <- append(Statistic_Value, 
                                  stat_table[stat_table$Gene_Name == geneGOs$Gene[g], 
                                             Statistic_Value_prelim[s]])
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
                          Statistic_Value
  )
  return(output_df)
}

# Takes table from StatByChrTable as input and outputs a unique plot for each
# scaffold-statistic pair
StatByChrGraph <- function(Full_Table){
  library(ggplot2)
  library(gridExtra)
  
  # Count the number of unique scaffolds in the table
  unique_scaffs <- levels(as.factor(Full_Table$Scaffold))
  num_scaff <- length(unique_scaffs)
  
  # Find number of unique statistics so you can create a distinct plot for each
  # statistic-scaffold pair
  unique_stats <- levels(as.factor(Full_Table$Statistic_Type))
  num_stats <- length(unique_stats)
  
  # Initialize a list in which to store all plots
  plist <- vector(mode = "list", length = num_scaff*num_stats)
  
  # Start at the 0th entry in the list
  p = 0
  # For each scaffold-statistic combo, create a new plot
  for(s in 1:num_stats){
    # For each scaffold represented by the selected genes, create a new plot
    for(c in 1:num_scaff){
      # Each time you reach a new scaffold-statistic combination, move to a new 
      # position in the plot list 
      p = p + 1
      # Find all rows in the table whose basic statistic and scaffold match the
      # current statistic and scaffold
      plot_dat <- data.frame(
        Locus = Full_Table$Start_Position[(Full_Table$Scaffold == unique_scaffs[c])
                                          & (Full_Table$Statistic_Type == unique_stats[s])],
        Value = Full_Table$Statistic_Value[(Full_Table$Scaffold == unique_scaffs[c])
                                           & (Full_Table$Statistic_Type == unique_stats[s])],
        Populations = Full_Table$Population[(Full_Table$Scaffold == unique_scaffs[c])
                                            & (Full_Table$Statistic_Type == unique_stats[s])]
      )
      
      # Plot the statistic values for all genes which occur on the scaffold
      # of interest, colored by population or population pair
      plist[[p]] <- ggplot(plot_dat, aes(x=Locus, y=Value, color=Populations)) + 
        ylab(stat_vec[s]) +
        xlab("") +
        ggtitle(paste(c("Scaffold:", unique_scaffs[c]), collapse = " ")) +
        geom_point() +
        theme_bw()
      #Store the plot in a vector
    }
  }
  return(plist)
}

# Creates a table of the top or bottom quantile of a stat-of-interest OR finds 
# all values above/below the stat-of-interest
StatDistTable <- function(in_type, UL, stat, thresh, stat_table, pops){
  library(tibble)
  # Create a two vectors of the indices corresponding to each population or pop
  # pair's statistic values
  indices <- c()
  pop_strings <- c()
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
        print(paste(c("Note: Statistic",stat,
                      "is not present for the populations",two_pops[1, pair],"and",
                      two_pops[2, pair]),collapse = " "))
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
        print(paste(c("Note: Statistic",stat,
                      "is not present for the population",pops[p]),
                    collapse = " "))
      }else{
        pop_strings <- append(pop_strings,pops[p])
        indices <- append(indices,val)
      }
    }
  }else{
    return("ERROR: Invalid statistic name")
  }
  
  # Initialize vectors of gene names, populations, and statistic values
  genes <- c()
  DF_pops <- c()
  stat_vals <- c()
  
  
  # Check if statistic value or gene count was entered
  # If value was entered...
  if(in_type == "Statistic Value"){
    # Check whether top or bottom proportion was requested
    # If higher proportion was requested, iterate through each index
    if(UL == "top"){
      for(i in 1:length(indices)){
        # First, remove all NA values for this index
        stat_table <- stat_table[!is.na(stat_table[,indices[i]]),]
        # For each index, collect all genes, scaffolds, populations, and values
        # whose stat values fall above the entered value
        genes <- append(genes,stat_table$Gene_Name[stat_table[,indices[i]] >= thresh])
        DF_pops <- append(DF_pops,rep(pop_strings[i], 
                                      length(stat_table$Gene_Name[stat_table[,indices[i]] >= thresh])))
        stat_vals <- append(stat_vals,stat_table[stat_table[,indices[i]] >= thresh,indices[i]])
      }
      # If lower tail was requested, iterate through each index
    }else if(UL == "bottom"){
      for(i in 1:length(indices)){
        # First, remove all NA values for this index
        stat_table <- stat_table[!is.na(stat_table[,indices[i]]),]
        # For each index, collect all genes, populations, and statistic values
        # whose values fall in lowest tail
        genes <- append(genes,stat_table$Gene_Name[stat_table[,indices[i]] <= thresh])
        DF_pops <- append(DF_pops,rep(pop_strings[i], 
                                      length(stat_table$Gene_Name[stat_table[,indices[i]] <= thresh])))
        stat_vals <- append(stat_vals,stat_table[stat_table[,indices[i]] <= thresh,indices[i]])
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
          top_genes <- order(stat_table[,indices[i]], decreasing = T)[1:thresh]
          # Collect the genes with the highest values, as well as the associated 
          # populations and values
          genes <- append(genes,stat_table$Gene_Name[top_genes])
          DF_pops <- append(DF_pops,rep(pop_strings[i],length(top_genes)))
          stat_vals <- append(stat_vals,stat_table[top_genes,indices[i]])
        }
        # If statistic is a one-population statistic, output collect the N genes 
        # with the HIGHEST stat value, regardless of pop
      }else if(stat_type == "One Pop"){
        # Collect stat values for ALL indices into a 3 vectors: row
        # in one column, population name in another,and stat value in the other
        all_stats <- c()
        all_pops <- c()
        all_genes <- c()
        for(i in 1:length(indices)){
          # First, remove all NA values for this index
          stat_table <- stat_table[!is.na(stat_table[,indices[i]]),]
          all_stats <- append(all_stats, stat_table[,indices[i]])
          all_pops <- append(all_pops,rep(pop_strings[i],
                                          length(stat_table[,indices[i]])))
          all_genes <- append(all_genes, stat_table$Gene_Name)
        }
        # Organize vectors into a dataframe
        temp_df <- data.frame(
          all_stats,
          all_pops,
          all_genes
        )
        # Retreat the parallel indeices for the N highest genes
        par_indices <- order(temp_df$all_stats, decreasing = T)[1:thresh]
        # Retrieve the gene names, population names
        genes <- append(genes,temp_df$all_genes[par_indices])
        DF_pops <- append(DF_pops,temp_df$all_pops[par_indices])
        stat_vals <- append(stat_vals,temp_df$all_stats[par_indices])
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
          # Collect the genes with the highest values, as well as the associated 
          # populations and values
          genes <- append(genes,stat_table$Gene_Name[bottom_genes])
          DF_pops <- append(DF_pops,rep(pop_strings[i],length(bottom_genes)))
          stat_vals <- append(stat_vals,stat_table[bottom_genes,indices[i]])
        }
        # If statistic is a one-population statistic, output collect the N genes 
        # with the HIGHEST stat value, regardless of pop
      }else if(stat_type == "One Pop"){
        # Collect stat values for ALL indices into a 3 vectors: row
        # in one column, population name in another,and stat value in the other
        all_stats <- c()
        all_pops <- c()
        all_genes <- c()
        for(i in 1:length(indices)){
          # First, remove all NA values for this index
          stat_table <- stat_table[!is.na(stat_table[,indices[i]]),]
          all_stats <- append(all_stats, stat_table[,indices[i]])
          all_pops <- append(all_pops,rep(pop_strings[i],
                                          length(stat_table[,indices[i]])))
          all_genes <- append(all_genes, stat_table$Gene_Name)
        }
        # Organize vectors into a dataframe
        temp_df <- data.frame(
          all_stats,
          all_pops,
          all_genes
        )
        # Retreat the parallel indeices for the N highest genes
        par_indices <- order(temp_df$all_stats, decreasing = F)[1:thresh]
        # Retrieve the gene names, population names
        genes <- append(genes,temp_df$all_genes[par_indices])
        DF_pops <- append(DF_pops,temp_df$all_pops[par_indices])
        stat_vals <- append(stat_vals,temp_df$all_stats[par_indices])
      }
    }
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
    if((genes[g] %in% position_table$V16) & (genes[g] %in% GeneToGO$Gene.names)){
      scaffs[g] = position_table$V1[position_table$V16 == genes[g]]
      DF_GOs[g] = GeneToGO$Gene.ontology.IDs[GeneToGO$Gene.names == genes[g]]
      # If gene is present in position table but NOT GO table, output NA for GO
      # term but output real scaffold
    }else if((genes[g] %in% position_table$V16) & !(genes[g] %in% GeneToGO$Gene.names)){
      scaffs[g] = position_table$V1[position_table$V16 == genes[g]]
      DF_GOs[g] = NA
      # If gene is present in GO table but NOT position table, output NA for scaffold
      # term but output real GO
    }else if(!(genes[g] %in% position_table$V16) & (genes[g] %in% GeneToGO$Gene.names)){
      scaffs[g] = NA
      DF_GOs[g] = GeneToGO$Gene.ontology.IDs[GeneToGO$Gene.names == genes[g]]
      # If gene is present in neither GO nor position tables, output NA for scaff
      # and GO
    }else if(!(genes[g] %in% position_table$V16) & !(genes[g] %in% GeneToGO$Gene.names)){
      scaffs[g] = NA
      DF_GOs[g] = NA
    }
  }
  # Create a dataframe of scaffolds, gene names, GO terms, statistic types, and 
  # stat values
  prelim_df <- data.frame(
    DF_pops,
    stat_vals,
    scaffs,
    genes,
    DF_GOs
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
        new_rows <- rows[order(stat_vals, decreasing = T),]
        new_rows <- rows[1:nrow(rows),]
        # Collect number of genes found above threshhold for THIS population
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
        new_rows <- rows[order(stat_vals, decreasing = F),]
        new_rows <- rows[1:nrow(rows),]
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
    "GO Term(s)"
  )
  
  # Output df
  if(is.na(final_df[1,2])){
    final_df <- final_df[-1,]
  }
  return(final_df)
}

# Plots statistic values for all genes across genome
StatDistPlot <- function(stat, UL, thresh, stat_table, pops){
  library("WVPlots")
  # EDIT: Check if threshold is within appropriate range for statistic-of-interest. If
  # not, return an error.
  
  # Check if statistic of interest makes comparisons between TWO populations
  if((stat == "Fst") | (stat == "Dxy")){
    # Check if pops is a vector
    if(is.vector(pops)){
      # If pops is a vector, read the strings, find the column corresponding 
      # to the stat of interest for the populations in the vector, and set 
      # "indx" equal to the column housing this statistic
      indx <- which(grepl(pops[1], names(stat_table))
                    & grepl(pops[2], names(stat_table))
                    & grepl(stat, names(stat_table)))
      # If statistic for populations-of-interest is not present, return error
      if(length(indx) == 0){
        return(paste(c("ERROR: Statistic",stat,
                       "is not present for the populations",pops[1],"and",
                       pops[2]),collapse = " "))
      }
      # If pops is NOT a vector, return an error
    }else if(!is.vector(pops)){
      return("ERROR: Only one population supplied for a two-population statistic
             \n Either 'pops' is not a vector or 'pops' has only one entry.\n")
      # Check if statistic of interest makes comparisons between ONE population
    }
  }else if((stat == "TajimasD") | (stat == "Pi")){
    # Ensure that only one population was entered
    if(is.character(pops) & length(pops) == 1){
      # If pops is a string, set indx equal to the column housing the stat of
      # interest for this population
      indx <- which(grepl(pops, names(stat_table))
                    & grepl(stat, names(stat_table)))
      # If statistic for populations-of-interest is not present, return error
      if(length(indx) == 0){
        return(paste(c("ERROR: Statistic",stat,
                       "is not present for the population",pops),
                     collapse = " "))
      }
      # If more than one population was entered, return an error
    }else if(!is.character(pops) | length(pops) != 1){
      return("ERROR: Detected inappropriate number of populations for a 
              one-population statistic\n Either 'pops' is not a string or 'pops' 
             was not supplied.\n")
    }
  }else{
    return("ERROR: Invald statistic name")
  }
  # Remove all rows where stat-of-interest is NA
  filt_table <- stat_table[!is.na(stat_table[,indx]),]
  # Create dataframe with gene names, appropriate stat values, and GO terms
  # EDIT: Include GO terms
  df <- data.frame(
    Genes = filt_table$Gene_Name,
    Stats = filt_table[,indx]
  )
  # Rename df satistics column with specific statistic name
  names(df)[2] <- stat
  # Pick a tail based on whether "top" or "bottom" was specified
  if(UL == "top"){
    t = "right"
  }else if(UL == "bottom"){
    t = "left"
  }
  # Create density plot. Title must be different depending on whether you have 1 
  # or two populations
  if(is.vector(pops)){
    ShadedDensity(frame = df, 
                  xvar = stat,
                  threshold = thresh,
                  title = paste(c(stat, " values for ", pops[1], " and ", pops[2]), collapse = ""),
                  tail = t)
  }else{
    ShadedDensity(frame = df, 
                  xvar = stat,
                  threshold = thresh,
                  title = paste(c(stat, " values for ", pops), collapse = ""),
                  tail = t)
  }
}

# Function: Find percent of genes with most extreme value for statistic of interest
TranscTable <- function(morph1, morph2, condition, direction, percent, in_table, GOTable = GeneToGO){
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
      ROIs$ï..Gene_stable_ID[1:n.rows],
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

# Input is a string of gene names, string of gene stable IDs, or phrase 
# associated with a gene
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


# Test examples
setwd("~/Fall 2021/Capstone")
dat <- read.table("Astyanax_mexicanus.Astyanax_mexicanus-2.0.104.gtf", fill = TRUE, skip = 5)
position_table <- dat[dat$V3 == "gene",c(1,4,5,10,16)]

condition_control <- read.csv("ShinyInputData/Morph-Control_TranscData.csv")
morph1.morph2 <- read.csv("ShinyInputData/Toy_RioChoyPachon.csv")

GeneToGO <- read.csv("ShinyInputData/AMexGOTerms.csv", fill = T)
GeneToGO <- GeneToGO[GeneToGO$Gene.names != "",]
# Convert gene names in GeneToGO to lowercase to make compatible with other tables
GeneToGO$Gene.names <- tolower(GeneToGO$Gene.names)

GoIDToNames <- read.table("ShinyInputData/GOIDs_and_Names.txt", fill = T, sep = "\t", header = T)
UpperLower <- read.table("ShinyInputData/GOTermAssociations.txt", fill = T, sep = "\t", header = T)

stat_table <- read.csv("ShinyInputData/AMexicanus_Genes_and_Stats.csv")
stat_table <- stat_table[,(names(stat_table) != "X")]

# Find bottom 30 values for Tajima's D across all populations
test1 <- StatDistTable(in_type = "Gene Count", UL = "bottom", stat = "TajimasD",
                   thresh = 30,stat_table,
                   pops = c("Tinaja","Rascon","Pachon","Molino","RioChoy"))
# Find all genes with Fst values above .95
test2 <- StatDistTable(in_type = "Statistic Value", UL = "top", stat = "Fst",
                   thresh = 0.95,stat_table,
                   pops = c("Tinaja","Rascon","Pachon","Molino","RioChoy"))


test3 <- StatByChrTable(GOTerm = "GO:0035639", GeneToGO, GoIDToNames, UpperLower,
                        stat_vec <- c("Fst","Pi"), position_table, stat_table, 
                        all_pops = c("Tinaja","Molino","Pachon","RioChoy","Rascon"))

################################# 10-07 New Functions ###########################

############################## Transcription Function ##########################
# NEW:Find 30% of genes with greatest downregulation in Pachon relative to Rio Choy
test4 <- TranscTable("Pachon","Rio Choy", "Between morph", "Downregulated", 30)
# NEW:Find 20% of genes with greatest upregulation in response to sleep deprivation
# in Rio Choy
test5 <- TranscTable("Rio Choy","Control", "Sleep deprivation", "Upregulated", 20)

############################ Gene-Centered Function ############################
# NEW: Find transcription and pop-gen data on a comma-separated group of genes
test6 <- GeneCentered("gna14, cyp26a1", stat_table, GeneToGO, condition_control,
                      position_table)
# NEW: Find transcription and pop-gen data on a comma-separated group of IDs
test7 <- GeneCentered("ENSAMXG00000014853, ENSAMXG00000006096", stat_table, GeneToGO, condition_control,
                      position_table)
# NEW: Find transcription and pop-gen data on genes associated with a phrase
test8 <- GeneCentered("melatonin receptor", stat_table, GeneToGO, condition_control,
                      position_table)
