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

if (interactive()) {

    # User Interface
    ui <- fluidPage(
        theme = "sandstone.css",
        tags$head(tags$link(rel="shortcut icon", href="favicon.ico")),
        tags$style(type="text/css",
                   ".shiny-output-error { visibility: hidden; }",
                   ".shiny-output-error:before { visibility: hidden; }"),
        radioButtons("type",
                     label = "Search for Top/Bottom Number of Genes or Genes Above/Below a Statistical Value?",
                     choices = c("Number of Genes" = "Gene Count", "Statistic Value")
        ),
        radioButtons("statist",
                     label = "Statistic of Interest",
                     choices = c("Fst","Dxy","Tajima's D" = "TajimasD","Pi")),
        checkboxGroupInput("pops",
                     label = "Population(s) of Interest",
                     choices = c("Molino", "Pachon", "Rascon", "Rio Choy", "Tinaja")),

        # Only show this panel if the user wants to find the number of genes
        # with the greatest or smallest values for the stat of interest
        conditionalPanel(
            condition = "input.type == 'Gene Count'",
            sliderInput("gene_count", "How many genes would you like to see?",
                        min = 1, max = 1000, value = 10),
            radioButtons("TB",
                         label = "Output genes with largest or smallest values of the desired statistic?",
                         choices = c("Largest" = "top", "Smallest" ="bottom")),
            actionButton("GCDistTable_enter","Find Genes"),
            tableOutput("GCdist_tab"),
            textOutput("GCdist_wrnings"),
        ),
        # Only show this panel if the user wants to find genes with stat values
        # above or below a specific threshhold
        conditionalPanel(
            condition = "input.type == 'Statistic Value'",
            sliderInput("thrsh", "Threshhold statistical value: ",min = -3,
                        max = 3, value = 0, step = 0.05),
            radioButtons("TB",
                         label = "Output genes whose value is above or below the threshhold?",
                         choices = c("Above" = "top", "Below" ="bottom")),
            actionButton("SVDistTable_enter","Find Genes"),
            actionButton("SVDistPlot_enter", "Visualize"),
            plotOutput("SVdist_plot"),
            tableOutput("SVdist_tab"),
            textOutput("SVdist_wrnings"),
        )
    )


    # Server Logic
    server <- function(input, output) {
        s_table <- read.csv("data/AMexicanus_Genes_and_Stats.csv")
        s_table <- s_table[,(names(s_table) != "X")]
        dat <- read.table("data/Astyanax_mexicanus.Astyanax_mexicanus-2.0.104.gtf", fill = TRUE, skip = 5)
        position_table <- dat[dat$V3 == "gene",c(1,4,5,10,16)]
        GeneToGO <- read.csv("data/AMexGOTerms.csv", fill = T)
        GeneToGO <- GeneToGO[GeneToGO$Gene.names != "",]
        GeneToGO$Gene.names <- tolower(GeneToGO$Gene.names)


        StatDistTable <- function(in_type, UL, stat, thresh, stat_table, pops){
            library(tibble)
            # Create a two vectors of the indices corresponding to each population or pop
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
            # If no population pairs were found, output and error
            if(length(indices) == 0){
                null.df <- data.frame(matrix(nrow = 1, ncol = 6))
                names(null.df) <- c("Rank","Population(s)",stat,"Scaffold",
                                     "Gene Name","GO Term(s)")
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
                if((genes[g] %in% position_table$V16) & (genes[g] %in% GeneToGO$Gene.names)){
                    scaffs[g] = position_table$V1[position_table$V16 == genes[g]]
                    DF_GOs[g] = GeneToGO$Gene.ontology.IDs[GeneToGO$Gene.names == genes[g]]
                    # If gene is present in position table but NOT GO table, output NA for GO
                    # term but output real scaffold
                }else if((genes[g] %in% position_table$V16) & !(genes[g] %in% GeneToGO$Gene.names)){
                    scaffs[g] = position_table$V1[position_table$V16 == genes[g]]
                    DF_GOs[g] = "Not applicable"
                    # If gene is present in GO table but NOT position table, output NA for scaffold
                    # term but output real GO
                }else if(!(genes[g] %in% position_table$V16) & (genes[g] %in% GeneToGO$Gene.names)){
                    scaffs[g] = "Not applicable"
                    DF_GOs[g] = GeneToGO$Gene.ontology.IDs[GeneToGO$Gene.names == genes[g]]
                    # If gene is present in neither GO nor position tables, output NA for scaff
                    # and GO
                }else if(!(genes[g] %in% position_table$V16) & !(genes[g] %in% GeneToGO$Gene.names)){
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

            # Output df and warnings
            if(is.na(final_df[1,2])){
                final_df <- final_df[-1,]
            }
            return(list(wrnings, final_df))
        }
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

        SVDT <- eventReactive(input$SVDistTable_enter, valueExpr = {
            StatDistTable(input$type,
                              input$TB,
                              input$statist,
                              thresh = input$thrsh,
                              stat_table = s_table,
                              input$pops)
        }
        )
        output$SVdist_tab <- renderTable(SVDT()[[2]])
        output$SVdist_wrnings <- renderText(SVDT()[[1]])

        SVDP <- eventReactive(input$SVDistPlot_enter, valueExpr = {
            StatDistPlot(stat = input$statist,
                         UL = input$TB,
                         thresh = input$thrsh,
                         stat_table = s_table,
                         pops = input$pops)
        }
        )
        output$SVdist_plot <- renderPlot(SVDP())

        GCDT <- eventReactive(input$GCDistTable_enter, valueExpr = {
            StatDistTable(input$type,
                          input$TB,
                          input$statist,
                          thresh = input$gene_count,
                          stat_table = s_table,
                          input$pops)
        }
        )
        output$GCdist_tab <- renderTable(GCDT()[[2]])
        output$GCdist_wrnings <- renderText(GCDT()[[1]])
    }
}

# Run the application
shinyApp(ui = ui, server = server)
