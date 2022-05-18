#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

############################### Source Functions ###############################
source("functions/CaveCrawler_functions.R")

  ui = fluidPage(
    setBackgroundColor("white"),
    chooseSliderSkin("Flat", color = "#e8c4c2"),
    theme = "style.css",
    tags$head(tags$link(rel="shortcut icon", href="favicon.ico")),
    tags$style(type="text/css",
               ".shiny-output-error { visibility: hidden; }",
               ".shiny-output-error:before { visibility: hidden; }",
               'body {color:black;}'),
    # Change background color of tabs
    tags$style(HTML("
    .tabbable > .nav > li > a                  {background-color: white; border-color: grey;}
    .tabbable > .nav > li[class=active]    > a {background-color: #e8c4c2; border-color: #e4867e}
    .tabbable > .nav > li > a:hover {background-color: #e8c4c2; border-color: #e4867e}
  ")),
    tabsetPanel(
      tabPanel(h2("Home"), fluid = TRUE,
               h1("Welcome to CaveCrawler"),
               img(src="Astyanax_Evolution_GIF.gif", align = "center",height='250px',width='500px'),
               br(),
               br(),
               "CaveCrawler is a reactive web interface for bioinformatic 
               analysis of data in the Mexican tetra (Astyanax mexicanus), an 
               emerging evolutionary model organism.",
               br(),
               br(),
               "CaveCrawler consists of 4 modules: a Gene Search page for 
               querying data about specific genes, a Transcription page for 
               finding genes whose transcriptional levels differ between samples,
               a Population Genetics page for investigating statistics on 
               diversity and selection, and a GO Term Info page for identifying 
               and obtaining information on GO terms-of-interest.",
               br(),
               br(),
               "Finally, there is the Data Sources module, which describes the 
               publications from which each dataset was obtained, as well as the 
               dates for Ensembl and UniProt-derived information.",
               br(),
               br(),
               "To request that new data be integrated into CaveCrawler, please 
               email Heath Blackmon at hblackmon@bio.tamu.edu or Alex Keene at 
               akeene@bio.tamu.edu",
               br(),
               br(),
               "To cite CaveCrawler, please cite our paper, currently available on BioRXiv: 'CaveCrawler: An interactive analysis suite for cavefish bioinformatics'",
               br(),
               plotOutput("home_plot")
      ),
      tabPanel(h2("Gene Search"), fluid = TRUE,
               sidebarLayout(
                 sidebarPanel(id = "sidebar",
                  checkboxGroupInput("GSbools",
                     label = "What data would you like to see?",
                     choices = c("Position", "Transcription","Population Genetics","GO")),
                   searchInput(
                     inputId = "Gene_search",
                     label = "Gene name, gene stable ID, GO ID, or phrase",
                     placeholder = "mtnr1al, ENSAMXG00000010894, GO:0016021, melatonin receptor, etc...",
                     btnSearch = icon("search"),
                     width = "450px"
                   ),
                   conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                    tags$div("Crawling through the data...",id="loadmessage"))
                 ),
                 mainPanel(id = "main",
                   tags$head(tags$style(".download{background-color:#c8feca;} .download{color: #71c596 !important;} .download{border-color: #71c596 !important;}")),
                   # If position data is requested, output position data
                   conditionalPanel(
                     condition = "input.GSbools.includes('Position')",
                     fluidRow(class = "text-center",
                              column(width = 4,
                                     h1("Position Data"),
                                     radioButtons("GSPwhich_sort",
                                                  label = "Sort by...",
                                                  choices = c(
                                                    "Scaffold" = "GSPsc",
                                                    "Start Locus" = "GSPsl",
                                                    "End Locus" = "GSPel"
                                                  )
                                     ),
                                     conditionalPanel(
                                       condition = "input.GSPwhich_sort == 'GSPsc'",
                                       radioButtons("GSPsc_dir",
                                                    label = "... from...",
                                                    choices = c(
                                                      "High to Low" = "GSPsc_hl",
                                                      "Low to High" = "GSPsc_lh"
                                                    )
                                       )
                                     ),
                                     conditionalPanel(
                                       condition = "input.GSPwhich_sort == 'GSPsl'",
                                       radioButtons("GSPsl_dir",
                                                    label = "",
                                                    choices = c(
                                                      "High to Low" = "GSPsl_hl",
                                                      "Low to High" = "GSPsl_lh"
                                                    )
                                       )
                                     ),
                                     conditionalPanel(
                                       condition = "input.GSPwhich_sort == 'GSPel'",
                                       radioButtons("GSPel_dir",
                                                    label = "",
                                                    choices = c(
                                                      "High to Low" = "GSPel_hl",
                                                      "Low to High" = "GSPel_lh"
                                                    )
                                       )
                                     ),
                                     downloadButton("GSPosDL", "Download Position Data", class = "download")   
                              ),
                              column(width = 6,
                                     tableOutput("GSPos_table")
                              )
                              )
                     ),
                   br(),
                   # If Transcription data is requested, output Transcription data
                   conditionalPanel(
                     condition = "input.GSbools.includes('Transcription')",
                     fluidRow(
                       column(width = 5,
                              h1("Transcription Data"),
                              radioButtons("GSTwhich_sort",
                                           label = "Sort by...",
                                           choices = c("logFC", "p")
                              ),
                              conditionalPanel(
                                condition = "input.GSTwhich_sort == 'logFC'",
                                radioButtons("GSTFC_dir",
                                             label = "... from...",
                                             choices = c(
                                               "High to Low" = "GSTFC_hl",
                                               "Low to High" = "GSTFC_lh"
                                             )
                                )
                              ),
                              conditionalPanel(
                                condition = "input.GSTwhich_sort == 'p'",
                                radioButtons("GSTp_dir",
                                             label = "... from...",
                                             choices = c(
                                               "High to Low" = "GSTp_hl",
                                               "Low to High" = "GSTp_lh"
                                             )
                                )
                              ),
                              downloadButton("GSTranscDL", "Download Transcription Data", class = "download")
                       ),
                       column(width = 6,
                              tableOutput("GSTransc_table")
                              )
                     )
                     ),
                   br(),
                   # If Population Genetics data is requested, output Population
                   # Genetics data
                   conditionalPanel(
                     condition = "input.GSbools.includes('Population Genetics')",
                     fluidRow(class = "text-center",
                              column(width = 4,
                                     h1("Popgen Data"),
                                     radioButtons("GSPG_dir",
                                                  label = "Sort statistic values from...",
                                                  choices = c(
                                                    "High to Low" = "GSPG_hl",
                                                    "Low to High" = "GSPG_lh"
                                                  )
                                     ),
                                     downloadButton("GSPopgenDL", "Download Popgen Data", class = "download")
                                     
                                     ),
                              column(width = 6,
                                     tableOutput("GSPopgen_table")
                              )
                          ),
                     ),
                   br(),
                   # If GO data is requested, output GO data
                   conditionalPanel(
                     condition = "input.GSbools.includes('GO')",
                     fluidRow(
                       column(width = 4,
                              h1("GO Data"),
                              downloadButton("GSGODL", "Download GO Data", class = "download")
                       ),
                       column(width = 6,
                         tableOutput("GSGO_table")
                       )
                     )
                   ),
                   # Regardless of what was inputted, output warnings
                   br(),
                   fluidRow(
                     h1("Warnings:"),
                     textOutput("GSwarnings")
                   )
                   
               )
          )
      ),
      tabPanel(h2("Transcription"), fluid = TRUE,
               sidebarLayout(
                 sidebarPanel(id = "sidebar",
                   selectInput(
                     inputId = "morph1",
                     label = "Compare...",
                     choices = c("Pachon","Molino","Tinaja","Rascon","Rio Choy")
                   ),
                   selectInput(
                     inputId = "morph2",
                     label = "to...",
                     choices = c("")
                   ),
                   radioButtons("trstat",
                                label = textOutput("trstat_label"),
                                choices = c("logFC", "p")
                   ),
                   radioButtons("direction",
                                label = textOutput("dir_label"),
                                choices = c("Above", "Below")
                   ),
                   
                   # If morph2 IS set to control, populate a drop-down of available conditions
                   conditionalPanel(
                     condition = "input.morph2 == 'Control'",
                     selectInput(
                       inputId = "condition",
                       label = "Condition?",
                       choices = c("Sleep deprivation")
                     )
                   ),
                   
                   conditionalPanel(
                     condition = "input.trstat == 'logFC'",
                     uiOutput("FCthresh_updater")
                   ),
                   conditionalPanel(
                     condition = "input.trstat == 'p'",
                     sliderInput(inputId = "pthresh",
                                 label = textOutput("pthresh_label"),
                                 min = 0, max = 1, value = 0.5, step = 0.000001)
                     
                   ),
                   actionButton("Transc_enter","Find Genes"),
                   # Change sorting of output table
                   radioButtons("TRwhich_sort",
                                label = "Sort by...",
                                choices = c("logFC", "p")
                   ),
                   conditionalPanel(
                     condition = "input.TRwhich_sort == 'p'",
                     radioButtons("TRsort_p",
                                  label = "...from...",
                                  choices = c("High to Low" = "TRp_hl", 
                                              "Low to High" = "TRp_lh")
                     )
                   ),
                   conditionalPanel(
                     condition = "input.TRwhich_sort == 'logFC'",
                     radioButtons("TRsort_FC",
                                  label = "...from...",
                                  choices = c("High to Low" = "TRFC_hl", 
                                              "Low to High" = "TRFC_lh")
                     )
                  ),
                  conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                    tags$div("Crawling through the data...",id="loadmessage"))
                  ),
                 mainPanel(fluidRow(
                   conditionalPanel(condition = "Transc_enter",
                                    downloadButton("TranscDL", "Download", class = "download"),
                                    tableOutput("transc_table_out"),
                                    textOutput("transc_wrnings_out")
                   )
                 )
                 )
               )
      ),
        tabPanel(h2("Population Genetics"), fluid = TRUE,
                 sidebarLayout(
                   sidebarPanel(id = "sidebar2",
                     radioButtons("which_function",
                                  label = "Would you like to find genes which are outliers with respect to a statistic-of-interest or find the statistic values for all genes related to a GO-term-of-interest?",
                                  choices = c("Outliers" = "distr_func",
                                              "GO-term-of-interest" = "stat_by_chr_func")
                     ),
                     conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                      tags$div("Crawling through the data...",id="loadmessage"))
                   ),
                   mainPanel(
                     br(),
                     br(),
                     br(),
                     br(),
                     br(),
                     br(),
                     br(),
                     br(),
                     # Display download button appropriate to specific conditions
                     conditionalPanel(
                       condition = "input.which_function == 'stat_by_chr_func'",
                       downloadButton("SBCDL", "Download", class = "download")
                     ),
                     conditionalPanel(
                       condition = "input.which_function == 'distr_func'",
                       conditionalPanel(
                         condition = "input.type == 'Gene Count'",
                         downloadButton("GCDistDL", "Download", class = "download")
                       ),
                       conditionalPanel(
                         condition = "input.type == 'Statistic Value'",
                         downloadButton("SVDistDL", "Download", class = "download")
                       )
                     )
                   )
                 ),
                 conditionalPanel(
                   condition = "input.which_function == 'distr_func'",
                   sidebarLayout(
                     sidebarPanel(id = "sidebar",
                       radioButtons("type",
                                    label = "Search for Top/Bottom Number of Genes or Genes Above/Below a Statistical Value?",
                                    choices = c("Number of Genes" = "Gene Count", "Statistic Value")
                       ),
                       radioButtons("dist_statist",
                                    label = "Statistic of Interest",
                                    choices = c("Fst","Dxy","Tajima's D" = "TajimasD","Pi")),
                       checkboxGroupInput("dist_pops",
                                          label = "Population(s) of Interest",
                                          choices = c("Molino", "Pachon",
                                                      "Rascon", "Rio Choy" = "RioChoy",
                                                      "Tinaja",
                                                      "Chica 1" = "Chica1",
                                                      "Chica 2" = "Chica2")),

                       # Only show this panel if the user wants to find the number of genes
                       # with the greatest or smallest values for the stat of interest
                       conditionalPanel(
                         condition = "input.type == 'Gene Count'",
                         sliderInput("gene_count", "How many genes would you like to see?",
                                     min = 1, max = 1000, value = 10),
                         radioButtons("gcTB",
                                      label = "Output genes with largest or smallest values of the desired statistic?",
                                      choices = c("Largest" = "top", "Smallest" ="bottom")),
                         actionButton("GCDistTable_enter","Find Genes"),
                         conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                          tags$div("Crawling through the data...",id="loadmessage"))
                       ),
                       # Only show this panel if the user wants to find genes with stat values
                       # above or below a specific threshhold
                       conditionalPanel(
                         condition = "input.type == 'Statistic Value'",
                         uiOutput("thrsh_slider"),
                         radioButtons("svTB",
                                      label = "Output genes whose value is above or below the threshhold?",
                                      choices = c("Above" = "top", "Below" ="bottom")),
                         actionButton("SVDistTable_enter","Find Genes"),
                         uiOutput("SVDistPlot_enter"),
                         conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                          tags$div("Crawling through the data...",id="loadmessage"))
                       )
                     ),
                     mainPanel(
                       conditionalPanel(
                         condition = "input.type == 'Gene Count'",
                         tableOutput("GCdist_tab"),
                         textOutput("GCdist_wrnings")
                       ),
                       # Only show this panel if the user wants to find genes with stat values
                       # above or below a specific threshhold
                       conditionalPanel(
                         condition = "input.type == 'Statistic Value'",
                         plotOutput("SVdist_plot"),
                         textOutput("SVdist_plot_wrnings"),
                         tableOutput("SVdist_tab"),
                         textOutput("SVdist_wrnings")
                       )
                     )
                   )
                 ),
                 conditionalPanel(
                   condition = "input.which_function == 'stat_by_chr_func'",
                   sidebarLayout(
                     sidebarPanel(id = "sidebar",
                       checkboxGroupInput("sbc_statist",
                                          label = "Statistic of Interest",
                                          choices = c("Fst","Dxy","Tajima's D" = "TajimasD","Pi")),
                       checkboxGroupInput("sbc_pops",
                                          label = "Population(s) of Interest",
                                          choices = c("Molino", "Pachon",
                                                      "Rascon", "Rio Choy",
                                                      "Tinaja",
                                                      "Chica 1" = "Chica1",
                                                      "Chica 2" = "Chica2")),
                       searchInput(
                         inputId = "GO_search", label = "GO ID or Term",
                         placeholder = "GO:0000001, mitochondrion, etc...",
                         btnSearch = icon("search"),
                         width = "450px"
                       ),
                       # If stat-by-chr is checked, enable sorting
                       radioButtons("SBCwhich_sort",
                                    label = "Sort table by...",
                                    choices = c("Start Position" = "SBCsp_sort",
                                                "End Position" = "SBCep_sort",
                                                "Statistic Value" = "SBCsv_sort")
                       ),
                       # Create ascending or descending sort based on initial input
                       conditionalPanel(
                         condition = "input.SBCwhich_sort == 'SBCsp_sort'",
                         radioButtons("SBCsp_dir", label = "Sort order:",
                                      choices = c("High to Low" = "sp_asc", 
                                                  "Low to High" ="sp_desc"))
                       ),
                       conditionalPanel(
                         condition = "input.SBCwhich_sort == 'SBCep_sort'",
                         radioButtons("SBCep_dir",
                                      label = "Sort order:",
                                      choices = c("High to Low" = "ep_asc", 
                                                  "Low to High" ="ep_desc"))
                       ),
                       conditionalPanel(
                         condition = "input.SBCwhich_sort == 'SBCsv_sort'",
                         radioButtons("SBCsv_dir",
                                      label = "Sort order:",
                                      choices = c("High to Low" = "sv_asc", 
                                                  "Low to High" ="sv_desc"))
                       ),
                       uiOutput("stat_PlotSelect"),
                       uiOutput("scaff_PlotSelect"),
                       uiOutput("SBCP_enter"),
                       conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                        tags$div("Crawling through the data...",id="loadmessage"))
                     ),
                     
                     mainPanel(
                       # EDIT: Remove once fixed
                       textOutput("test1"),
                       textOutput("test2"),
                       plotOutput("SBC_plot"),
                      tableOutput("SBC_table"),
                      textOutput("SBC_wrnings")
                     )
                   )
                 )
        ),
      tabPanel(h2("GO Term Info"), fluid = TRUE,
               sidebarLayout(
                 sidebarPanel(id = "sidebar",
                   searchInput(
                     inputId = "GO_info_search",
                     label = "Phrase or comma-separated list of GO IDs",
                     placeholder = "GO:0000001, mitochondrion, etc...",
                     btnSearch = icon("search"),
                     width = "450px"
                   ),
                   conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                    tags$div("Crawling through the data...",id="loadmessage"))
                 ),
                 mainPanel(
                   fluidRow(class = "text-center",
                     column(width = 4,
                            h1("GO Data"),
                            "From Gene Ontology Consortium, 2021-09 release",
                            br(),
                            br(),
                            downloadButton("GOInfoDL", "Download", class = "download"),
                     ),
                     column(width = 8,
                            tableOutput("GOinfo_table"),
                            br(),
                            textOutput("GOinfo_wrnings")
                            )
                   )
                 )
               )
      ),
      tabPanel(h2("Data Sources"), fluid = TRUE, align="left",
               h1("CaveCrawler Data Sources"),
               "1. Herman, A., Brandvain, Y., Weagley, J., Jeffery, W. R., 
               Keene, A. C., Kono, T., Bilandzija, H., Borowsky, R., Espinasa, 
               L., O'Quin, K., Ornelas-Garcia, C. P., Yoshizawa, M., Carlson, B., 
               Maldonado, E., Gross, J. B., Cartwright, R. A., Rohner, N., 
               Warren, W. C., and McGaugh, S. E. (2018) The role of gene flow in 
               rapid and repeated evolution of cave related traits in Mexican 
               tetra, Astyanax mexicanus. Molecular ecology, 27, 4397-4416.",
               br(),
               br(),
               "Description:",
               br(),
               "Fst, Pi, Dxy, and Tajima's D calculated using whole genomes 
               sequenced from fin clips of wild-caught Pachon (N= 9; collected 
               in 2013), Tinaja (N = 10; collected in 2002 and 2009), Molino 
               (N = 9; collected in 1994 and 2004), Rascon (N = 6; collected in 
               2013), Rio Choy (N = 9; collected in 2013). Filters were appleid 
               to SNPs and indels to remove low confidence calls. Due to high 
               heterozygosity in all individuals, the authors excluded all SNP/
               indel sites where 100% of individuals were heterozygous. See 
               publication for more details on sequencing and read filtration.
               All data was mapped to Astyanax mexicanus 1.02 assembly, Ensembl 
               93 release.",
               br(),
               br(),
               br(),
               "2. Moran, R.L., Jaggard, J.B., Roback, E.Y., Kenzior, A., 
               Rohner, N., Kowalko, J.E., Ornelas-Garcia, P., McGaugh, S.E. and 
               Keene, A.C. (2022) Hybridization underlies localized trait 
               evolution in cavefish. iScience",
               br(),
               br(),
               "Description:",
               br(),
               "The authors calculated Dxy and Fst for Chica cave fish from two 
               pools: Pool 1 (n = 5; referred as Chica1 in tables) 91m from entry 
               to cave and Pool 2 (n = 14; referred as Chica2 in tables) 10m 
               further into cave relative to Pool 1.",
               br(),
               "The authors calculated Dxy for Rascon-Pachon, Rascon-Tinaja. 
               Pachon, Tinaja, and Rascon short-read files from Herman et. al 
               2018 were downloaded from SRA. 1 Rascon and 2 Tinaja individuals 
               were excluded from D-statistic calculations due to recent hybrid 
               ancestry, and additional samples were added to bring the final 
               sample sizes up to Tinaja = 10, Pachon = 9, and Rascon = 8.",
               br(),
               "Filters were applied to SNPs and indels to remove low confidence 
               calls. See publication for more details on sequencing and read 
               filtration, as well as the websites from which gene information 
               was derived. All data was mapped to Astyanax mexicanus 1.02 
               assembly, Ensembl 93 release.",
               br(),
               br(),
               br(),
               "3. Bradic, M., Beerli, P., Garcia-de Leon, F.J., 
               Esquivel-Bobadilla, S. and Borowsky, R.L. (2012) Gene flow and 
               population structure in the Mexican blind cavefish complex 
               (Astyanax mexicanus). BMC evolutionary biology, 12, 1-17.",
               br(),
               br(),
               "Description:",
               br(),
               "Geographic coordinates of 11 cave (Pachon, Yerbaniz, Japonis, 
               Arroyo, Tinaja, Curva, Toro, Chica, Molino, Caballo Moro) and 10 
               surface populations (Subterraneo, Rio Frio, Arroyo Sarco, Chamal,
               Rio Meco, Rio Tantaon, Rio Florido, Rio Tampaon, Nacimiento del 
               Rio Santa Clara, San Rafael Los Castros, Rio Subterraneo Valley)
               were identified. These coordinates were used, in conjunction with
               the Google Maps coordinates for the Rascon and Rio Choy surface 
               populations (accessed December 2021), to generate the map on the 
               Home module.",
               br(),
               br(),
               br(),
               "4. Mack, K.L., Jaggard, J.B., Persons, J.L., Roback, E.Y., 
               Passow, C.N., Stanhope, B.A., Ferrufino, E., Tsuchiya, D., Smith,
               S.E. and Slaughter, B.D. (2021) Repeated evolution of circadian 
               clock dysregulation in cavefish populations. PLoS genetics, 17, 
               e1009642.",
               br(),
               br(),
               "Description:",
               br(),
               "Lab-born Molino, Pachon, Tinaja, and Rio Choy fish were reared
               on a 14:10 light:dark cycle. At 30 dpf, 6 whole organisms from 
               each population were sampled for RNAseq at 6 time points (144 
               samples total; average of 14,197,772 reads per sample). Genes 
               were considered rhythmic if their JTK_cycle 24 hr periodicity was 
               below an FDR cutoff of 5%. logFC was calculated for 
               Molino-Rio Choy, Pachon-Rio Choy, and Tinaja-Rio Choy. All data 
               was mapped to Astyanax mexicanus 1.02 assembly, Ensembl 93 release.",
               br(),
               br(),
               br(),
               "5. McGaugh, S.E., Passow, C.N., Jaggard, J.B., Stahl, B.A. and 
               Keene, A.C. (2020) Unique transcriptional signatures of sleep 
               loss across independently evolved cavefish populations. Journal 
               of Experimental Zoology Part B: Molecular and Developmental 
               Evolution, 334, 497-510.",
               br(),
               br(),
               "Description:",
               br(),
               "Tinaja, Molino, Pachon, and Rio Choy fish were raised on 14:10 hr
               light-dark cycle.",
               "At 29 dpf, fry were sleep deprived by shaking in Erlenmeyer 
               flasks at random intervals < 60 seconds apart throughout a period 
               of 10 hours during the night. Control fish were housed under 
               identical conditions but were not shaken.",
               "At 30 dpf, 6 whole fish per population (Tinaja, Pachon, Molino, 
               and Rio Choy) and experimental group (sleep-deprived and control) 
               were sampled for RNA-seq.",
               "The logFC in response to sleep deprivation was calculated for 
               17,187 genes in each population. Genes were labeled as 
               differentially expressed if the Benjamini-Hochberg adjusted 
               p-value was less than 0.05.",
               "All data was mapped to Astyanax mexicanus 1.02 assembly, Ensembl
               93 release.",
               br(),
               br(),
               br(),
               "6. Warren, W.C., Boggs, T.E., Borowsky, R., Carlson, B.M., 
               Ferrufino, E., Gross, J.B., Hillier, L., Hu, Z., Keene, A.C. and 
               Kenzior, A. (2021) A chromosome-level genome of Astyanax mexicanus 
               surface fish for comparing population-specific genetic differences 
               contributing to trait evolution. Nature communications, 12, 1-12.",
               br(),
               br(),
               "Description:",
               br(),
               "The authors of this study created a fully-annotated surface fish
               genome assembly. The genome assembly from this study was acquired
               for CaveCrawler using Ensembl Genome Browser, release 104.",
               br(),
               br(),
               br(),
               "7. The UniProt Consortium. (2020) UniProt: the universal protein 
               knowledgebase in 2021. Nucleic Acids Research, 49, D480-D489.",
               br(),
               br(),
               "Description:",
               br(),
               "All Mexican tetra Gene Ontology information was obtained from
               UniProtKB (Feb. 2 2021 release) using the following search phrase: ",
               br(),
               "organism:","'Astyanax mexicanus (Blind cave fish) (Astyanax fasciatus mexicanus) [7994]'","AND proteome:up000018467"
               )#,
      #tabPanel(h2("Community Resources"), fluid = TRUE, align="left",
      #         h1("Astyanax Mexicanus Community Resources", align = "center"),
      #         br(),
      #         textOutput("community_summary"),
      #         br(),
      #         tableOutput("community_table")
      #)
    )
  )

  server = function(input, output) {
    observe({
      transc_morph_choices <- c("Control", "Rio Choy")
      # Transcription Page: Update morph-selection widget to only enable
      # comparisons between current morph and morph which is NOT morph1
      updateSelectInput(session = getDefaultReactiveDomain(),
                        "morph2",
                        choices = transc_morph_choices[transc_morph_choices != input$morph1],
                        selected = tail(transc_morph_choices, 1)
      )
      # Transcription Page: Update min and max values for logFC based on data
      output$FCthresh_updater <- renderUI({
        if((input$trstat == "logFC") & (input$morph2 == "Control")){
          sec.table <- condition_control$logFC[
                      (paste(c(input$morph1, "-", input$morph2), collapse = "")
                       %in% condition_control$Comparison) & 
                        (condition_control$Condition == input$condition)]
          sliderInput("FCthresh", "...the following threshold:", 
                      min = round(min(sec.table), 2), 
                      max = round(max(sec.table), 2), value = 0, step = 0.0005)
          
        }else if((input$trstat == "logFC") & (input$morph2 != "Control")){
          sec.table <- morph1.morph2$logFC[
            grepl(input$morph1, morph1.morph2$Comparison) &
              grepl(input$morph2, morph1.morph2$Comparison)
          ]
          sliderInput("FCthresh", "...the following threshold:", 
                      min = round(min(sec.table), 2), 
                      max = round(max(sec.table), 2), value = 0, step = 0.0005)
          
        }else{
          sliderInput(inputId = "FCthresh",
                      label = "Threshold selector will appear here once you have inputted enough data to determine min and max possible logFC values",
                      min = 0, max = 0, value = 0)
        }
      })
      # Population Genetics (Stat-By-Chr Suppage): If a valid table has been
      # created, update widget for selecting which statistic to plot
      if(length(SBCT()) == 2){
        if(sum(is.na(SBCT()[[2]])) != 9){
          updateSelectInput(session = getDefaultReactiveDomain(),
                            inputId = "stat_PlotSelect",
                            label = "Visualize statistic...",
                            choices = input$sbc_statist)
        }else{
          updateSelectInput(session = getDefaultReactiveDomain(),
                            inputId = "stat_PlotSelect",
                            label = "Visualize statistic...",
                            choices = "")
      }
      }else{
        updateSelectInput(session = getDefaultReactiveDomain(),
                          inputId = "stat_PlotSelect",
                          label = "Visualize statistic...",
                          choices = "")
      }
      # Population Genetics (Stat-By-Chr Suppage): If a valid table has been
      # created, update the widget for selecting which scaffold to plot
      if(length(SBCT()) == 2){
          if(sum(is.na(SBCT()[[2]])) != 9){
          all_plots <- StatByChrGraph(SBCT()[[2]], stat_vec = input$sbc_statist)
          available_scaffs <- c()
          for(i in 1:length(all_plots)){
            available_scaffs <- append(available_scaffs,
                                       str_split(names(all_plots)[[i]], ":")[[1]][2])
          }
          updateSelectInput(session = getDefaultReactiveDomain(),
                            inputId = "scaff_PlotSelect",
                            label = "...plotted along scaffold:",
                            choices = available_scaffs)
          }else{
            updateSelectInput(session = getDefaultReactiveDomain(),
                              inputId = "scaff_PlotSelect",
                              label = "...plotted along scaffold:",
                              choices = "")
          }
      }else{
        updateSelectInput(session = getDefaultReactiveDomain(),
                          inputId = "scaff_PlotSelect",
                          label = "...plotted along scaffold:",
                          choices = "")
      }
      # Population Genetics (StatDist Subpage): Update min and max values for
      # threshold slider
      output$thrsh_slider <- renderUI({
        if(((length(input$dist_pops) >= 2) &
         ((input$dist_statist == "Fst") | (input$dist_statist == "Dxy"))) |
         ((length(input$dist_pops) < 2) &
         ((input$dist_statist != "Fst") & (input$dist_statist != "Dxy")))){
        min_max_list <- MinMax(mm_pops = input$dist_pops,
                               mm_stat = input$dist_statist,
                               in_table = stat_table)
        sliderInput("thrsh", "Threshhold statistical value: ",
                    min = round(min_max_list[[1]],2),
                    max = round(min_max_list[[2]],2), value = 0, step = 0.0005)
      }else{
        sliderInput("thrsh",
                    "Threshold selector will appear here once you have inputted enough populations for the selected statistic.",
                    min = 0, max = 0, value = 0,
                    step = 0.05)
      }
      })
    })
    
    # Home Page: Plot map of all populations
    # Plot map of Astyanax populations
    pop_map <- ggplot(data = world_map) +
      geom_polygon(aes(x = long,
                       y = lat,
                       group = group),
                   fill = "antiquewhite",
                   color = "black") +
      # get the portion of the map where your data is located at
      coord_fixed(xlim = c(min(Latit_Longit$Longitude) - 2, max(Latit_Longit$Longitude) + 2),
                  ylim = c(min(Latit_Longit$Latitude) - 1, max(Latit_Longit$Latitude) + 1),
                  ratio = 1.3) +
      # plot your points
      geom_point(data = Latit_Longit,
                 aes(y = Latitude,
                     x = Longitude,
                     colour = Morph),
                 cex = 3) +
      # change the axis labels
      xlab("Longitude") +
      ylab("Latitude") +
      # Add lines from points to data labels
      geom_label_repel(
        data = Latit_Longit,
        aes(x = Longitude, y = Latitude, label = Population),
        color = "black", fill= "white", box.padding = 2,
        max.overlaps = 50
      ) +
      # Plot title
      ggtitle("Locations of Astyanax mexicanus Populations in Mexico") +
      # Adjust title
      theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

    output$home_plot <- renderPlot(pop_map)
    # Gene Search Page: Output four tables of all data associated with the
    # entered gene search term
    GeneSearchOutput <- eventReactive(input$Gene_search, valueExpr = {
      # Tell function which tables to output based on values of logicals
      if("Position" %in% input$GSbools){
        PB <- T
      }else{
        PB <- F
      }
      if("Transcription" %in% input$GSbools){
        TB <- T
      }else{
        TB <- F
      }
      if("Population Genetics" %in% input$GSbools){
        PGB <- T
      }else{
        PGB <- F
      }
      if("GO" %in% input$GSbools){
        GOB <- T
      }else{
        GOB <- F
      }
      if(input$Gene_search == ""){
        data.frame()
      }else{
        GeneSearch(input = input$Gene_search, 
                   posBool = PB, 
                   transcBool = TB, 
                   popgenBool = PGB, 
                   GOBool = GOB, 
                   position_table = position_table, 
                   morph1.morph2 = morph1.morph2, 
                   condition_control = condition_control,
                   stat_table = stat_table, 
                   GeneToGO = GeneToGO)
      }
    })

    # Only output position table if gene search output has stuff AND position
    # table was requested
    output$GSPos_table <- renderTable({
      if(("Position" %in% input$GSbools) & (length(GeneSearchOutput()) == 5)){
        GSPos_sorted <- GeneSearchOutput()[[1]]
        # Sort table based on inputs
        if(input$GSPwhich_sort == "GSPsc"){
          if(input$GSPsc_dir == "GSPsc_hl"){
            GSPos_sorted <- GSPos_sorted[order(GSPos_sorted$Scaffold, 
                                               decreasing = T),]
          }else{
            GSPos_sorted <- GSPos_sorted[order(GSPos_sorted$Scaffold),]
          }
        }else if(input$GSPwhich_sort == "GSPsl"){
          if(input$GSPsl_dir == "GSPsl_hl"){
            GSPos_sorted <- GSPos_sorted[order(GSPos_sorted$`Start Locus`, 
                                               decreasing = T),]
          }else{
            GSPos_sorted <- GSPos_sorted[order(GSPos_sorted$`Start Locus`),]
          }
        }else if(input$GSPwhich_sort == "GSPel"){
          if(input$GSPel_dir == "GSPel_hl"){
            GSPos_sorted <- GSPos_sorted[order(GSPos_sorted$`End Locus`, 
                                               decreasing = T),]
          }else{
            GSPos_sorted <- GSPos_sorted[order(GSPos_sorted$`End Locus`),]
          }
        }
        GSPos_sorted
      }else{
        data.frame()
      }
    })
  
    # Only output transcription table if gene search output has stuff AND 
    # transcription table was requested
    output$GSTransc_table <- renderTable({
      if(("Transcription" %in% input$GSbools) &
         (length(GeneSearchOutput()) == 5)){
        unformattedTransc <- GeneSearchOutput()[[2]]
        reformattedTransc <- data.frame(
          unformattedTransc[,1:5],
          format(unformattedTransc[,6], digits = 5),
          format(unformattedTransc[,7], digits = 5),
          unformattedTransc[,8:11]
        )
        names(reformattedTransc) <- names(unformattedTransc)
        if(input$GSTwhich_sort == "logFC"){
          if(input$GSTFC_dir == "GSTFC_hl"){
            reformattedTransc <- reformattedTransc[order(reformattedTransc$logFC, 
                                               decreasing = T),]
          }else{
            reformattedTransc <- reformattedTransc[order(reformattedTransc$logFC),]
          }
        }else if(input$GSTwhich_sort == "p"){
          if(input$GSTp_dir == "GSTp_hl"){
            reformattedTransc <- reformattedTransc[order(
              as.numeric(reformattedTransc$`p-value`), 
                                               decreasing = T),]
          }else{
            reformattedTransc <- reformattedTransc[order(
              as.numeric(reformattedTransc$`p-value`)),]
          }
        }
        reformattedTransc
      }else{
        data.frame()
      }
    })
  
    # Only output population genetics data if gene search output is valid and
    # popgen data was requested
    output$GSPopgen_table <- renderTable({
      if(("Population Genetics" %in% input$GSbools) & 
         (length(GeneSearchOutput()) == 5)){
        unformattedPopgen <- GeneSearchOutput()[[3]]
        reformattedPopgen <- data.frame(
          unformattedPopgen[,1:5],
          format(unformattedPopgen[,6], digits = 5),
          unformattedPopgen[,7]
        )
        names(reformattedPopgen) <- names(unformattedPopgen)
        # Sort based on statistic values
        if(input$GSPG_dir == "GSPG_hl"){
          reformattedPopgen <- reformattedPopgen[
            order(reformattedPopgen$`Statistic Type`,
                  reformattedPopgen$`Statistic Value`, decreasing = T),]
        }else{
          reformattedPopgen <- reformattedPopgen[
            order(reformattedPopgen$`Statistic Type`,
                  reformattedPopgen$`Statistic Value`),]
        }
        reformattedPopgen
      }else{
        data.frame()
      }
    })
  
    # Only output GO data if gene search is valid and GO data was requested
    output$GSGO_table <- renderTable({
      if(("GO" %in% input$GSbools) & 
         (length(GeneSearchOutput()) == 5)){
        GeneSearchOutput()[[4]]
      }else{
        data.frame()
      }
    })
    # Gene Search Page: Output any warnings associated with the searched gene
    output$GSwarnings <- renderText({
      if(("Position" %in% input$GSbools) & (length(GeneSearchOutput()) == 5)){
        GeneSearchOutput()[[5]]
      }else{
        "No data to display yet"
      }
    })

    # Gene Search Page: Enable downloading of Gene Search Position table
    output$GSPosDL <- downloadHandler(
      filename = function() {
        paste("CaveCrawler-GeneSearch-Position-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        # If a valid table was outputted. enable downloading of that table
        if(typeof(GeneSearchOutput()) == "list"){
          GSPosDLTable <- GeneSearchOutput()[[1]]
        # If a valid table was not outputted, enable downloading of an empty df
        }else{
          GSPosDLTable <- data.frame()
        }

        write.csv(GSPosDLTable, file, row.names = F)
      }
    )
    
    # Gene Search Page: Enable downloading of Gene Search Transcription table
    output$GSTranscDL <- downloadHandler(
      filename = function() {
        paste("CaveCrawler-GeneSearch-Transcription-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        # If a valid table was outputted. enable downloading of that table
        if(typeof(GeneSearchOutput()) == "list"){
          unformattedTransc <- GeneSearchOutput()[[2]]
          GSTranscDLTable <- data.frame(
            unformattedTransc[,1:5],
            format(unformattedTransc[,6], digits = 5),
            format(unformattedTransc[,7], digits = 5),
            unformattedTransc[,8:11]
          )
          names(GSTranscDLTable) <- names(unformattedTransc)
          # If a valid table was not outputted, enable downloading of an empty df
        }else{
          GSTranscDLTable <- data.frame()
        }
        
        write.csv(GSTranscDLTable, file, row.names = F)
      }
    )
    
    # Gene Search Page: Enable downloading of Gene Search Population Genetics table
    output$GSPopgenDL <- downloadHandler(
      filename = function() {
        paste("CaveCrawler-GeneSearch-Popgen-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        # If a valid table was outputted. enable downloading of that table
        if(typeof(GeneSearchOutput()) == "list"){
          unformattedPopgen <- GeneSearchOutput()[[3]]
          GSPopgenDLTable <- data.frame(
            unformattedPopgen[,1:5],
            format(unformattedPopgen[,6], digits = 5),
            unformattedPopgen[,7]
          )
          names(GSPopgenDLTable) <- names(unformattedPopgen)
          # If a valid table was not outputted, enable downloading of an empty df
        }else{
          GSPopgenDLTable <- data.frame()
        }
        
        write.csv(GSPopgenDLTable, file, row.names = F)
      }
    )

    # Gene Search Page: Enable downloading of Gene Search GO table
    output$GSGODL <- downloadHandler(
      filename = function() {
        paste("CaveCrawler-GeneSearch-GO-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        # If a valid table was outputted. enable downloading of that table
        if(typeof(GeneSearchOutput()) == "list"){
          GSGODLTable <- GeneSearchOutput()[[4]]
          # If a valid table was not outputted, enable downloading of an empty df
        }else{
          GSGODLTable <- data.frame()
        }
        
        write.csv(GSGODLTable, file, row.names = F)
      }
    )

    # Transcription Page labels
    output$trstat_label <- renderText("Display genes whose logFC/p-value...")
    output$dir_label <- renderText("... is above/below....")
    output$pthresh_label <- renderText("... the following threshold:")
    
    # Transcription Page: Output a table with specified transcription data
    transc_table <- eventReactive(input$Transc_enter, valueExpr = {
      if(input$morph2 == "Control"){
        condit <- input$condition
      }else{
        condit <- "Between morph"
      }
      if(input$trstat == "logFC"){
        TranscTable(morph1 = input$morph1,
                    morph2 = input$morph2,
                    condition = condit,
                    direction = input$direction,
                    tr.stat = input$trstat,
                    tr.thresh = input$FCthresh,
                    GOTable = GeneToGO)
        
      }else if(input$trstat == "p"){
        TranscTable(morph1 = input$morph1,
                    morph2 = input$morph2,
                    condition = condit,
                    direction = input$direction,
                    tr.stat = input$trstat,
                    tr.thresh = input$pthresh,
                    GOTable = GeneToGO)
      }
    })
    output$transc_table_out <- renderTable({
      # Change log fold changes to have 5 decimal points
      reformattedTranscT <- data.frame(
        transc_table()[[1]][,1:4],
        format(transc_table()[[1]][,5], digits = 5),
        format(transc_table()[[1]][,6], digits = 5),
        transc_table()[[1]][,7:10]
      )
      names(reformattedTranscT) <- names(transc_table()[[1]])
      # Sort table based on p or logFC value
      if(input$TRwhich_sort == "p"){
        if(input$TRsort_p == "TRp_hl"){
          reformattedTranscT <- reformattedTranscT[order(
            as.numeric(reformattedTranscT$`p-value`), decreasing = T),]
        }else{
          reformattedTranscT <- reformattedTranscT[order(
            as.numeric(reformattedTranscT$`p-value`)),]
        }
      }else if(input$TRwhich_sort == "logFC"){
        if(input$TRsort_FC == "TRFC_hl"){
          reformattedTranscT <- reformattedTranscT[order(reformattedTranscT$logFC, 
                                                         decreasing = T),]
        }else{
          reformattedTranscT <- reformattedTranscT[order(reformattedTranscT$logFC),]
        }
      }
      
      reformattedTranscT
    })
    output$transc_wrnings_out <- renderText({transc_table()[[2]]})

    # Transcription Page: Enable downloading of Transcription table
    output$TranscDL <- downloadHandler(
      filename = function() {
        paste("CaveCrawler-Transcription-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(transc_table(), file, row.names = F)
      }
    )
    
    # Population Genetics (Distribution Suppage): If Statistic Value was
    # specified, output a table of all genes within the specified range of
    # statistic values
    SVDT <- eventReactive(input$SVDistTable_enter, valueExpr = {
      StatDistTable(input$type,
                    input$svTB,
                    input$dist_statist,
                    thresh = input$thrsh,
                    stat_table,
                    input$dist_pops)
    }
    )
    # Render table with 5 units per statistic value
    output$SVdist_tab <- renderTable(
      if(length(SVDT()) == 2){
        SVtemp_df <- data.frame(
          SVDT()[[2]][,1:2],
          format(SVDT()[[2]][,3], digits = 5),
          SVDT()[[2]][,4:7]
        )
        names(SVtemp_df) <- names(SVDT()[[2]])
        SVtemp_df
      }
    )
    output$SVdist_wrnings <- renderText(SVDT()[[1]])
    
    # Population Genetics (Distribution Suppage): If Visualize was pressed,
    # output a plot of the number of genes with each value of the specified
    # statistic
    SVDP <- eventReactive(input$SVDistPlot_enter, valueExpr = {
      StatDistPlot(stat = input$dist_statist,
                   UL = input$svTB,
                   thresh = input$thrsh,
                   stat_table,
                   pops = input$dist_pops)
    }
    )
    output$SVdist_plot_wrnings <- renderText(SVDP()[[1]])
    output$SVdist_plot <- renderPlot(SVDP()[[2]])

    # Statistic Value SubPage: Output visualize button only if data has been inputted
    observeEvent(input$SVDistTable_enter, {
      output$SVDistPlot_enter <- renderUI({
        actionButton("SVDistPlot_enter", "Visualize")
      })
    })
    
    # Statistic Value SubPage: Enable downloading of Statistic Value table
    output$SVDistDL <- downloadHandler(
      filename = function() {
        paste("CaveCrawler-Outliers-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(SVDT()[[2]], file, row.names = F)
      }
    )

    # Population Genetics (Distribution Subpage): If Gene Count was specified,
    # output the statistics associated with the indicated number of genes
    GCDT <- eventReactive(input$GCDistTable_enter, valueExpr = {
      StatDistTable(input$type,
                    input$gcTB,
                    input$dist_statist,
                    thresh = input$gene_count,
                    stat_table,
                    input$dist_pops)
    }
    )
    # Output table with 5 decimals for the statistic column
    output$GCdist_tab <- renderTable(
      if(length(GCDT()) == 2){
        GCtemp_df <- data.frame(
          GCDT()[[2]][,1:2],
          format(GCDT()[[2]][,3], digits = 5),
          GCDT()[[2]][,4:7]
        )
        names(GCtemp_df) <- names(GCDT()[[2]])
        GCtemp_df
      }
    )
    output$GCdist_wrnings <- renderText(GCDT()[[1]])

    # Gene Count SubPage: Enable downloading of Gene Count table
    output$GCDistDL <- downloadHandler(
      filename = function() {
        paste("CaveCrawler-Outliers-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(GCDT()[[2]], file, row.names = F)
      }
    )


    # Population Genetics (Stat-By-Chr Subpage): If population, statistic, and
    # GO term have been entered, output a table of associated genes
    SBCT <- eventReactive(input$GO_search, valueExpr = {
      if((length(input$sbc_statist) == 0) & (length(input$sbc_pops) != 0)){
        "Please choose at least one statistic of interest"
      }else if((length(input$sbc_statist) != 0) & (length(input$sbc_pops) == 0)){
        "Please choose at least one population of interest"
      }else if((length(input$sbc_statist) == 0) & (length(input$sbc_pops) == 0)){
        "Please choose at least one statistic and population of interest"
      }else if((length(input$sbc_statist) != 0) & (length(input$sbc_pops) != 0)){
        StatByChrTable(GOTerm = input$GO_search,
                       GeneToGO,
                       MasterGO,
                       UpperLower,
                       stat_vec = input$sbc_statist,
                       position_table,
                       stat_table,
                       pops = input$sbc_pops
        )
      }
    }
    )

    # Population Genetics (Stat-By-Chr Subpage): Only enable visualization once
    # data is entered
    observeEvent(input$GO_search, {
      if(input$GO_search != ""){
        output$stat_PlotSelect <- renderUI({
          selectInput(
            inputId = "stat_PlotSelect",
            label = "Visualize statistic...",
            choices = "",
            selected = NULL,
            multiple = FALSE
          )
        })
        output$scaff_PlotSelect <- renderUI({
          selectInput(
            inputId = "scaff_PlotSelect",
            label = "...plotted along scaffold:",
            choices = "",
            selected = NULL,
            multiple = FALSE
          )
        })
        output$SBCP_enter <- renderUI({
          actionButton("SBCP_enter","Visualize")
        })
      }
    })
    
    # Population Genetics (Stat-By-Chr Subpage): If visualize was pressed and
    # input table is NOT full of NAs, output a plot of the appropriate statistic x
    # scaffold pair
    SBCP <- eventReactive(input$SBCP_enter, valueExpr = {
      if(sum(is.na(SBCT()[[2]])) != 9){
        SBC_plots <- StatByChrGraph(SBCT()[[2]], stat_vec = input$sbc_statist)
        for(i in 1:length(SBC_plots)){
          if((str_split(names(SBC_plots)[[i]], ":")[[1]][1] == input$stat_PlotSelect)
             & (str_split(names(SBC_plots)[[i]], ":")[[1]][2] == input$scaff_PlotSelect)){
            plot_out <- SBC_plots[[i]]
          }
        }
        plot_out
      }
    })
    output$SBC_table <- renderTable(
      # First, check if anything has been inputted into function
      if(length(SBCT()) == 2){
        # If function has inputs, check if inputs yielded a valid table
        if(sum(is.na(SBCT()[[2]])) != 9){
          # If inputs did yield a valid table, output the table with adjusted
          # decimal places
          temp_df <- data.frame(
            SBCT()[[2]][,1:7],
            format(SBCT()[[2]][,8], digits = 5),
            SBCT()[[2]][,9]
          )
          names(temp_df) <- names(SBCT()[[2]])
          # EDIT: Remove after fixing lack of sorting
          output$test1 <- renderText(input$SBCwhich_sort)
          
          # Order the table based on parameters specified
          if(input$SBCwhich_sort == "SBCsp_sort"){
            # EDIT: Remove after fixing lack of sorting
            output$test2 <- renderText(input$SBCsp_dir)
            
            if(input$SBCsp_dir == "sp_asc"){
              temp_df <- temp_df[order(temp_df$Start_Position, decreasing = T),]
            }else{
              temp_df <- temp_df[order(temp_df$Start_Position),]
            }
          }else if(input$SBCwhich_sort == "SBCep_sort"){
            # EDIT: Remove after fixing lack of sorting
            output$test2 <- renderText(input$SBCep_dir)
            
            if(input$SBCep_dir == "ep_asc"){
              temp_df <- temp_df[order(temp_df$End_Position, decreasing = T),]
            }else{
              temp_df <- temp_df[order(temp_df$End_Position),]
            }
          }else if(input$SBCwhich_sort == "SBCsv_sort"){
            # EDIT: Remove after fixing lack of sorting
            output$test2 <- renderText(input$SBCsv_dir)
            
            if(input$SBCsv_dir == "sv_desc"){
              temp_df <- temp_df[order(temp_df$Statistic_Value),]
            }else{
              temp_df <- temp_df[order(temp_df$Statistic_Value, decreasing = T),]
            }
          }
          
          temp_df
        }else{
          # If inputs yielded an erroneous table, simply output an empty table
          data.frame(Data = "NA")
        }
      }else{
        data.frame(Data = c("No data to display; see explanation below"))
      })
    output$SBC_wrnings <- renderText(SBCT()[[1]])
    output$SBC_plot <- renderPlot(SBCP())

    # Statistic-by-Chromosome SubPage: Enable downloading of SBC table
    output$SBCDL <- downloadHandler(
      filename = function() {
        paste("CaveCrawler-Outliers-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(SBCT()[[2]], file, row.names = F)
      }
    )


    # GO Term Info: If GO ID or phrase was inputted, output class, lower-level
    # GO IDs, and GO terms associated with all relevant GO IDs
    GOInfoOutWarnings <- eventReactive(input$GO_info_search, valueExpr = {
      if(input$GO_info_search == ""){
        "Warnings will populate here once GO ID or phrase is inputted."
      }else{
        GOInfo(
          GO_input = input$GO_info_search,
           MasterGO,
           UpperLower
          )[[1]]
      }
    })
    GOInfoOutTable <- eventReactive(input$GO_info_search, valueExpr = {
      if(input$GO_info_search == ""){
        data.frame(`Column Name` = "Data will populate here once GO ID or phrase is inputted.")
      }else{
        GOInfo(
          GO_input = input$GO_info_search,
          MasterGO,
          UpperLower
        )[[2]]
      }
    })
    output$GOinfo_table <- renderTable(
        GOInfoOutTable()
    )
    output$GOinfo_wrnings <- renderText(
        GOInfoOutWarnings()
    )
    # GO Info Page: Enable downloading of GO Info table
    output$GOInfoDL <- downloadHandler(
      filename = function() {
        paste("CaveCrawler-GOInfo-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        # Do not write row numbers to CSV
        write.csv(GOInfoOutTable(), file, row.names = F)
      }
    )
    
    # Text summary and table for community resources page
    #output$community_summary <- renderText("The table below describes molecular tools, genetic resources, stock populations, etc. available from different labs in the Mexican tetra research community, as well as contact info for each lab")
    #Community_Data <- read.csv("data/CommunityData.csv")
    # Replace periods with spaces
    #colnames(Community_Data) <- gsub(pattern = ".", replacement = " ", 
    #                                 x = colnames(Community_Data), fixed = T)
    #output$community_table <- renderTable(Community_Data, align = 'c')
    }


# Run the application
shinyApp(ui = ui, server = server)
