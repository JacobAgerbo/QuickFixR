tabPanel("Upload",
         useShinyjs(),
         tags$style(appCSS),
         tags$div(
           class = "jumbotron",
           tags$div(
             class = "container",
             fluidRow(
               column(7, h1("QuickFixR"))
             ),
             p("Interactive Omics Explorative Toolkit"),
             uiOutput("tab")
             
           )
         ),
         sidebarLayout(
           sidebarPanel(
             tags$span(style="color:#72bcd4", "Application Settings"),
             checkboxInput("global_adv", "Always Show Advanced"),
             checkboxInput("upload_adv", "Alternative Upload"),
             
             conditionalPanel(condition = "input.upload_adv == false",
                              tags$span(style="color:#72bcd4", "Upload"),
                              radioButtons("uploadChoice", "",
                                           c("Count file" = "count",
                                             "Example data" = "example",
                                             "R data object" = "R data.file"
                                           ))
             ),
             
             conditionalPanel(
               condition = "input.upload_adv == true",
               tags$span(style="color:#72bcd4", "Upload"),
               
               radioButtons("uploadChoiceAdv", "",
                            c(
                              "BIOM file" = "biom",
                              "Count file with tax id" = 'countTi',
                              "PathoScope file" = "pathofiles",
                              "R data.file" = "R-id"
                              
                            ))
             ),
             conditionalPanel(
               condition = "input.upload_adv == false",
               conditionalPanel(condition = sprintf("input['%s'] == 'example'", "uploadChoice"),
                                selectInput("example_data", "Example dataset",
                                            c(
                                              "Simulated dataset" = "toy",
                                              "Disease Challenge 16S profiling" = "Disease_Challenge",
                                              "Norwegian Wild Salmon Metagenomes" = "NWS"
                                              
                                            )),
                                withBusyIndicatorUI(
                                  actionButton("upload_example",
                                               "Upload",
                                               class = "btn-primary")
                                )
               ),
               conditionalPanel(condition = sprintf("input['%s'] == 'R data.file'", "uploadChoice"),
                                fileInput("rdfile", ".rds file (required):",
                                          accept = c(
                                            ".rds"
                                          )
                                ),
                                radioButtons("rdtype", "Filetype",
                                             choices = c(
                                               rds = "rds"
                                             ),
                                             selected = "rds"
                                ),
                                withBusyIndicatorUI(
                                  actionButton("upload_animalcules",
                                               "Upload",
                                               class = "btn-primary")
                                )
                                
               ),
               
               conditionalPanel(condition = sprintf("input['%s'] == 'count'", "uploadChoice"),
                                fileInput("countsfile", "Counts file (required):",
                                          accept = c(
                                            "text/csv",
                                            "text/comma-separated-values",
                                            "text/tab-separated-values",
                                            "text/plain",
                                            ".csv",
                                            ".tsv"
                                          )
                                ),
                                fileInput("taxon.table", "Classification table file (required):",
                                          accept = c(
                                            "text/csv",
                                            "text/comma-separated-values",
                                            "text/tab-separated-values",
                                            "text/plain",
                                            ".csv",
                                            ".tsv"
                                          )
                                ),
                                fileInput("annotfile.count", "Annotation file (required):",
                                          accept = c(
                                            "text/csv",
                                            "text/comma-separated-values",
                                            "text/tab-separated-values",
                                            "text/plain",
                                            ".csv",
                                            ".tsv"
                                          )
                                ),
                                numericInput("metadata_sample_name_col_count", "Which column in metadata is sample name?",
                                             value = 1),
                                # Input: Checkbox if file has header ----
                                checkboxInput("header.count", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep.count", "Separator",
                                             choices = c(Tab = "\t",
                                                         Comma = ",",
                                                         Semicolon = ";"
                                             ),
                                             selected = ","),
                                withBusyIndicatorUI(
                                  actionButton("uploadDataCount",
                                               "Upload",
                                               class = "btn-primary")
                                )
               )
             ),
             conditionalPanel(
               condition = "input.upload_adv == true",
               conditionalPanel(condition = sprintf("input['%s'] == 'animalcules-id'", "uploadChoiceAdv"),
                                fileInput("rdfile_id", ".rds file (required):",
                                          accept = c(
                                            ".rds"
                                          )
                                ),
                                radioButtons("mae_data_type", "Choose count type",
                                             choices = c(
                                               "EM count" = "em",
                                               "Best hit" = 'hit'
                                             )
                                ),
                                withBusyIndicatorUI(
                                  actionButton("upload_mae",
                                               "Upload",
                                               class = "btn-primary")
                                )
                                
               ),
               conditionalPanel(condition = sprintf("input['%s'] == 'biom'", "uploadChoiceAdv"),
                                fileInput("biom_id", "biom file (required):",
                                          accept = c(
                                            ".biom"
                                          )
                                ),
                                helpText('Make sure the .biom file has sample metadata included.'),
                                withBusyIndicatorUI(
                                  actionButton("upload_biom",
                                               "Upload",
                                               class = "btn-primary")
                                )
                                
               ),
               conditionalPanel(condition = sprintf("input['%s'] == 'countTi'", "uploadChoiceAdv"),
                                fileInput("countsfileTi", "Counts file (required):",
                                          accept = c(
                                            "text/csv",
                                            "text/comma-separated-values",
                                            "text/tab-separated-values",
                                            "text/plain",
                                            ".csv",
                                            ".tsv"
                                          )
                                ),
                                fileInput("annotfile.countTi", "Annotation file (required):",
                                          accept = c(
                                            "text/csv",
                                            "text/comma-separated-values",
                                            "text/tab-separated-values",
                                            "text/plain",
                                            ".csv",
                                            ".tsv"
                                          )
                                ),
                                numericInput("metadata_sample_name_col_countTi", "Which column in metadata is sample name?",
                                             value = 1),
                                # Input: Checkbox if file has header ----
                                checkboxInput("header.countTi", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep.countTi", "Separator",
                                             choices = c(Tab = "\t",
                                                         Comma = ",",
                                                         Semicolon = ";"
                                             ),
                                             selected = ","),
                                withBusyIndicatorUI(
                                  actionButton("uploadDataCountTi",
                                               "Upload",
                                               class = "btn-primary")
                                )
               ),
               conditionalPanel(condition = sprintf("input['%s'] == 'pathofiles'", "uploadChoiceAdv"),
                                h5("Upload PathoScope generated .tsv files:"),
                                fileInput("countsfile.pathoscope", "PathoScope outputs (required):",
                                          multiple = TRUE,
                                          accept = c(
                                            "text/csv",
                                            "text/comma-separated-values",
                                            "text/tab-separated-values",
                                            "text/plain",
                                            ".csv",
                                            ".tsv"
                                          )
                                ),
                                fileInput("annotfile.ps", "Annotation file (required):",
                                          accept = c(
                                            "text/csv",
                                            "text/comma-separated-values",
                                            "text/tab-separated-values",
                                            "text/plain",
                                            ".csv",
                                            ".tsv"
                                          )
                                ),
                                textInput("report_suffix", "Report suffix", value = "-sam-report.tsv"),
                                numericInput("metadata_sample_name_col", "Which column in metadata is sample name?",
                                             value = 1),
                                # Input: Checkbox if file has header ----
                                checkboxInput("header.ps", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep.ps", "Separator",
                                             choices = c(Tab = "\t",
                                                         Comma = ",",
                                                         Semicolon = ";"
                                             ),
                                             selected = "\t"),
                                withBusyIndicatorUI(
                                  actionButton("uploadDataPs",
                                               "Upload",
                                               class = "btn-primary")
                                ),
                                helpText("This might take a few seconds to upload.")
                                
               )
             )
           ),
           mainPanel(
             conditionalPanel(
               condition = "input.upload_adv == true",
               conditionalPanel(condition = "input.uploadChoiceAdv === 'pathofiles'",
                                h4("Note: please click \"open in browser\" for enabling functions like multiple files upload."),
                                helpText("Counts Table: column names must be sample name"),
                                DT::dataTableOutput("contents.count"),
                                helpText("Annotation table"),
                                DT::dataTableOutput("contents.meta")
               ),
               conditionalPanel(condition = "input.uploadChoiceAdv === 'countTi'",
                                
                                tags$img(src='data_example.png', height = 120, width = 1200),
                                helpText("Counts Table: "),
                                helpText("1. Column names must be sample name"),
                                helpText("2. The first column must be Tax id"),
                                
                                DT::dataTableOutput("contents.count.2Ti"),
                                
                                helpText("Annotation table: "),
                                helpText("1. Row names must be sample name"),
                                helpText("2. The first row must sample attribute labels"),
                                
                                DT::dataTableOutput("contents.meta.2Ti")
               ),
               
               conditionalPanel(condition = "input.uploadChoiceAdv === 'biom'",
                                h5('Note: Please check http://biom-format.org/documentation/adding_metadata.html 
                                 to add sample metadate into .biom file if sample metadata is missing.'),
                                helpText("Counts Table"),
                                DT::dataTableOutput("biom.count"),
                                helpText("Annotation table"),
                                DT::dataTableOutput("biom.meta"),
                                helpText("Classification table"),
                                DT::dataTableOutput("biom.tax")
               )),
             
             conditionalPanel(
               condition = "input.upload_adv == false",
               conditionalPanel(condition = "input.uploadChoice === 'example'",
                                
                                conditionalPanel(
                                  condition = "input.example_data == 'Disease_Challenge'",
                                  h4("16S profiling"),
                                  h5("16S profiling is a real dataset containing X samples and Y microbes from gut content of rainbow trout. "),
                                  h6("Reference: Rasmussen et al. 2021. Integrative analyses of probiotics, pathogenic infections, and host immunity highlight the importance of gut microbiota in understanding disease resilience in rainbow trout (Oncorhynchus mykiss)")
                                ),
                                conditionalPanel(
                                  condition = "input.example_data == 'NWS'",
                                  h4("Norwegian Wild Salmon Metagenomes"),
                                  h5("Gut content metagenomic data from 70 adult wild Atlantic salmons sampled from five different location in northern Norway"),
                                  h6("Rasmussen et al. 2021. The hologenomic biogeography of Norwegian wild Atlantic salmon (Salmo salar)")
                                ),
                                conditionalPanel(
                                  condition = "input.example_data == 'toy'",
                                  h4("Simulated dataset"),
                                  h5("Simulated dataset is a small synthetic microbiome dataset containing 50 samples and 100 microbes. 
                     It's automatically loaded, so you could simply continue and play with the software from this point!")
                                  
                                )
                                
               ),
               conditionalPanel(condition = "input.uploadChoice === 'count'",
                                
                                tags$img(src='data_example.png', height = 120, width = 1200),
                                helpText("Counts Table: "),
                                helpText("1. Column names must be sample name"),
                                helpText("2. The first column must be microbe name"),
                                
                                DT::dataTableOutput("contents.count.2"),
                                helpText("Classification Table: "),
                                helpText("1. Column names must be Classification levels, like family, genus, species for 16S or metagenomics. 
                                         This could also be metabolic annotation or transcriptomic gene calls, right?"),
                                helpText("2. The first column must be classification ID, like ASV ID, gene calls, COG20 categories, KEGG pathways, MAG ID, feature ID etc."),
                                
                                DT::dataTableOutput("contents.taxonomy"),
                                helpText("Annotation table: "),
                                helpText("1. Row names must be sample name"),
                                helpText("2. The first row must sample attribute labels"),
                                helpText("3. Add all variables needed for your analysis"),
                                
                                DT::dataTableOutput("contents.meta.2")
               )
             ))
         )
         
)