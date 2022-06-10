##
tabPanel("Network Analysis",
         tabsetPanel(
           tabPanel("Community Network",
                    tags$br(),
                    sidebarLayout(
                      sidebarPanel(
                        numericInput('NW_max_dist', 'Threshold for maximal distance between nodes', value = 0.3, min=0, max=1),
                        selectizeInput('NW_color', 'Node color', choices = covariates.colorbar),
                        selectizeInput('NW_shape', 'Node shape', choices =covariates.colorbar),
                        selectizeInput('NW_label', 'Node label', choices = c('Sample ID' = '',covariates), selected = NULL),
                        
                        checkboxInput("NW_adv", "Advanced Options"),
                        
                        conditionalPanel(
                          condition = "input.NW_adv == true | input.global_adv == true",
                          selectInput("NW_distance", "Select distance type", c("Bray-Curtis" = "bray",
                                                                             "Binary"  = "jaccard",
                                                                             "Unifrac" = "unifrac"), selected = "bray"),
                          selectInput("NW_palette", "Select palette", c("Dark2" = "Dark2",
                                                                               "Set1"  = "Set1",
                                                                               "Set2" = "Set2",
                                                                               "Set3" = "Set3",
                                                                               "Paired" = "Paired",
                                                                        "Pastel" = "Pastel"), selected = "Dark2"),
                          selectInput("NW_type", "Select Network type", c("Sample based" = "samples",
                                                                               "Taxonomy based"  = "taxa"), selected = "samples")
                        ),
                        
                        # Do plot button
                        actionButton("NW_plot_btn", "Go!", class = "btn-primary"),
                        width=5
                      ),
                      
                      mainPanel(
                        fluidRow(
                          column(9,
                                 plotlyOutput("NW_plot", height="500px")
                          )
                        ),
                        width=9)
                    )
           )
         )
)
