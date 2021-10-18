##
tabPanel("Testing",
         tabsetPanel(
           tabPanel("Prevalence Testing",
                    tags$br(),
                    sidebarLayout(
                      sidebarPanel(
                        
                        selectizeInput('prev_taxlev', 'Classification Level', choices = tax.name, selected=tax.default),
                        selectizeInput("prev_variable", "Test prevalence of features by:", covariates.two.levels),
                        
                        checkboxInput("prev_adv", "Advanced Options"),
                        
                        conditionalPanel(
                          condition = "input.prev_adv == true | input.global_adv == true",
                        numericInput('prev_p_threshold', 'Threshold for significant prevalence', 0.05, min=0, max=1),
                        numericInput('prev_filtering_cutof', 'Filtering cutoff relative abundance', 0.0001, min=0, max=1),
                        selectInput("prev_datatype", "Select data type", c("Relative Abundance" = "relabu",
                                                                                   "Counts"  = "counts",
                                                                                   "log(CPM)" = "logcpm"), selected = "relabu")
                        ),
                        
                        # Do plot button
                        actionButton("prev_plot_btn", "Go!", class = "btn-primary"),
                        actionButton("prev_table_btn", "Table"),
                        actionButton("prev_print_btn", "Print"),
                        width=3
                      ),
                      
                    mainPanel(
                        fluidRow(
                          column(9,
                                 plotlyOutput("prev_plot", height="500px"),
                                 dataTableOutput("prev_table")
                          )
                        ),
                        width=9)
                    )
                    ), 
           tabPanel("Differential Testing",
                    tags$br(),
                    sidebarLayout(
                      sidebarPanel(
                        selectizeInput('diff_taxlev', 'Classification Level', choices = tax.name, selected=tax.default),
                        selectizeInput("input_diff_condition", "Test difference of features by:", covariates.two.levels),
                        checkboxInput("diff_adv", "Advanced Options"),
                        
                        conditionalPanel(
                          condition = "input.diff_adv == true | input.global_adv == true",
                          numericInput('diff_padj_cutoff', 'Threshold for visualising adjusted p-values', 0.05, min=0, max=1),
                          numericInput('diff_min_num_filter', 'Filtering of minimal count', 5, min=0, max=10^9),
                          selectInput("diff_sig_only", "Show only significant results", c("Yeah why not?" = "Yes", "Nope, show me all!"  = "No"), selected = "No")
                        ),
                        
                        # Do plot button
                        actionButton("diff_plot_btn", "Go!", class = "btn-primary"),
                        actionButton("diff_table_btn", "Table"),
                        downloadButton('diff_print_btn', 'Print Plot'),
                        width=3
                      ),
                      
                      mainPanel(
                        fluidRow(
                          column(9,
                                 plotlyOutput("diff_plot", height="500px"),
                                 dataTableOutput("diff_table")
                          )
                        ),
                        width=9)
                    )
                    ),
             tabPanel("Linear Mixed Effect Models",
                      tags$br(),
                      sidebarLayout(
                        sidebarPanel(
                          selectizeInput('lme_taxlev', 'Classification Level', choices = tax.name, selected=tax.default, multiple=FALSE),
                          uiOutput("lme_feature"),
                          selectizeInput("lme_exp_var", "Variable from metadata which should be tested:", covariates),
                          selectizeInput("lme_random_var", "Random variable from metadata which should be taken into account:", covariates),
                          checkboxInput("lme_adv", "Advanced Options"),
                          conditionalPanel(
                            condition = "input.lme_adv == true | input.global_adv == true",
                            numericInput('lme_tolerance', 'Tolerance of overfitting for R adjustment', 10e-5, min=0, max=1),
                            selectInput("lme_datatype", "Select data type", c("Relative Abundance" = "relabu",
                                                                              "Counts"  = "counts",
                                                                              "log(CPM)" = "logcpm"), selected = "relabu")
                          ),
                          
                          # Do plot button
                          actionButton("lme_plot_btn", "Go!", class = "btn-primary"),
                          actionButton("lme_stats_btn", "Adjust", class = "btn-primary"),
                          width=3
                        ),
                        
                        mainPanel(
                          fluidRow(
                            column(9,
                                   plotOutput("lme_plot")
                                   ),
                              column(5,
                                     tableOutput("lme_stats"))
                                  ),
                          width=9)
                        )
                      ),
           tabPanel("Bayesian Ordination and Regression AnaLysis",
                    tags$br(),
                    sidebarLayout(
                      sidebarPanel(
                        selectizeInput('boral_taxlev', 'Classification Level', choices = tax.name, selected=tax.default, multiple=FALSE),
                        selectizeInput("boral_covariates", "Covariate(s) from metadata which should be tested:", covariates, multiple=TRUE),
                        selectizeInput("boral_family", "Distribution of count data which should be used for modelling", c("Normal Distribution"="normal",
                                                                                                                          "Binomial Distribution" = "binomial", 
                                                                                                                          "Counts: Poisson Distribution"  = "poisson", 
                                                                                                                          "Negative Binomial Distribution" = "negative.binomial", 
                                                                                                                          "Log Normal Distribution"="lnormal", 
                                                                                                                          "Tweedie Distribution" = "tweedie") , selected = "normal"),
                        checkboxInput("boral_adv", "Advanced Options"),
                        conditionalPanel(
                          condition = "input.boral_adv == true | input.global_adv == true",
                          selectInput("boral_datatype", "Select data type", c("Relative Abundance" = "relabu",
                                                                            "Counts"  = "counts",
                                                                            "log(CPM)" = "logcpm"), selected = "logcpm"),
                          selectInput("boral_MCMC_control", "Select size of MCMC sampling control", c("Recommended for publishing" = "High",
                                                                              "Stratch the surface"  = "Low",
                                                                              "Sneakpeak" = "DryRun"), selected = "DryRun")
                        ),
                        
                        # Do plot button
                        actionButton("boral_plot_btn", "Go!", class = "btn-primary"),
                        actionButton("boral_stats_btn", "Table", class = "btn-primary"),
                        width=3
                      ),
                      
                      mainPanel(
                        fluidRow(
                          splitLayout(cellWidths = c("25%", "25%","50%"), plotlyOutput("boral_coef_plot", height="400px"), plotlyOutput("boral_var_plot", height="400px"), plotOutput("boral_ord_plot"),
                          ),
                          column(7,
                                 dataTableOutput("boral_stats"))
                        ),
                        width=9)
                    )
           )
         )
)
           