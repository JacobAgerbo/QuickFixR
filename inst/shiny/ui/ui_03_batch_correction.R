tabPanel("Correction and Normalisation",
  tabsetPanel(
    tabPanel("Batch Correction",
             tags$br(),
             sidebarLayout(
               sidebarPanel(
                 selectizeInput("bnorm_batch", "Batch variable to be corrected for:", covariates, multiple=FALSE),
                 
                 # Do plot button
                 actionButton("batch_correct_btn", "Go!", class = "btn-primary"),
                 width=3
               ),
               
               mainPanel(
                 fluidRow(
                   splitLayout(cellWidths = c("25%", "75%"), plotOutput("bnorm_Rsq"),plotOutput("bnorm_PCA"),
                   )
                 ),
                 width=9)
             )
    )
  )
)