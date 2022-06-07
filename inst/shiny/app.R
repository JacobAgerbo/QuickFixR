dependencies <- c("devtools","BiocManager","shiny","shinyjs",
                  "MultiAssayExperiment","ggplot2","plotly",
                  "vegan","dplyr","magrittr","biomformat",
                  "shinythemes","RColorBrewer","decontam",
                  "animalcules","limma", "broom.mixed", "lmerTest", "performance", "gt", "gtExtras","boral","ggboral", "tidyverse", "pbkrtest", "ggiraph", "hilldiv","sva","factoextra")


# Install packages not yet installed
#CRAN
installed_packages <- dependencies %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(dependencies[!installed_packages])}
#BiocManager
installed_packages <- dependencies %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
BiocManager::install(dependencies[!installed_packages])}
#Github
installed_packages <- dependencies %in% rownames(installed.packages())
if (installed_packages[21] == FALSE) {
  remotes::install_github("jthomasmock/gtExtras")}
installed_packages <- dependencies %in% rownames(installed.packages())
if (installed_packages[23] == FALSE) {
  remotes::install_github("mbedward/ggboral")}

# Packages loading
invisible(lapply(dependencies, library, character.only = TRUE))
# source helper for ui
source(file.path("helpers", "helpers.R"),  local = TRUE)
source(file.path("helpers", "ui_helpers.R"),  local = TRUE)

## 
ui <- navbarPage(
  title = paste("QuickFixR - helps scientist making stats in R"),
  id="QuickFixR",
  fluid=TRUE,
  theme = shinytheme("sandstone"),
  source(file.path("ui", "ui_01_upload.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_02_filter.R"),  local = TRUE)$value,
  #source(file.path("ui", "ui_03_batch_correction.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_04_relabu.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_05_diversity.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_06_ordination.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_07_testing.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_08_network.R"),  local = TRUE)$value
)

server <- function(input, output, session) {
  source(file.path("server", "server_01_upload.R"),  local = TRUE)$value
  source(file.path("server", "server_02_filter.R"),  local = TRUE)$value
  #source(file.path("server", "server_03_batch_correction.R"),  local = TRUE)$value
  source(file.path("server", "server_04_relabu.R"),  local = TRUE)$value
  source(file.path("server", "server_05_diversity.R"),  local = TRUE)$value
  source(file.path("server", "server_06_ordination.R"),  local = TRUE)$value
  source(file.path("server", "server_07_testing.R"),  local = TRUE)$value
  source(file.path("server", "server_08_network.R"),  local = TRUE)$value
}

shinyApp(ui = ui, server = server)
