#' Run QucikFixR shiny app
#'

#' @import assertthat
#' @import covr
#' @import lattice
#' @import DT
#' @import animalcules
#' @import limma
#' @import decontam
#' @import MultiAssayExperiment
#' @import ggplot2
#' @import plotly
#' @import vegan
#' @import magrittr
#' @import tidyverse
#' @import biomformat
#' @import shinythemes
#' @import RColorBrewer
#' @import lmerTest
#' @import performance
#' @import gt
#' @import gtExtras
#' @import boral
#' @import phyloseq
#' @import ape
#' @importFrom shinyjs addClass

#' @return The shiny app will open
#'
#' @param dev Run the applicaiton in developer mode
#'
#' @examples
#' \dontrun{
#' QuickFix()
#' }
#' @export

QuickFix <- function(dev=FALSE) {
  appDir <- system.file("shiny", package="QuickFixR")
  if (appDir == "") {
    stop("Could not find myapp. Try re-installing `mypackage`.",
         call. = FALSE)
  }
  if (dev) {
    options(shiny.autoreload=TRUE)
  }
  shiny::runApp(appDir, display.mode="normal")
}
