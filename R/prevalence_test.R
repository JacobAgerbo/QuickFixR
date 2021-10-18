#' Prevalence test for a variable of a given feature 
#'
#' @param MAE A multi-assay experiment object
#' @param tax_level The classification level used for feature
#' @param variable Compare groups by binary variable e.g. 'Disease State'
#' @param prev_threshold Threshold for prevalence testing. Defaults to 0.05
#' @param filtering_cutof Cutoff for filtering prior the test. Defaults to 0.001
#' @param datatype counts, relative abundance,logcpm
#' @return A plotly object and a table of statistics
#'
#' @examples
#' data_dir = system.file('inst/extdata/Disease_Challenge.rds', package = 'QuickFixR')
#' df <- readRDS(data_dir)
#' p <- prevalence_test(df,
#'                     tax_level='Genus',
#'                     variable='Disease_State',
#'                     datatype='relabu')
#' p
#'
#' @import plotly
#' @import animalcules
#' @import decontam
#' @import dplyr
#' @import reshape2
#' @import MultiAssayExperiment
#'
#' @export
prevalence_test <- function (MAE,variable, tax_level, prev_threshold=0.05, filtering_cutof=0.001,
                            datatype = c("logcpm", "relabu", "counts"))
{
  datatype <- match.arg(datatype)
  #datatype <- "relabu"
  variable <- variable
  #variable <- "DISEASE"
  tax_level <- tax_level
  #tax_level <- "genus"
  threshold <- prev_threshold
  #threshold <- 0.05
  cutoff <- filtering_cutof
  #cutoff <- 0.001
  
  microbe <- MAE[["MicrobeGenetics"]]
  tax_table <- as.data.frame(rowData(microbe))
  sam_table <- as.data.frame(colData(microbe))
  voi <- sam_table[,variable]
  null_hypothesis = voi[1]
  counts_table <- as.data.frame(assays(microbe))[, rownames(sam_table)]
  df <- counts_table %>% upsample_counts(tax_table, tax_level) %>%  
    {
      if (datatype == "relabu") {
        counts_to_relabu(.)
      }
      else if (datatype == "logcpm") {
        counts_to_logcpm(.)
      }
      else {
        .
      }
    } %>% {
      if (sum(base::rowSums(as.matrix(.)) == 0) > 0) {
        . <- .[-which(base::rowSums(as.matrix(.)) == 0), 
        ]
      }
      else {
        .
      }
    } %>%  as.matrix() %>%  
    t()
  #
  is.null <- sam_table[,variable] == null_hypothesis
  prevalence <- decontam::isContaminant(df, method="prevalence", neg=is.null, threshold=threshold)
  colnames(prevalence) <- c("Frequence","Prevalence", "p.Frequence", "p.Prevalence","p-value", "Prevalent_to_value")
  df <- t(df)
  df.no <- df[,is.null == FALSE]
  df.no[df.no < cutoff] <- 0
  df.no[df.no >0] <- 1
  
  df.yes <- df[,is.null == TRUE]
  df.yes[df.yes < cutoff] <- 0
  df.yes[df.yes >0] <- 1
  
  df.plot <- data.frame(No=rowSums(df.no), Yes=rowSums(df.yes),
                        Prevalence=prevalence$Prevalent_to_value,
                        p_value = prevalence$`p-value`)
  
  No <- df.plot$No
  Yes <- df.plot$Yes
  
  p <- plot_ly(df.plot, x = No, 
               y = Yes, 
               mode = "markers", 
               type = "scatter", 
               color = df.plot$Prevalence,
               text = paste(rownames(df.plot), 
                            paste("p-value: ", 
                                  round(df.plot$p_value, 3), sep = ""), sep = "\n"), 
               marker = list(size = 6))
  return(list(plot = p, table = df.plot))
}
