#' Linear Mixed Effect Models with random variable 
#'
#' @param MAE A multi-assay experiment object
#' @param tax_level The classification level used for feature
#' @param variable Numeric variable for correlation, eg. 'weight'
#' @param feature A given feature from count data, which should be correlated with variable
#' @param random_var Random variable, which should be taken into account for correlation, eg. 'Location' or 'batch effect'.
#' @param tolerance Tolerance of dispersion for correlation adjustment
#' @param datatype Select datatype, like relative abundance (relabu), counts, or logcpm
#' @return A plotly object, a ggplot object for pdf saving, and a table with statistics
#'
#' @examples
#' data_dir = system.file('inst/extdata/NWS_MAG_Profiling.rds', package = 'QuickFixR')
#' df <- readRDS(data_dir)
#' p <- lme_cor(df,
#'                     tax_level='Genus',
#'                     feature = 'Mycoplasma'
#'                     variable='Weight_g',
#'                     random_var="Location",
#'                     datatype = "relabu")
#' p
#'
#' @import plotly
#' @import animalcules
#' @import performance
#' @import lmerTest
#' @import broom.mixed
#' @import ggplot2
#' @import tidyverse
#' @import MultiAssayExperiment
#'
#' @export
lme_cor <- function(MAE,tax_level,feature, exp_var, random_var, datatype = c("logcpm", "relabu", "counts"), tolerance){
  tax_level = tax_level
  exp_var = exp_var
  random_var = random_var
  datatype = match.arg(datatype)
  feature = feature
  tolerance = tolerance
  
  if (datatype == "relabu"){
    lab = "Relative Abundance"
  } else if (datatype != "relabu"){
    if (datatype == "logcpm"){
      lab = "Log Counts Per Million"
    } else {
      lab = "Counts"
    }
  }
  
  
  ## ggplot theme
  theme_shannon <- function(){
    # theme minimal creates white background and removes ticks on axes ----
    theme_minimal() +
      theme(
        # removes grid lines from plot ----
        panel.grid = element_blank(),
        # moves legend to top instead of side ----
        legend.position = "top",
        # removes title from legend, often superfluous ----
        legend.title = element_blank(),
        # creates the light gray box around the plot ----
        panel.background = element_rect(color = "#F2F2F2"),
        # creates the gray background boxes for faceting labels ----
        strip.background = element_rect(
          color = "#F2F2F2",
          fill = "#F2F2F2"
        ),
        # if using facet grid, this rotates the y text to be more readable ----
        strip.text.y = element_text(angle = 0),
        # this produces a fully left justified title ----
        plot.title.position = "plot"
      )
  }
  
  microbe <- MAE[["MicrobeGenetics"]]
  tax_table <- as.data.frame(rowData(microbe))
  sam_table <- as.data.frame(colData(microbe))
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
    } %>%  t() %>%  
    as_tibble()
  
  rownames(df) <- colnames(counts_table)
  #
  
  df <- df %>% 
    dplyr::select(feature) %>% 
    # duplicate species variable for coloring & grouping ---
    mutate(exp_var = as.numeric(sam_table[,exp_var])) %>% 
    mutate(categories = sam_table[,random_var]) %>%
    mutate(random_var = sam_table[,random_var]) %>%
    drop_na()
  
  colnames(df) <- c("feature", colnames(df[2:4]))
  level_df <- c(levels(as.factor(df$categories)))
  
  plot <- df %>% 
    # add stack for all species to be analyzed together ----
  bind_rows(df %>% mutate(categories = "All")) %>% 
    # now examine by 3 species plus all ----
  group_by(categories) %>% 
    nest() %>%
    # within each group, compute base n and correlation ----
  mutate(
    base_n = map_int(data, nrow),
    corr = map(data, ~ cor.test(x = .x$exp_var, y = .x$feature) %>% broom::tidy())
  ) %>% 
    ungroup() %>% 
    # bring results back to raw data ----
  unnest(c(data, corr)) %>% 
    mutate(
      # create ordered facet label for plotting ----
      categories = fct_relevel(categories, level_df, "All"),
      corr_label =  glue::glue("{categories}\nn = {base_n}\n r = {scales::number(estimate, 0.01)}"),
      corr_label = fct_reorder(as.character(corr_label), as.numeric(categories))
    ) %>% 
    # create scatter plots ----
  ggplot(aes(x = exp_var, y = feature)) +
    geom_point(aes(color = random_var), alpha = 0.85, show.legend = FALSE) +
    geom_smooth(method = "lm", color = "darkgray", se = FALSE) +
    facet_wrap(. ~ corr_label, ncol = 4) +
    scale_color_brewer(palette = "Set2") + 
    theme_shannon() + 
    xlab(paste(exp_var)) + 
    ylab(paste(lab," of ",feature))
  
  
  plot <- ggplotly(plot)
  
  # estimate mixed model ----
  mixed_model <- lmerTest::lmer(feature ~ exp_var + (1 | random_var), df)
  
  # retrieve sign of coefficient ----
  coef_sign <- mixed_model %>% 
    broom.mixed::tidy() %>% 
    filter(term == "exp_var") %>% 
    pull(estimate) %>% 
    sign()
  
  # retrieve r2 measure ----
  r2_by_group <- performance::r2_nakagawa(mixed_model, by_group = TRUE, tolerance = tolerance)$R2[1]
  # compute adjusted correlation ----
  adj_corr <- coef_sign * sqrt(r2_by_group)
  aov <- anova(mixed_model, ddf="Kenward-Roger")
  aov$`Pr(>F)`
  stats <- c(adj_corr,aov$`Pr(>F)`)
  names(stats) <- c("Random Effect Adjusted Correlation", "p-value")
  return(list(plot = plot, stats = stats))  
  
}