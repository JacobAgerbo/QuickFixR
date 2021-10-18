#
# Prevalence
#
# Plot
#function 
#Prevalence Testing
prev_test_plot <- function (MAE,variable, tax_level, prev_threshold, filtering_cutof,
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
  
  t <- list(
    family = "Helvetica Neue",
    size = 12,
    color = "grey40")
  
  p <- plot_ly(df.plot, x = No, 
               y = Yes, 
               mode = "markers", 
              type = "scatter", 
              color = df.plot$Prevalence,
              text = paste(rownames(df.plot), 
                           paste("p-value: ", 
                                 round(df.plot$p_value, 3), sep = ""), sep = "\n"),
              marker = list(size = 6)) %>%
  layout(title = paste("Prevalence of ", variable, sep = ""), plot_bgcolor = "#ffffff", 
         xaxis = list(title = paste0("Not ",null_hypothesis)), 
         yaxis = list(title = paste0(null_hypothesis)), font=t)
  
  return(list(plot = p))
  }

prev_test_tbl <- function (MAE,variable, tax_level, prev_threshold, filtering_cutof,
                       datatype = c("logcpm", "relabu", "counts"))
{
  datatype <- match.arg(datatype)
  variable <- variable
  tax_level <- tax_level
  threshold <- prev_threshold
  cutoff <- filtering_cutof
  microbe <- MAE[["MicrobeGenetics"]]
  tax_table <- as.data.frame(rowData(microbe))
  sam_table <- as.data.frame(colData(microbe))
  voi <- as.factor(sam_table[,variable])
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
  
  return(list(table = df.plot))
}

prev_test_print <- function (MAE,variable, tax_level, prev_threshold, filtering_cutof,
                           datatype = c("logcpm", "relabu", "counts"))
{
  datatype <- match.arg(datatype)
  variable <- variable
  tax_level <- tax_level
  threshold <- prev_threshold
  cutoff <- filtering_cutof
  microbe <- MAE[["MicrobeGenetics"]]
  tax_table <- as.data.frame(rowData(microbe))
  sam_table <- as.data.frame(colData(microbe))
  voi <- as.factor(sam_table[,variable])
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
  
  print.plot <- ggplot(data=df.plot, aes(x=No, y=yes, color=Prevalence, size=2)) +
    xlab(paste("Prevalence (Not",null_hypothesis,")", sep = "")) + 
           ylab(paste("Prevalence (",null_hypothesis,")", sep = "")) + 
                  geom_point() + theme_minimal() + 
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2")

  return(list(print = print.plot))
  
}

## Tests
#prev_test_plot(MAE, variable = "DISEASE", tax_level = "genus", prev_threshold =  0.1, filtering_cutof =  0, "relabu")
#prev_test_tbl(MAE, variable = "DISEASE", tax_level = "genus", prev_threshold =  0.1, filtering_cutof =  0, "relabu")


#Differential Abundance testing

differential_test <- function(MAE,
                                   tax_level,
                                   input_da_condition = c(),
                                   min_num_filter = 5,
                                   input_da_padj_cutoff = 0.05, sig_only = c("Yes", "No")) 
  {
  
  tax_level = tax_level
  input_da_condition = input_da_condition
  # test presets
  #tax_level = "genus"
  #input_da_condition ="DISEASE"
  input_da_condition_covariate=NULL
  #tax_level = "Genus"
  #min_num_filter = 5
  min_num_filter = min_num_filter
  #input_da_padj_cutoff = 0.05
  input_da_padj_cutoff = input_da_padj_cutoff
  sig_only <- match.arg(sig_only)
  
  ## tables from MAE
  microbe <- MAE[['MicrobeGenetics']] #double bracket subsetting is easier
  # organism x taxlev
  tax_table <- 
    as.data.frame(SummarizedExperiment::rowData(microbe)) 
  # sample x condition
  sam_table <- 
    as.data.frame(SummarizedExperiment::colData(microbe)) 
  # organism x sample
  counts_table <- 
    as.data.frame(SummarizedExperiment::assays(microbe))[,rownames(sam_table)] 
  
    # Sum counts by taxon level
    count_table_tax <- counts_table %>%
      upsample_counts(tax_table, tax_level) %>%
      counts_to_logcpm()
    colnames_tmp <- colnames(count_table_tax)
    count_table_tax <- t(apply(count_table_tax, 1, as.integer))
    colnames(count_table_tax) <- colnames_tmp
    # sam table
    sam_table %<>% df_char_to_factor()
    # filter low count microbes
    count_table_tax <- 
      count_table_tax[base::rowSums(count_table_tax) >= log10(min_num_filter),]
    
    #
    if (is.null(input_da_condition_covariate)){
      dds_formula <- 
        stats::as.formula(paste("~",input_da_condition, sep = " "))
      design <- model.matrix(dds_formula, sam_table)
    } else{
      dds_formula <- stats::as.formula(paste("~",
                                             paste(
                                               paste(input_da_condition_covariate,
                                                     collapse = " + "),
                                               input_da_condition,
                                               sep = " + "),
                                             sep = " "))
      design <- model.matrix(dds_formula, sam_table)
    }
    #print(design)
    #print(str(count_table_tax))
    fit <- limma::lmFit(count_table_tax, design)
    ebayes <- limma::eBayes(fit)
    sigtab <- limma::topTable(ebayes, 
                              adjust = "BH",
                              number = nrow(count_table_tax),
                              p.value=input_da_padj_cutoff)
    p.tab <- limma::topTable(ebayes, 
                    adjust = "BH",
                    number = nrow(count_table_tax))
  
    ## make volcano plot
    if (sig_only == "Yes"){
      de = sigtab} else {
        de = p.tab}  
    
  
    
    mycolors <- c("#1C366B","#44AB65" ,"#D1362F", "#AEA8A8")
    names(mycolors) <- c("Significant", "FoldChange", "Significant&FoldChange", "Insignificant")
    
    
    # change the grouping for the entries with significance but not a large enough Fold change
    de["group"] <- "Insignificant"
    de[which(de['adj.P.Val'] < input_da_padj_cutoff & abs(de['logFC']) < 0.5 ),"group"] <- "Significant"
    
    # change the grouping for the entries a large enough Fold change but not a low enough p value
    de[which(de['adj.P.Val'] > input_da_padj_cutoff & abs(de['logFC']) > 0.5 ),"group"] <- "FoldChange"
    
    # change the grouping for the entries with both significance and large enough fold change
    de[which(de['adj.P.Val'] < input_da_padj_cutoff & abs(de['logFC']) > 0.5 ),"group"] <- "Significant&FoldChange"
    
    de["delabel"] <- NA
    de$delabel[de$group != "Insignificant"] <- rownames(de)[de$group != "Insignificant"]
    
    de.plot <- ggplot(data=de, aes(x=logFC, y=-log10(adj.P.Val), col=group, label=delabel)) + 
      geom_point() + 
      theme_minimal() +
      geom_text() + scale_colour_manual(values = mycolors) +
      geom_vline(xintercept=c(-0.5, 0.5), col="#f2856d", linetype="dotdash") +
      geom_hline(yintercept=-log10(0.05), col="#f2856d", linetype="dotdash") +
      ylab("-log10 Adjusted P-value") +
      xlab("Log2 Fold Change") +
      theme(legend.position = "none") 
    
    # make the Plot.ly plot
    # Find and label the top peaks..
    top_peaks <- de[with(de, order(logFC, adj.P.Val)),][1:5,]
    top_peaks <- rbind(top_peaks, de[with(de, order(-logFC, adj.P.Val)),][1:5,])
    
    a <- list()
    for (i in seq_len(nrow(top_peaks))) {
      m <- top_peaks[i, ]
      a[[i]] <- list(
        x = m[["logFC"]],
        y = -log10(m[["adj.P.Val"]]),
        text = rownames(m),
        xref = "x",
        yref = "y",
        showarrow = TRUE,
        arrowhead = 0.05,
        ax = 15,
        ay = -30
      )
    }
    
    p <- plot_ly(data = de, x = de$logFC, y = -log10(de$adj.P.Val), 
                 text = paste(rownames(de), 
                              paste(paste("Adj. p-value: ", 
                                    round(de$adj.P.Val, 5), sep = ""), 
                                    sep = "\n")),
                 mode = "markers", 
                 color = de$group,
                 colors = mycolors) %>% 
      layout(title ="Volcano Plot") %>%
      layout(annotations = a) %>%
      layout(xaxis = list(autotypenumbers = 'strict', title = 'Log2 Fold Change'),
             yaxis = list(title = '-log10 Adjusted P-value'),
             xaxis = list(
               zerolinecolor = '#ffff',
               zerolinewidth = 2,
               gridcolor = 'ffff'),
             yaxis = list(
               zerolinecolor = '#ffff',
               zerolinewidth = 2,
               gridcolor = 'ffff'))
    
    colnames(de)[which(colnames(de) == "adj.P.Val")] <- "padj"
    colnames(de)[which(colnames(de) == "P.Value")] <- "pValue"
    de <- de[,which(colnames(de) %in% c("padj", "pValue"))]
    de$microbe <- rownames(de)
    rownames(de) <- seq_len(nrow(de))
    de %<>% dplyr::select(microbe, padj, pValue)
    
    # Return results
    output<-list(table = de, 
                 print = de.plot, 
                 plot = p )
    return(output) 
}

# Test
#test <- differential_test(MAE, "Genus", "Infection_State", 5, 0.05, "No")
#test$plot
#test$print
#test$table

# Correlations
## linear mixed effect models (lmes) ----

lme_cor <- function(MAE,tax_level,feature=c(), exp_var, random_var, datatype = c("logcpm", "relabu", "counts"), tolerance){
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
  
  plot.data <- df %>% 
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
    )
  
  ##
  
  # create scatter plots ----
  plot <- ggplot(plot.data, aes(x = exp_var, y = feature)) +
    geom_point(aes(color = random_var), alpha = 0.95, show.legend = FALSE, size = 3) +
    geom_smooth(method = "lm", color = "darkgray", se = FALSE) +
    facet_wrap(. ~ corr_label, nrow = 1) +
    scale_color_brewer(palette = "Set2") + 
    theme_shannon() + 
    xlab(paste(exp_var)) + 
    ylab(paste(lab," of ",feature))
  
  # return outputs ----
  return(plot)
  
}

lme_cor_adj <- function(MAE,tax_level,feature=c(), exp_var, random_var, datatype = c("logcpm", "relabu", "counts"), tolerance){
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
  stats <- data.frame(adj_corr,aov$`Pr(>F)`)
  colnames(stats) <- c("Random Effect Adjusted Correlation", "p-value")
  
  # return outputs ----
  return(stats)
  
}

#lme_cor(MAE,tax_level = "genus", feature = "Fusarium", exp_var = "AGE",random_var = "GROUP", datatype = "relabu", tolerance = 10e-5)
#lme_cor_adj(MAE,tax_level = "genus", feature = "Fusarium", exp_var = "AGE",random_var = "GROUP", datatype = "relabu", tolerance = 10e-5)


# bayesian ordination and regression analysis ----
doBORAL <- function(MAE, tax_level, boral_covariates=c(), family, datatype = c("logcpm", "relabu", "counts"), MCMC.control=c("High","Low","DryRun"))
  {
  tax_level = tax_level
  #tax_level = "genus"
  
  #boral_covariates = match.arg(boral_covariates, several.ok = TRUE)
  #boral_covariates = c("AGE","DISEASE","GROUP")
  family = family
  #family = "normal"
  datatype = match.arg(datatype)
  #datatype = "logcpm"
  MCMC.control = match.arg(MCMC.control)
  #MCMC.control = "DryRun"

  # get data from MAE ----
  microbe <- MAE[["MicrobeGenetics"]]
  tax_table <- as.data.frame(rowData(microbe))
  sam_table <- as.data.frame(colData(microbe))
  counts_table <- as.data.frame(assays(microbe))[, rownames(sam_table)]
  counts_table <- counts_table
  y <- counts_table %>% upsample_counts(tax_table, tax_level) %>%  
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
    } %>%  t() 
  
  X <- sam_table %>%
    dplyr::select(boral_covariates)


  # change chars to nums ----
  # Identify and change character columns to numeric in the covariate matrix, else BORAL will fails. Strings are alphabetically. So "A" becomes 1, "B" becomes 2, etc. 
  char_columns <- sapply(X, is.character)
  char_columns <- names(char_columns[char_columns==TRUE])
  X[ , char_columns] <- apply(X[ , char_columns,drop=F], 2,           
                              function(x) as.numeric(as.factor(x)))
  
  # define testpath ----
  testpath <- file.path(tempdir(), "jagsboralmodel.txt")
  
  # define MCMC control based on shiny settings ----
  if (MCMC.control == "High"){
    control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123)
  } else if (MCMC.control == "Disease_Challenge"){
    control = list(n.burnin = 1000, n.iteration = 4000, n.thin = 30, seed = 123)
  } else if (MCMC.control == "DryRun"){
    control = list(n.burnin = 10, n.iteration = 400, n.thin = 30, seed = 1)
  }
  
  # run boral ----
  fit_traits <- boral(y, X = X, family = family,
                      mcmc.control = control, model.name = testpath,
                      lv.control = list(num.lv = 2, type = "independent", distmat = NULL), save.model = TRUE)
  # plot model ----
    var_plot <- ggplotly(gg_varpart(fit_traits, as.percent = TRUE, label.means = FALSE))
    coef_plot <- ggplotly(gg_coefsplot(fit_traits, palette = "Greens"))
    ord_plot <- gg_lvsplot(fit_traits, include = "both") + scale_color_brewer(palette = "Dark2")
    coef_stats <- as.data.frame(summary(fit_traits)$coefficients)
    coef_stats$Feature <- rownames(coef_stats)
    return(list(ord_plot= ord_plot, 
                var_plot = var_plot, 
                coef_plot = coef_plot, 
                coef_stats = coef_stats))
}



doBORAL_ord <- function(MAE, tax_level, boral_covariates = c(), family, datatype = c("logcpm", "relabu", "counts"), MCMC.control=c("High","Low","DryRun"))
{
  tax_level = tax_level
  #tax_level = "genus"
  #boral_covariates = boral_covariates
  #boral_covariates = match.arg(boral_covariates, several.ok = TRUE)
  #boral_covariates = c("AGE","DISEASE","GROUP")
  family = family
  #family = "normal"
  datatype = match.arg(datatype)
  #datatype = "logcpm"
  MCMC.control = match.arg(MCMC.control)
  #MCMC.control = "DryRun"
  
  # get data from MAE ----
  microbe <- MAE[["MicrobeGenetics"]]
  tax_table <- as.data.frame(rowData(microbe))
  sam_table <- as.data.frame(colData(microbe))
  counts_table <- as.data.frame(assays(microbe))[, rownames(sam_table)]
  counts_table <- counts_table
  y <- counts_table %>% upsample_counts(tax_table, tax_level) %>%  
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
    } %>%  t() 
  
  X <- sam_table %>%
    dplyr::select(boral_covariates)
  
  
  # change chars to nums ----
  # Identify and change character columns to numeric in the covariate matrix, else BORAL will fails. Strings are alphabetically. So "A" becomes 1, "B" becomes 2, etc. 
  char_columns <- sapply(X, is.character)
  char_columns <- names(char_columns[char_columns==TRUE])
  X[ , char_columns] <- apply(X[ , char_columns,drop=F], 2,           
                              function(x) as.numeric(as.factor(x)))
  
  # define testpath ----
  testpath <- file.path(tempdir(), "jagsboralmodel.txt")
  
  # define MCMC control based on shiny settings ----
  if (MCMC.control == "High"){
    control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123)
  } else if (MCMC.control == "Disease_Challenge"){
    control = list(n.burnin = 1000, n.iteration = 4000, n.thin = 30, seed = 123)
  } else if (MCMC.control == "DryRun"){
    control = list(n.burnin = 10, n.iteration = 400, n.thin = 30, seed = 1)
  }
  
  # run boral ----
  fit_traits <- boral(y, X = X, family = family,
                      mcmc.control = control, model.name = testpath,
                      lv.control = list(num.lv = 2, type = "independent", distmat = NULL), save.model = TRUE)
  # plot model ----
  ord_plot <- gg_lvsplot(fit_traits, include = "both") + scale_color_brewer(palette = "Dark2")
  
  return(ord_plot)
}


# test
#doBORAL_ord(MAE, "genus", boral_covariates=c("AGE","DISEASE","SEX"), family = "normal", datatype ="logcpm", MCMC.control="DryRun")


# shiny ----

## Next to do make shiny ui -server connection for BORAL (rendering)

do_boral_var_plot <- eventReactive(input$boral_plot_btn, {
  result <- doBORAL(MAE = vals$MAE,
                           tax_level = input$boral_taxlev,
                           boral_covariates = input$boral_covariates,
                           family = input$boral_family,
                           MCMC.control = input$boral_MCMC_control,
                           datatype = input$boral_datatype)
  return(suppressWarnings(result$var_plot))
})

output$boral_var_plot <- renderPlotly({
  p <- suppressWarnings(do_boral_var_plot())
  return(suppressWarnings(p))
})

do_boral_coef_plot <- eventReactive(input$boral_plot_btn, {
  result <- doBORAL(MAE = vals$MAE,
                           tax_level = input$boral_taxlev,
                           boral_covariates = input$boral_covariates,
                           family = input$boral_family,
                           MCMC.control = input$boral_MCMC_control,
                           datatype = input$boral_datatype)
  return(suppressWarnings(result$coef_plot))
})

output$boral_coef_plot <- renderPlotly({
  p <- suppressWarnings(do_boral_coef_plot())
  return(suppressWarnings(p))
})

do_boral_ord_plot <- eventReactive(input$boral_plot_btn, {
  result <- doBORAL_ord(MAE = vals$MAE,
                           tax_level = input$boral_taxlev,
                           boral_covariates = input$boral_covariates,
                           family = input$boral_family,
                           MCMC.control = input$boral_MCMC_control,
                           datatype = input$boral_datatype)
  return(suppressWarnings(result))
})

output$boral_ord_plot <- renderPlot({
  p <- suppressWarnings(do_boral_ord_plot())
  return(suppressWarnings(p))
},height = 400,width = 600)

do_boral_stats <- eventReactive(input$boral_stats_btn, {
  result <- doBORAL(MAE = vals$MAE,
                           tax_level = input$boral_taxlev,
                           boral_covariates = input$boral_covariates,
                           family = input$boral_family,
                           MCMC.control = input$boral_MCMC_control,
                           datatype = input$boral_datatype)
  return(suppressWarnings(result$coef_stats))
})

output$boral_stats <- renderDataTable({
  p <- suppressWarnings(do_boral_stats())
  return(suppressWarnings(p))
})







  #############
 ##  SHINY  ##
#############

## Shiny call for function
#plot
do_prev_plot <- eventReactive(input$prev_plot_btn, {
  result <- prev_test_plot(MAE = vals$MAE,
                      tax_level = input$prev_taxlev,
                      variable = input$prev_variable,
                      prev_threshold = input$prev_p_threshold,
                      filtering_cutof = input$prev_filtering_cutof,
                      datatype = input$prev_datatype)
  return(suppressWarnings(result$plot))
})

output$prev_plot <- renderPlotly({
    p <- suppressWarnings(do_prev_plot())
    return(suppressWarnings(p))
})



# Table
do_prev_table <- eventReactive(input$prev_table_btn, {
  result <- prev_test_tbl(MAE = vals$MAE,
                      tax_level = input$prev_taxlev,
                      variable = input$prev_variable,
                      prev_threshold = input$prev_p_threshold,
                      filtering_cutof = input$prev_filtering_cutof,
                      datatype = input$prev_datatype)
    return(result$table)
})
output$prev_table <- renderDataTable({
    t <- do_prev_table()
    return(t)
})



## print plot

# Print ggplot pdf
do_prev_print <- eventReactive(input$prev_print_btn, {
  result <- prev_test_plot(MAE = vals$MAE,
                           tax_level = input$prev_taxlev,
                           variable = input$prev_variable,
                           prev_threshold = input$prev_p_threshold,
                           filtering_cutof = input$prev_filtering_cutof,
                           datatype = input$prev_datatype)
  return(suppressWarnings(result$print))
})


output$prev_print = downloadHandler(
  filename = function() {paste("Prevalence_plot_",variable,".pdf", sep = "")},
  content = function(file) {
    pdf(file)
    result$print
    dev.off()
  }
)



### 

# Differential testing
do_diff_plot <- eventReactive(input$diff_plot_btn, {
  result <- differential_test(MAE = vals$MAE,
                           tax_level = input$diff_taxlev,
                           input_da_condition = input$input_diff_condition,
                           #if (is.null(input$input_diff_condition_covariate)) {input_da_condition_covariate <- NULL} else {input_da_condition_covariate <- input$input_diff_condition_covariate},
                           min_num_filter = input$diff_min_num_filter,
                           input_da_padj_cutoff = input$diff_padj_cutoff,
                           sig_only = input$diff_sig_only)
  return(suppressWarnings(result$plot))
})

output$diff_plot <- renderPlotly({
  p <- suppressWarnings(do_diff_plot())
  return(suppressWarnings(p))
})

# Table
do_diff_table <- eventReactive(input$diff_table_btn, {
  result <- differential_test(MAE = vals$MAE,
                           tax_level = input$diff_taxlev,
                           input_da_condition = input$input_diff_condition,
                           #input_da_condition_covariate = input$input_diff_condition_covariate,
                           min_num_filter = input$diff_min_num_filter,
                           input_da_padj_cutoff = input$diff_padj_cutoff,
                           sig_only = input$diff_sig_only)
  return(result$table)
})
output$diff_table <- renderDataTable({
  t <- do_diff_table()
  return(t)
})

# Print ggplot pdf
do_diff_print <- eventReactive(input$diff_print_btn, {
  result <- differential_test(MAE = vals$MAE,
                              tax_level = input$diff_taxlev,
                              input_da_condition = input$input_diff_condition,
                              #if (is.null(input$input_diff_condition_covariate)) {input_da_condition_covariate <- NULL} else {input_da_condition_covariate <- input$input_diff_condition_covariate},
                              min_num_filter = input$diff_min_num_filter,
                              input_da_padj_cutoff = input$diff_padj_cutoff,
                              sig_only = input$diff_sig_only)
  return(suppressWarnings(result$print))
})


output$diff_print = downloadHandler(
  filename = function() {paste("Differential_plot_",variable,".pdf", sep = "")},
  content = function(file) {
    pdf(file)
    result$print
    dev.off()
  }
)

## 

do_lme_adj <- eventReactive(input$lme_stats_btn, {
  withBusyIndicatorServer("lme_stats_btn", {
    tavlevs <- as.list(input$lme_taxlev)
    result <- lapply(tavlevs, function(x) {
      id <- paste("lme_feature", x, sep="_")
      organisms <- input[[id]]
      lme_cor_adj(MAE = vals$MAE,
              tax_level = input$lme_taxlev,
              feature = organisms,
              exp_var = input$lme_exp_var,
              random_var = input$lme_random_var,
              datatype = input$lme_datatype,
              tolerance = input$lme_tolerance)
      
    })
    return(result)
  })
})


# Plot when button is pressed
do_lme_plot <- eventReactive(input$lme_plot_btn, {
  withBusyIndicatorServer("lme_plot_btn", {
    taxlevs <- as.list(input$lme_taxlev)
    result <- lapply(taxlevs, function(x) {
      id <- paste("lme_feature", x, sep="_")
      organisms <- input[[id]]
      
      # One organism one plot
        lme_cor(MAE = vals$MAE,
                tax_level = input$lme_taxlev,
                feature = organisms,
                exp_var = input$lme_exp_var,
                random_var = input$lme_random_var,
                datatype = input$lme_datatype,
                tolerance = input$lme_tolerance)
    })
    return(result)
  })
})




# Reaction to button pressing
output$lme_plot <- renderPlot({
  p <- do_lme_plot()
  return(p)
},height = 400,width = 1000)


output$lme_stats <- renderTable({
  t <- do_lme_adj()
  return(t)
})

#renderPlot

# Return a dynamic number of organism choices for each tax level selected
output$lme_feature <- renderUI({
  taxlevs <- as.list(input$lme_taxlev)
  inputs <- lapply(taxlevs, function(x) {
    id <- paste("lme_feature", x, sep="_")
    organisms <- unique(as.data.frame(rowData(experiments(vals$MAE)[[1]]))[,x])
    selectizeInput(id, label=x, choices=organisms, selected=organisms[1], multiple=TRUE)
  })
  return(inputs)
})