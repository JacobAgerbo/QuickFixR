#
# Batch Correction and Normalisation
#

batch_norm <- function(MAE, batch=c()){
  
  MAE=MAE
  
  #batch="GROUP"
  # source some repos ----
  repos = getOption("repos")
  repos[c("sva","ber", "NormalizeMets")] = "http://cran.us.r-project.org"
  options(repos = repos)
  invisible(repos)
  
  # get data ----
  microbe <- MAE[["MicrobeGenetics"]]
  tax_table <- as.data.frame(rowData(microbe))
  sam_table <- as.data.frame(colData(microbe))
  counts_table <- as.data.frame(assays(microbe))[, rownames(sam_table)]
  counts_table <- t(counts_table)
  counts_table <- as.data.frame(counts_table)
  Batch <- sam_table[,batch]
  m <- cbind(Batch,counts_table)
  
  # batch normalise, using ComBat, np-ComBat, Ber, and Ber-Bagging ----
  colnames(m)[1] <- "Batch"
  com_p <- ComBat(t(as.matrix(m[, -1])), as.factor(m[, 1]), 
                  par.prior = TRUE, prior.plots = F)
  com_n <- ComBat(t(as.matrix(m[, -1])), as.factor(m[, 1]), 
                  par.prior = FALSE, prior.plots = F)
  y <- ber(as.matrix(m[, -1]), as.factor(m[, 1]))
  y.bag <- ber_bg(as.matrix(m[, -1]), as.factor(m[, 1]), partial = TRUE, 
                  nSim = 150)
  g <- data.frame(m[, 1], y)
  g.bag <- data.frame(m[, 1], y.bag)
  pcom <- data.frame(m[1], t(com_p))
  npcom <- data.frame(m[1], t(com_n))
  corr <- sapply(2:ncol(m), function(i) {
    fit <- lm(m[, i] ~ as.factor(m[, 1]), drop.unused.levels = TRUE)
    return(summary(fit)[["adj.r.squared"]])
  })
  corr.b <- sapply(2:ncol(g), function(i) {
    fit.b <- lm(g[, i] ~ as.factor(g[, 1]), drop.unused.levels = TRUE)
    return(summary(fit.b)[["adj.r.squared"]])
  })
  corr.bag <- sapply(2:ncol(g.bag), function(i) {
    fit.b.bag <- lm(g.bag[, i] ~ as.factor(g.bag[, 1]), drop.unused.levels = TRUE)
    return(summary(fit.b.bag)[["adj.r.squared"]])
  })
  corr.cp <- sapply(2:ncol(pcom), function(i) {
    fit.cp <- lm(pcom[, i] ~ as.factor(m[, 1]), drop.unused.levels = TRUE)
    return(summary(fit.cp)[["adj.r.squared"]])
  })
  corr.ncp <- sapply(2:ncol(npcom), function(i) {
    fit.ncp <- lm(npcom[, i] ~ as.factor(m[, 1]), drop.unused.levels = TRUE)
    return(summary(fit.ncp)[["adj.r.squared"]])
  })
  
  a1 <- abs(corr[1:length(corr)])
  am1 <- max(a1)
  a2 <- abs(corr.b[1:length(corr.b)])
  am2 <- max(a2)
  a2.bag <- abs(corr.bag[1:length(corr.bag)])
  am2.bag <- max(a2.bag)
  a3 <- abs(corr.cp[1:length(corr.cp)])
  am3 <- max(a3)
  a4 <- abs(corr.ncp[1:length(corr.ncp)])
  am4 <- max(a4)
  datm <- data.frame(Dataset = c("Raw", "nonpara-ComBat", "para-ComBat", 
                                 "ber", "ber-bagging"), adjRSq = c(am1, am4, am3, am2, 
                                                                   am2.bag))
  
  # find corrected dataset with lowest adjusted R^2 value ----
  min <- min(datm$adjRSq)
  p <- prcomp(m[, -1], scale = T)
  if (min(datm$adjRSq[datm$Dataset=="Raw"]) == min){
    data.norm <- m[,-c(1)]
    p.c <- prcomp(m[, -1], scale = T)
    lab = paste("Raw data")
  } else if (min(datm$adjRSq[datm$Dataset=="nonpara-ComBat"]) == min){
    data.norm <- com_n
    ncom <- data.frame(m[1], t(com_n))
    ncom[is.na(ncom)] <- 0
    p.c <- prcomp(ncom[-1], scale = T)
    lab = paste("nonpara-ComBat")
  } else if (min(datm$adjRSq[datm$Dataset=="para-ComBat"]) == min){
    data.norm <- com_p
    pcom <- data.frame(m[1], t(com_p))
    pcom[is.na(pcom)] <- 0
    p.c <- prcomp(pcom[-1], scale = T)
    lab = paste("para-ComBat")
  } else if (min(datm$adjRSq[datm$Dataset=="ber"]) == min){
    data.norm <- y
    ycom <- data.frame(m[1], y)
    ycom[is.na(ycom)] <- 0
    p.c <- prcomp(ycom[-1], scale = T)
    lab = paste("ber")
  } else if (min(datm$adjRSq[datm$Dataset=="ber-bagging"]) == min){
    data.norm <- y.bag
    y.bag.com <- data.frame(m[1], y.bag)
    y.bag.com[is.na(y.bag.com)] <- 0
    p.c <- prcomp(y.bag.com[-1], scale = T)
    lab = paste("ber-bagging")
  }
  
  #return stuff ----
  se_mgx <-
    data.norm %>%
    base::data.matrix() %>%
    S4Vectors::SimpleList() %>%
    magrittr::set_names("MGX")
  
  se_colData <-
    sam_table %>%
    S4Vectors::DataFrame()
  
  se_rowData <-
    tax_table %>%
    base::data.frame() %>%
    dplyr::mutate_all(as.character) %>%
    #dplyr::select(superkingdom, phylum, class, order, family, genus) %>%
    S4Vectors::DataFrame()
  
  microbe_se <-
    SummarizedExperiment::SummarizedExperiment(assays = se_mgx,
                                               colData = se_colData,
                                               rowData = se_rowData)
  mae_experiments <-
    S4Vectors::SimpleList(MicrobeGenetics = microbe_se)
  
  MAE.norm <-
    MultiAssayExperiment::MultiAssayExperiment(experiments = mae_experiments,
                                               colData = se_colData)
    
  return(MAE.norm)
}


batch_norm_plot <- function(MAE, batch=c()){
  
  MAE=MAE
  
  #batch="GROUP"
  # source some repos ----
  repos = getOption("repos")
  repos[c("sva","ber", "NormalizeMets")] = "http://cran.us.r-project.org"
  options(repos = repos)
  invisible(repos)
  
  # get data ----
  microbe <- MAE[["MicrobeGenetics"]]
  tax_table <- as.data.frame(rowData(microbe))
  sam_table <- as.data.frame(colData(microbe))
  counts_table <- as.data.frame(assays(microbe))[, rownames(sam_table)]
  counts_table <- t(counts_table)
  counts_table <- as.data.frame(counts_table)
  Batch <- sam_table[,batch]
  m <- cbind(Batch,counts_table)
  
  # batch normalise, using ComBat, np-ComBat, Ber, and Ber-Bagging ----
  colnames(m)[1] <- "Batch"
  com_p <- ComBat(t(as.matrix(m[, -1])), as.factor(m[, 1]), 
                  par.prior = TRUE, prior.plots = F)
  com_n <- ComBat(t(as.matrix(m[, -1])), as.factor(m[, 1]), 
                  par.prior = FALSE, prior.plots = F)
  y <- ber(as.matrix(m[, -1]), as.factor(m[, 1]))
  y.bag <- ber_bg(as.matrix(m[, -1]), as.factor(m[, 1]), partial = TRUE, 
                  nSim = 150)
  g <- data.frame(m[, 1], y)
  g.bag <- data.frame(m[, 1], y.bag)
  pcom <- data.frame(m[1], t(com_p))
  npcom <- data.frame(m[1], t(com_n))
  corr <- sapply(2:ncol(m), function(i) {
    fit <- lm(m[, i] ~ as.factor(m[, 1]), drop.unused.levels = TRUE)
    return(summary(fit)[["adj.r.squared"]])
  })
  corr.b <- sapply(2:ncol(g), function(i) {
    fit.b <- lm(g[, i] ~ as.factor(g[, 1]), drop.unused.levels = TRUE)
    return(summary(fit.b)[["adj.r.squared"]])
  })
  corr.bag <- sapply(2:ncol(g.bag), function(i) {
    fit.b.bag <- lm(g.bag[, i] ~ as.factor(g.bag[, 1]), drop.unused.levels = TRUE)
    return(summary(fit.b.bag)[["adj.r.squared"]])
  })
  corr.cp <- sapply(2:ncol(pcom), function(i) {
    fit.cp <- lm(pcom[, i] ~ as.factor(m[, 1]), drop.unused.levels = TRUE)
    return(summary(fit.cp)[["adj.r.squared"]])
  })
  corr.ncp <- sapply(2:ncol(npcom), function(i) {
    fit.ncp <- lm(npcom[, i] ~ as.factor(m[, 1]), drop.unused.levels = TRUE)
    return(summary(fit.ncp)[["adj.r.squared"]])
  })
  
  a1 <- abs(corr[1:length(corr)])
  am1 <- max(a1)
  a2 <- abs(corr.b[1:length(corr.b)])
  am2 <- max(a2)
  a2.bag <- abs(corr.bag[1:length(corr.bag)])
  am2.bag <- max(a2.bag)
  a3 <- abs(corr.cp[1:length(corr.cp)])
  am3 <- max(a3)
  a4 <- abs(corr.ncp[1:length(corr.ncp)])
  am4 <- max(a4)
  datm <- data.frame(Dataset = c("Raw", "nonpara-ComBat", "para-ComBat", 
                                 "ber", "ber-bagging"), adjRSq = c(am1, am4, am3, am2, 
                                                                   am2.bag))
  
  pm <- ggplot(data = datm, aes(x = reorder(Dataset, -adjRSq), 
                                y = adjRSq, fill = Dataset)) + geom_bar(stat = "identity", 
                                                                        width = 1) + theme_classic()
  pm <- pm + labs(x = "Dataset (Raw and Corrected)")
  pm <- pm + labs(y = "Maximum adj.R2")
  pm <- pm + coord_flip() + scale_fill_brewer(palette = "Dark2") 
  
  # find corrected dataset with lowest adjusted R^2 value ----
  min <- min(datm$adjRSq)
  p <- prcomp(m[, -1], scale = T)
  if (min(datm$adjRSq[datm$Dataset=="Raw"]) == min){
    data.norm <- m[,-c(1)]
    p.c <- prcomp(m[, -1], scale = T)
    lab = paste("Raw data")
  } else if (min(datm$adjRSq[datm$Dataset=="nonpara-ComBat"]) == min){
    data.norm <- com_n
    ncom <- data.frame(m[1], t(com_n))
    ncom[is.na(ncom)] <- 0
    p.c <- prcomp(ncom[-1], scale = T)
    lab = paste("nonpara-ComBat")
  } else if (min(datm$adjRSq[datm$Dataset=="para-ComBat"]) == min){
    data.norm <- com_p
    pcom <- data.frame(m[1], t(com_p))
    pcom[is.na(pcom)] <- 0
    p.c <- prcomp(pcom[-1], scale = T)
    lab = paste("para-ComBat")
  } else if (min(datm$adjRSq[datm$Dataset=="ber"]) == min){
    data.norm <- y
    ycom <- data.frame(m[1], y)
    ycom[is.na(ycom)] <- 0
    p.c <- prcomp(ycom[-1], scale = T)
    lab = paste("ber")
  } else if (min(datm$adjRSq[datm$Dataset=="ber-bagging"]) == min){
    data.norm <- y.bag
    y.bag.com <- data.frame(m[1], y.bag)
    y.bag.com[is.na(y.bag.com)] <- 0
    p.c <- prcomp(y.bag.com[-1], scale = T)
    lab = paste("ber-bagging")
  }
  
  # plot composition of raw and corrected data ----
  theme_nice <- function(){
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
  
  PCA.raw <- autoplot(p, data = m, colour = 'Batch', frame = TRUE, frame.type = 'norm') + theme_nice() + scale_color_brewer(palette = "Dark2") + 
    ggtitle("Raw Data") + 
    geom_hline(yintercept=0, linetype="dotted", 
               color = "grey40", size=0.5) + 
    geom_vline(xintercept=0, linetype="dotted", color = "grey40", size=0.5)
  
  PCA.corrected <- autoplot(p.c, data = m, colour = 'Batch', frame = TRUE, frame.type = 'norm') + theme_nice() + scale_color_brewer(palette = "Dark2") + 
    ggtitle(paste("Corrected Data",lab, sep = "\n")) +
    geom_hline(yintercept=0, linetype="dotted", 
               color = "grey40", size=0.5) + 
    geom_vline(xintercept=0, linetype="dotted", color = "grey40", size=0.5)
  
  PCA.plots <- cowplot::plot_grid(PCA.raw, PCA.corrected, nrow = 1)
  
  return(PCA.plots)
}

#batch_norm(MAE, "Location")

# Filter by average relative abundance
observeEvent(input$batch_correct_btn,{
  withBusyIndicatorServer("batch_correct_btn", {
    
    vals$MAE <- batch_norm(MAE = vals$MAE,
                             batch = input$bnorm_batch)
    update_inputs(session)
  })
})

do_BC_plot <- eventReactive(input$batch_correct_btn, {
  withBusyIndicatorServer("batch_correct_btn", {
      # One organism one plot
    results <- batch_norm_plot(MAE = vals$MAE,
               batch = input$bnorm_batch)
    return(results)
  })
})


output$lme_plot <- renderPlot({
  p <- do_BC_plot()
  return(p)
  },height = 400,width = 1000)