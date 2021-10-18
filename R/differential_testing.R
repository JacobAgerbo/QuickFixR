#' Differential testing of abundant/expression of a given feature for a variable with two levels
#'
#' @param MAE A multi-assay experiment object
#' @param tax_level The classification level used for feature
#' @param variable Compare groups by binary variable e.g. 'Disease State'
#' @param min_num_filter Filtering of minimal number of counts. Defaults to 5
#' @param padj_cutoff Cutoff for p value adjustment, using BH. Defaults to 0.05
#' @param sig_only Show only significant values in plot
#' @return A plotly object, a ggplot object for pdf saving, and a table with statistics
#'
#' @examples
#' data_dir = system.file('inst/extdata/Disease_Challenge.rds', package = 'QuickFixR')
#' df <- readRDS(data_dir)
#' p <- differential_test(df,
#'                     tax_level='Genus',
#'                     variable='Disease_State',
#'                     sig_only="No")
#' p
#'
#' @import plotly
#' @import animalcules
#' @import limma
#' @import ggplot2
#' @import dplyr
#' @import reshape2
#' @import MultiAssayExperiment
#'
#' @export
differential_test <- function(MAE,
                              tax_level,
                              variable = c(),
                              min_num_filter = 5,
                              padj_cutoff = 0.05, sig_only = c("Yes", "No")) 
{
  
  tax_level = tax_level
  input_da_condition = variable
  # test presets
  #tax_level = "Genus"
  #input_da_condition ="Infection_State"
  input_da_condition_covariate=NULL
  #tax_level = "Genus"
  #min_num_filter = 5
  min_num_filter = min_num_filter
  #input_da_padj_cutoff = 0.05
  input_da_padj_cutoff = padj_cutoff
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
  if (nrow(sigtab) == 0){
    sigtab[1,1] <- "No differentially abundant items found!"
    colnames(sigtab) <- "result"}
  
  
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
  
  colnames(sigtab)[which(colnames(sigtab) == "adj.P.Val")] <- "padj"
  colnames(sigtab)[which(colnames(sigtab) == "P.Value")] <- "pValue"
  sigtab <- sigtab[,which(colnames(sigtab) %in% c("padj", "pValue"))]
  sigtab$microbe <- rownames(sigtab)
  rownames(sigtab) <- seq_len(nrow(sigtab))
  sigtab %<>% select(microbe, padj, pValue)
  
  # Return results
  output<-list(table = de, 
               print = de.plot, 
               plot = p )
  return(output) 
}
