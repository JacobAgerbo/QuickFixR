



tax_level = "genus"
condition ="DISEASE"
qvalue = 0 

## alfa diversity Boxplot function
hildiv_boxplot <- function (MAE, tax_level, condition, qvalue = c(0,1,2,3,4,5)){
  microbe <- MAE[["MicrobeGenetics"]]
  tax_table <- as.data.frame(SummarizedExperiment::rowData(microbe))
  sam_table <- as.data.frame(SummarizedExperiment::colData(microbe))
  
  counts_table <- as.data.frame(SummarizedExperiment::assays(microbe))[,rownames(sam_table)]
  
  counts_table %<>% upsample_counts(tax_table, tax_level)
  sam_table$richness <- hilldiv::hill_div(counts_table, qvalue = qvalue)
  colnames(sam_table)[ncol(sam_table)] <- "richness"
  colnames(sam_table)[which(colnames(sam_table) == condition)] <- "condition"
  
  
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
  
  
  if (!is.character(sam_table$condition)) {
    g <- ggplot2::ggplot(sam_table, ggplot2::aes(condition, 
                                                 richness, text = rownames(sam_table), color = condition)) + 
      ggplot2::geom_point() + ggplot2::labs(title = paste("Hill diversity between ",condition, sep = "")) +
      ggplot2::xlab(condition) +
      ggplot2::ylab(paste("Diversity q-value (",qvalue,")", sep = "")) + theme_nice()
    
  } else {
    g <- ggplot2::ggplot(sam_table, ggplot2::aes(condition, 
                                                 richness, text = rownames(sam_table), color = condition)) + 
      ggplot2::geom_point() + ggplot2::geom_boxplot() + 
      ggplot2::labs(title = paste("Hill diversity between ",condition, sep = "")) +
      ggplot2::xlab(condition) +
      ggplot2::ylab(paste("Diversity q-value (",qvalue,")", sep = "")) + theme_nice() + scale_color_brewer(palette = "Set2")
    
  }
  
  g <- ggplotly(g, tooltip = "text")
  g$elementId <- NULL
  return(g)
}

hildiv_boxplot(MAE, "genus", "DISEASE", 0)
## plot alpha diversity
plotAlphaBoxplotButton <- eventReactive(input$alpha_boxplot,{
  hildiv_boxplot(MAE = vals$MAE,
                  tax_level = input$taxl.alpha,
                  condition = input$select_alpha_div_condition,
                  qvalue = input$select_hilldiv_order)
})
output$AlphaDiversity <- renderPlotly({
  plotAlphaBoxplotButton()
})


#
# # Alpha diversity table
 do_alpha_table <- function() {
   shinyInput <- vals$shiny.input
   physeq1 <- shinyInput$pstat
   if (input$taxl.alpha !="no rank")  {
     physeq1 <- tax_glom(physeq1, input$taxl.alpha)
   }
   meta.data <- physeq1@sam_data
   meta.data$sample.name <- rownames(meta.data)
   meta.data$richness <- suppressWarnings(estimate_richness(physeq = physeq1, split = T, measures = input$select_alpha_div_method)[,1])
   colnames(meta.data)[which(colnames(meta.data) == input$select_alpha_div_condition)] <- "condition"
   rownames(meta.data) <- seq_len(nrow(meta.data))
   DT::datatable(meta.data %>% dplyr::select(sample.name, condition, richness))}
#

 
 do_hilldiv_test <- function (MAE, tax_level, condition, qvalue = c(0,1,2,3,4,5), 
                              alpha_stat = c("Wilcoxon rank sum test", "T-test", "Kruskal-Wallis")) 
 {
   microbe <- MAE[["MicrobeGenetics"]]
   tax_table <- as.data.frame(rowData(microbe))
   sam_table <- as.data.frame(colData(microbe))
   counts_table <- as.data.frame(assays(microbe))[, rownames(sam_table)]
   counts_table %<>% upsample_counts(tax_table, tax_level)
   sam_table$richness <- hilldiv::hill_div(counts_table, qvalue = qvalue)
   colnames(sam_table)[ncol(sam_table)] <- "richness"
   colnames(sam_table)[which(colnames(sam_table) == condition)] <- "condition"
   return(alpha_div_test(sam_table = sam_table, alpha_stat = alpha_stat))
 }
 
## do alpha diversity statistical test
plotAlphaBoxplotButton3 <- eventReactive(input$alpha_boxplot,{
  do_hilldiv_test(MAE = vals$MAE,
                  tax_level = input$taxl.alpha,
                  condition = input$select_alpha_div_condition,
                  qvalue = input$select_hilldiv_order,
                  alpha_stat = input$select_alpha_stat_method)
})
output$alpha.stat.test <- DT::renderDataTable({
  plotAlphaBoxplotButton3()
}, options = list(paging = TRUE, 
                  scrollX = TRUE, 
                  pageLength = 5,
                  sDom  = '<"top">t<"bottom">ip'))


## beta diversity heatmap
plotBetaHeatmapServerButton <- eventReactive(input$beta_heatmap,{
  diversity_beta_heatmap(MAE = vals$MAE,
                       tax_level = input$taxl.beta,
                       input_beta_method = input$beta_method,
                       input_bdhm_select_conditions = input$bdhm_select_conditions,
                       input_bdhm_sort_by = input$bdhm_sort_by)
})
output$BetaDiversityHeatmap <- renderPlotly({
  plotBetaHeatmapServerButton()
})


## beta diversity boxplot
plotBetaBoxplotServerButton <- eventReactive(input$beta_boxplot,{
  diversity_beta_boxplot(MAE = vals$MAE,
                       tax_level = input$taxl.beta,
                       input_beta_method = input$beta_method,
                       input_select_beta_condition = input$select_beta_condition)
})
output$BetaDiversityBoxplot <- renderPlotly({
  plotBetaBoxplotServerButton()
})


#
# # Beta diversity table
# do_beta_table <- function() {
#   shinyInput <- vals$shiny.input
#   physeq1 <- shinyInput$pstat
#
#   if (input$taxl.beta=="no rank")  {
#       if (input$select_beta_div_method == "bray"){
#       #First get otu_table and transpose it:
#       dist.matrix <- t(data.frame(otu_table(physeq1)))
#       #Then use vegdist from vegan to generate a bray distance object:
#       dist.mat <- vegdist(dist.matrix, method = "bray")
#   }else{
#       dist.mat = phyloseq::distance(physeq1, method = input$select_beta_div_method)
#   }
#
#   } else{
#     physeq2 <- tax_glom(physeq1, input$taxl.beta)
#       if (input$select_beta_div_method == "bray"){
#       #First get otu_table and transpose it:
#       dist.matrix <- t(data.frame(otu_table(physeq2)))
#       #Then use vegdist from vegan to generate a bray distance object:
#       dist.mat <- vegdist(dist.matrix, method = "bray")
#       }else{
#       dist.mat = phyloseq::distance(physeq2, method = input$select_beta_div_method)
#       }
#   }
#   dist.mat <- as.matrix(dist.mat)
#   return(dist.mat)
# }
# plotBetaBoxplotServerButton2 <- eventReactive(input$beta_heatmap,{
#   do_beta_table()
# })
# output$table.beta <- DT::renderDataTable({
#   plotBetaBoxplotServerButton2()
# }, options=list(paging = TRUE, scrollX = TRUE))
#
# # Download beta diversity table
# output$download_table_beta <- downloadHandler(
#   filename = function() { paste('Beta_diversity_table', '.csv', sep='') },
#   content = function(file) {
#       shinyInput <- vals$shiny.input
#     physeq1 <- shinyInput$pstat
#
#     if (input$taxl.beta=="no rank")  {
#       dist.mat = phyloseq::distance(physeq1, method = input$select_beta_div_method)
#     } else{
#       physeq2 <- tax_glom(physeq1, input$taxl.beta)
#       dist.mat = phyloseq::distance(physeq2, method = input$select_beta_div_method)
#     }
#     dist.mat <- as.matrix(dist.mat)
#     write.csv(data.frame(dist.mat), file)
#   }
# )

plotBetaBoxplotServerButton3 <- eventReactive(input$beta_boxplot, {
  diversity_beta_test(MAE = vals$MAE,
                    tax_level = input$taxl.beta,
                    input_beta_method = input$beta_method,
                    input_select_beta_condition = input$select_beta_condition,
                    input_select_beta_stat_method = input$select_beta_stat_method,
                    input_num_permutation_permanova = input$num_permutation_permanova)
})
output$beta.stat.test <- DT::renderDataTable({
  plotBetaBoxplotServerButton3()
}, options = list(sDom  = '<"top">t<"bottom">ip'))
