library(phyloseq)
setwd("Documents/GitHub/QuickFixR/inst/shiny/extdata/")
MAE <- readRDS("MAE.rds")

# Make network
networking <- function(MAE,
                       NW_max_dist = c(),
                              NW_color = c(),
                              NW_shape = c(),
                              NW_label=c(),
                              NW_type = c("samples", "taxa"),
                              NW_distance = c("jaccard","unifrac", "bray")) 
{
  set.seed(1)
  # test presets
  max.dist = NW_max_dist
  color=NW_color
  #tax_level = "Genus"
  shape=NW_shape
  label=NW_label
  type = match.arg(NW_type)
  distance = match.arg(NW_distance)

  ## tables from MAE
  microbe <- MAE[['MicrobeGenetics']] #double bracket subsetting is easier
  # organism x taxlev
  tax_table <- 
    as.data.frame(SummarizedExperiment::rowData(microbe)) 
  # sample x condition
  sam_table <- 
    as.data.frame(SummarizedExperiment::colData(microbe)) 
  sam_table$ID <- rownames(sam_table)
  # organism x sample
  counts_table <- 
    as.data.frame(SummarizedExperiment::assays(microbe))[,rownames(sam_table)] 
  
  #
  physeq <- phyloseq::phyloseq(otu_table(as.matrix(counts_table),taxa_are_rows=TRUE),
                     tax_table(as.matrix(tax_table)),
                     sample_data(sam_table))
  
  
  random_tree = ape::rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
  physeq = phyloseq::merge_phyloseq(physeq, random_tree)
  
  Network <- phyloseq::make_network(physeq, max.dist=max.dist, distance = distance)
  if (length(label) == 0){
    label = c("ID")
  } else{
    label = label
  }
  plot.nw <- phyloseq::plot_network(Network, physeq, type = type,
                          color = color, 
                          shape=shape,
                          line_weight = 0.5,
                          label=label)
  
  plotly.nw <- plotly::ggplotly(plot.nw)
  
  # Return results
  output<-list(print = plot.nw, 
               plot = plotly.nw )
  return(output) 
}

# Test
#test <- networking(MAE, NW_type = "samples", NW_color = "SEX", NW_shape = "GROUP", NW_distance = "unifrac")
#test$plot
#test$print

## SHINY
# Network Analysis
do_Network <- eventReactive(input$NW_plot_btn, {
  result <- networking(MAE = vals$MAE,
                              NW_max_dist = input$NW_max_dist,
                              NW_color = input$NW_color,
                              NW_shape = input$NW_shape,
                              NW_label = input$NW_label,
                              NW_type = input$NW_type,
                              NW_distance = input$NW_distance,)
  return(suppressWarnings(result$plot))
})

output$Network_plot <- renderPlotly({
  p <- suppressWarnings(do_Network())
  return(suppressWarnings(p))
})
