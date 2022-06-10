# Make network
networking <- function(MAE,
                       NW_max_dist = c(),
                              NW_color = c(),
                              NW_shape = c(),
                              NW_label=c(),
                              NW_type = c("samples", "taxa"),
                              NW_distance = c("jaccard","unifrac", "bray"),
                              NW_palette = c("Dark2", "Set1", "Set2","Set3", "Paired", "Pastel")) 
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
  palette = match.arg(NW_palette)

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
  physeq <- phyloseq::phyloseq(phyloseq::otu_table(as.matrix(counts_table),taxa_are_rows=TRUE),
                               phyloseq::tax_table(as.matrix(tax_table)),
                               phyloseq::sample_data(sam_table))
  
  
  random_tree = ape::rtree(phyloseq::ntaxa(physeq), rooted=TRUE, tip.label=phyloseq::taxa_names(physeq))
  physeq = phyloseq::merge_phyloseq(physeq, random_tree)
  
  Network <- phyloseq::make_network(physeq, max.dist=max.dist, distance = distance)
  if (length(label) == 0){
    label = c("ID")
  } else{
    label = label
  }
  dark_mode <- source("https://raw.githubusercontent.com/nsgrantham/ggdark/master/R/dark_mode.R")

  dark_theme_classic <- function(base_size = 11, base_family = "", base_line_size = base_size/22,
                                 base_rect_size = base_size/22) {
    dark_mode(theme_classic(base_size, base_family, base_line_size, base_rect_size))
  }
  plot.nw <- phyloseq::plot_network(Network, physeq, type = type,
                          color = color, 
                          shape=shape,
                          line_weight = 0.5,
                          label=label) + scale_color_brewer(palette = palette) + scale_fill_brewer(palette = palette) + dark_theme_classic() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()
      )
  
  plotly.nw <- plotly::ggplotly(plot.nw, tooltip = c("colour","shape", "label"))
  
  
  return(list(plotly = plotly.nw))
}

# Test
networking(MAE,NW_max_dist = 0.3, NW_type = "samples", NW_color = "GROUP", NW_shape = "SEX", NW_distance = "unifrac", NW_palette = "Pastel")

## SHINY
## Shiny call for function
#plot
# Network Analysis
do_Network <- eventReactive(input$NW_plot_btn, {
  result <- networking(MAE = vals$MAE,
                       NW_max_dist = input$NW_max_dist,
                       NW_color = input$NW_color,
                       NW_shape = input$NW_shape,
                       NW_label = input$NW_label,
                       NW_type = input$NW_type,
                       NW_distance = input$NW_distance,
                       NW_palette = input$NW_palette)
  return(suppressWarnings(result$plotly))
})

output$NW_plot <- renderPlotly({
  p <- suppressWarnings(do_Network())
  return(suppressWarnings(p))
})