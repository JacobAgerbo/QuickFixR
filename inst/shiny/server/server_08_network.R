library(phyloseq)
setwd("Documents/GitHub/QuickFixR/inst/shiny/extdata/")
MAE <- readRDS("MAE.rds")

microbe <- MAE[["MicrobeGenetics"]]
tax_table <- as.data.frame(rowData(microbe))
sam_table <- as.data.frame(colData(microbe))
voi <- sam_table[,variable]
null_hypothesis = voi[1]
counts_table <- as.data.frame(assays(microbe))[, rownames(sam_table)]

physeq <- phyloseq(otu_table(as.matrix(counts_table),taxa_are_rows=TRUE),
                              tax_table(as.matrix(tax_table)),
                              sample_data(sam_table))

library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
physeq = merge_phyloseq(physeq, random_tree)

Network <- make_network(physeq, max.dist=0.3, distance = "unifrac")

plot.nw <- plot_network(Network, physeq, color = "SEX", 
                        shape="GROUP",
                        line_weight = 0.5,
                        label=NULL)

library(plotly)
ggplotly(plot.nw)
