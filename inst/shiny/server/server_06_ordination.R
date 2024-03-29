#
# PCA
#
# Plot
do_dimred_pca_plot <- eventReactive(input$dimred_pca_plot_btn, {
    if (input$dimred_pca_shape == "None") {shape <- NULL} else {shape <- input$dimred_pca_shape}
    if (is.na(input$dimred_pca_z)) {pcz <- NULL} else {pcz <- input$dimred_pca_z}
    result <- dimred_pca(MAE = vals$MAE,
                         tax_level = input$dimred_pca_taxlev,
                         color = input$dimred_pca_color,
                         shape = shape,
                         pcx = input$dimred_pca_x,
                         pcy = input$dimred_pca_y,
                         pcz = pcz,
                         datatype = input$dimred_pca_datatype)
    return(suppressWarnings(result$plot))
})
output$dimred_pca_plot <- renderPlotly({
    p <- suppressWarnings(do_dimred_pca_plot())
    return(suppressWarnings(p))
})
# Table
do_dimred_pca_table <- eventReactive(input$dimred_pca_table_btn, {
    if (input$dimred_pca_shape == "None") {shape <- NULL} else {shape <- input$dimred_pca_shape}
    result <- dimred_pca(MAE = vals$MAE,
                         tax_level = input$dimred_pca_taxlev,
                         color = input$dimred_pca_color,
                         shape = shape,
                         pcx = input$dimred_pca_x,
                         pcy = input$dimred_pca_y,
                         datatype = input$dimred_pca_datatype)
    return(result$table)
})
output$dimred_pca_table <- renderDataTable({
    t <- do_dimred_pca_table()
    return(t)
})

#
# PCoA
#
# Plot
do_dimred_pcoa_plot <- eventReactive(input$dimred_pcoa_plot_btn, {
    if (input$dimred_pcoa_shape == "None") {shape <- NULL} else {shape <- input$dimred_pcoa_shape}
    if (is.na(input$dimred_pcoa_z)) {axz <- NULL} else {axz <- input$dimred_pcoa_z}
    result <- dimred_pcoa(MAE = vals$MAE,
                          tax_level = input$dimred_pcoa_taxlev,
                          color = input$dimred_pcoa_color,
                          shape = shape,
                          axx = input$dimred_pcoa_x,
                          axy = input$dimred_pcoa_y,
                          axz = axz,
                          method = input$dimred_pcoa_method)
    return(result$plot)
})
output$dimred_pcoa_plot <- renderPlotly({
    p <- do_dimred_pcoa_plot()
    return(p)
})
# Table
do_dimred_pcoa_table <- eventReactive(input$dimred_pcoa_table_btn, {
    if (input$dimred_pcoa_shape == "None") {shape <- NULL} else {shape <- input$dimred_pcoa_shape}
    result <- dimred_pcoa(MAE = vals$MAE,
                          tax_level = input$dimred_pcoa_taxlev,
                          color = input$dimred_pcoa_color,
                          shape = shape,
                          axx = input$dimred_pcoa_x,
                          axy = input$dimred_pcoa_y,
                          method = input$dimred_pcoa_method)
    return(result$table)
})
output$dimred_pcoa_table <- renderDataTable({
    t <- do_dimred_pcoa_table()
    return(t)
})

#
# UMAP
#
# Plot
do_dimred_umap_plot <- eventReactive(input$dimred_umap_plot_btn, {
    if (input$dimred_umap_shape == "None") {shape <- NULL} else {shape <- input$dimred_umap_shape}
    if (is.na(input$dimred_umap_z)) {cz <- NULL} else {cz <- input$dimred_umap_z}
    result <- dimred_umap(MAE = vals$MAE,
                         tax_level = input$dimred_umap_taxlev,
                         color = input$dimred_umap_color,
                         shape = shape,
                         cx = input$dimred_umap_x,
                         cy = input$dimred_umap_y,
                         cz = cz,
                         n_neighbors = input$dimred_umap_n_neighbors,
                         metric = input$dimred_umap_metric,
                         n_epochs = input$dimred_umap_n_epochs,
                         init = input$dimred_umap_init,
                         min_dist = input$dimred_umap_min_dist,
                         datatype = input$dimred_umap_datatype)
    return(suppressWarnings(result$plot))
})
output$dimred_umap_plot <- renderPlotly({
    p <- suppressWarnings(do_dimred_umap_plot())
    return(suppressWarnings(p))
})