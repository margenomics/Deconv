
##' @name Heatmap_graph
##' @alias Heatmap_graph
##' @title Heatmap_graph
##'
##' @usage  Function that generate a heatmap with deconvolution resultant data.
##' @param df Parameter in which the matrix resulting from the deconvolution is entered, with the cell types in the rows and the samples analysed in the columns.
##' @param file Parameter to enter the name under which you want to save the resulting graphic ("EPIC_Heatmap_Exp_data_heatmap.png").
##' @param results_dir Parameter to enter the directory where the graphics will be generated (e.g. "/bicoh/nidia/Deconv")
##' @param number_format Parameter that allows to change the visualization of the numbers in the heatmap. "%.2f" for two decimals and "%.1e" for exponential notation.
##' @param width Parameter in which a value is entered that determines the width of the graph, by default it is 700.
##' @param height Parameter in which a value is entered that determines the height of the graph, by default it is 700.
##' @return Returns the generated heatmap graph.
##' @author Nidia Barco Armengol
##' @export
##'
##' @examples
##' Heatmap_graph(df= C, file = file, results_dir = results_dir, number_format = number_format, width = width_heatmap, height = height_heatmap)

Heatmap_graph<- function(df, file, results_dir, number_format, width=700, height=700){
  # A heatmap is made with the pheatmap package on an increasing color scale.
  require(pheatmap)
  df$cell_type<- NULL
  png((filename = paste(results_dir,file,sep="/")), width = width , height = height)
  pheatmap(df, display_numbers = T, color = colorRampPalette(c('white','purple'))(100), cluster_rows = F, cluster_cols = F, fontsize_number = 15, number_color = "Black", fontsize_row = 15, fontsize_col = 15, number_format = number_format)
  dev.off()
}
