
##' @name Deconvolution_graph
##' @alias Deconvolution_graph
##' @title Deconvolution_graph
##'
##' @usage Function that generate the staked bar chart of different deconvolution methods.
##' @param df Parameter in which the matrix resulting from the deconvolution is entered, with the cell types in the rows and the samples analysed in the columns.
##' @param file Parameter to enter the name under which you want to save the resulting graphic ("EPIC_Deconv_Exp_data_deconv.png").
##' @param results_dir Parameter to enter the directory where the graphics will be generated (e.g. "/bicoh/nidia/Deconv").
##' @param height Parameter in which a value is entered that determines the height of graph, by default it is 10.
##' @param width Parameter in which a value is entered that determines the width of the graph, by default it is 9.
##' @param x.size This parameter allows you to select the font size of the X-axis. The default is 10.
##' @param y.size This parameter allows you to select the font size of the Y-axis. The default is 10.
##' @param l.size This parameter allows you to select the font size of the legend. The default is 10.
##' @return Returns the generated graph of the deconvolution.
##' @author Nidia Barco Armengol
##' @export
##'
##' @examples
##' Deconvolution_graph(df= C, file = file, results_dir = results_dir, height = height_deconv, width = width_deconv)

Deconvolution_graph<- function(df, file, results_dir, height= 10, width= 9, x.size= 10, y.size=10, l.size=10){
  # An applied bar chart is generated with ggplot2 and two joined colorbrewer palettes are used.
  require(gplots)
  require(ggplot2)
  require(RColorBrewer)
  df<- as.data.frame(df)
  colors <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
  getPalette = colorRampPalette(colors)
  p=df %>%
    gather(sample, fraction, -cell_type) %>%
    ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    theme(axis.text.x = element_text(size = x.size)) +
    theme(axis.text.y = element_text(size = y.size)) +
    theme(legend.text = element_text(size = l.size)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_manual(values = getPalette(length(rownames(df)))) +
    scale_x_discrete(limits = rev(levels(df)))
  ggsave(filename=paste(results_dir,file,sep="/"),dpi=300,width = width, height = height)
}
