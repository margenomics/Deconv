
##' @name GEDIT_Deconv
##' @alias GEDIT_Deconv
##' @title GEDIT_Deconv
##'
##' @usage Function that generates the relative deconvolution graphs with the GEDIT method. This function generates the two graphs without fractionation and returns a df with the deconvolution result.
##' @param matrix Parameter in which a matrix (dataframe) is introduced, which in the first column has the names of the genes and in the rest of the columns the expressions for each gene of the different samples.
##' @param sig.matrix Parameter in which a matrix (dataframe) is introduced, with the first column containing the names of the genes and the rest of the columns containing the expression signatures of different cell types.
##' @param results_dir Parameter to enter the directory where the graphics will be generated (e.g. "/bicoh/nidia/Deconv").
##' @param height_deconv Parameter in which a value is entered that determines the height of the bar chart, by default it is 10.
##' @param width_deconv Parameter in which a value is entered that determines the width of the bar chart, by default it is 9.
##' @param height_heatmap Parameter in which a value is entered that determines the height of the heatmap, by default it is 700.
##' @param width_heatmap Parameter in which a value is entered that determines the width of the heatmap, by default it is 700.
##' @param name Parameter to enter the name you want to be included in the generated graphics (Ex. "LM22" -> "GEDIT_Deconv_LM22_plot.png"), default would be ("GEDIT_Deconv_plot.png").
##' @param number_format Parameter that allows to change the visualization of the numbers in the heatmap. ("\%.2f") for two decimals and ("\%.1e") for exponential notation.
##' @param byCond Parameter in which TRUE is introduced if we introduce a vector so that the function generates the graphs dividing the samples by its condition, by default it is FALSE.
##' @param cond Vector that assigns a condition to each sample of the data frame.
##' @return Returns the generated graphs and the df of the deconvolution.
##' @author Nidia Barco Armengol
##' @export
##'
##' @examples
##' c<- GEDIT_Deconv(matrix = matrix, sig.matrix = sig.matrix, results_dir = results_dir, byCond = TRUE, cond = fractions)

GEDIT_Deconv<- function(matrix, sig.matrix, results_dir, height_deconv= 10, width_deconv= 9, height_heatmap= 700, width_heatmap= 700, name= NULL, number_format= "%.2f", byCond= FALSE, cond, data4Tyers= NULL){
  require(usethis)
  require(devtools)
  require(gplots)
  require(data.table)
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(glmnet)
  require(RColorBrewer)

  if (!is.null(data4Tyers)){
    # This section allows you to take the mixed matrix from a document with a predefined structure.
    data4Tyers<- as.data.frame(data4Tyers)
    mean_pos<- grep("mean.", colnames(data4Tyers))
    mean_pos<- as.numeric(mean_pos)
    last_pos<- mean_pos[length(mean_pos)]
    pos<- last_pos+1
    matrix<- data4Tyers[,(pos:length(data4Tyers))]
    matrix<- 2^matrix
    gen_pos<- grep("Geneid", colnames(data4Tyers))
    gen_pos<- as.numeric(gen_pos)
    gene_names<- data4Tyers[,gen_pos]
    matrix<- cbind.data.frame(gene_names, matrix)
  }

  # Possible duplicate genes in the mixed matrix are eliminated.
  colnames(matrix)[1]<- "Gene.names"
  matrix$Gene.names<- as.character(matrix$Gene.names)
  matrix <- matrix %>%
    group_by(Gene.names) %>%
    summarise_if(is.numeric,mean) %>%
    ungroup()
  rownames(matrix)<- matrix$Gene.names
  mix<- matrix

  # Possible duplicate genes in the reference matrix are eliminated.
  colnames(sig.matrix)[1]<- "Gene.names"
  sig.matrix$Gene.names<- as.character(sig.matrix$Gene.names)
  sig.matrix <- sig.matrix %>%
    group_by(Gene.names) %>%
    summarise_if(is.numeric,mean) %>%
    ungroup()
  rownames(sig.matrix)<- sig.matrix$Gene.names
  signa<- sig.matrix

  # We take the two matrices and eliminate the genes that they do not have in common, to obtain the two matrices with the same genes to do the deconvolution.
  genmix<- rownames(mix)
  gensigna<- rownames(signa)
  genmix<- intersect(genmix,gensigna)
  gensigna<- intersect(gensigna, genmix)

  colnames(mix)[1]<- "Gene_symbol"
  mix<- as.data.frame(mix)
  mix<- mix[mix$Gene_symbol %in% genmix,] ##Cal mirar el nom de la columna amb els gens pot donar problemes
  mix<- as.matrix(mix)
  mix<- mix[,-1]
  mix<-as.data.frame(mix)
  colnames(signa)[1]<- "Gene_symbol"
  signa<- as.data.frame(signa)
  signa<- signa[signa$Gene_symbol %in% gensigna,] #Cal mirar el nom de la columna amb els gens pot donar problemes
  signa<- as.matrix(signa)
  signa<- signa[,-1]
  signa<- as.data.frame(signa)

  # The deconvolution is done with the functions that are in the GEDIT_functions.R document.
  G<- A_GEDITDecon(mix, signa)

  if(byCond==FALSE){
    # Without condition vector
    if (is.null(name)){
      file <- c("GEDIT_fractions_","plot.png")
    }else{
      title <- name
      file <- c("GEDIT_fractions_", title , "_","plot.png")
    }
    file <- paste(file, collapse = "")
    # The matrix resulting from the deconvolution is transformed so that it has the same structure as with the other methods.
    G<- t(G)
    G<- as.data.frame(G)
    G$cell_type=rownames(G)
    Deconvolution_graph(df= G, file = file, results_dir = results_dir, height = height_deconv, width = width_deconv)
    if (is.null(name)){
      file <- c("GEDIT_heatmap_","plot.png")
    }else{
      # With condition vector
      title <- name
      file <- c("GEDIT_heatmap_", title , "_","plot.png")
    }
    file <- paste(file, collapse = "")
    Heatmap_graph(df= G, file = file, results_dir = results_dir, number_format = number_format, height = height_heatmap, width = width_heatmap)
  }else{
    G<- t(G)
    G<- as.data.frame(G)
    G$cell_type=rownames(G)
    list_pos <- list()
    unics <- unique(cond)
    for (i in unics){
      positions <- c()
      for (x in 1:length(cond)){
        element <- cond[x]
        if (element==i){
          positions <- c(positions, x)
          list_pos[[i]] <- positions
        }
      }
    }

    for (x in 1:length(list_pos)){
      e <- list_pos[x]
      df <- G$cell_type
      df <- as.data.frame(df)
      df_names <- c("cell_type")
      for (i in e){
        c <- G[,i]
        df <- cbind(df, c)
        df_names <- c(df_names, colnames(G)[i])
      }
      colnames(df) <- df_names

      if (is.null(name)){
        file <- c("GEDIT_fractions_", unics[x], "_","plot.png")
      }else{
        title <- name
        file <- c("GEDIT_fractions_", title , "_", unics[x], "_","plot.png")
      }
      file <- paste(file, collapse = "")
      c$cell_type=rownames(c)
      Deconvolution_graph(df= c, file = file, results_dir = results_dir, height = height_deconv, width = width_deconv)
      if (is.null(name)){
        file <- c("GEDIT_heatmap_", unics[x], "_","plot.png")
      }else{
        title <- name
        file <- c("GEDIT_heatmap_", title , "_", unics[x], "_","plot.png")
      }
      file <- paste(file, collapse = "")
      Heatmap_graph(df= c, file = file, results_dir = results_dir, number_format = number_format, height = height_heatmap, width = width_heatmap)
    }
  }
  return(G)
}
