
##' @name QuanTIseq_Deconv
##' @alias QuanTIseq_Deconv
##' @title  QuanTIseq_Deconv
##'
##' @usage Function that generates the relative deconvolution graphs with the QuanTIseq method using the TIL10 reference matrix. This function generates the two graphs without fractionation and returns a df with the deconvolution result.
##' @param matrix Parameter in which a matrix (dataframe) is introduced, which in the first column has the names of the genes and in the rest of the columns the expressions for each gene of the different samples.
##' @param arrays Parameter in which TRUE is entered if our array is an array, by default FALSE.
##' @param tumor Parameter in which TRUE is entered if our matrix is of a tumour, by default FALSE.
##' @param results_dir Parameter to enter the directory where the graphics will be generated (e.g. "/bicoh/nidia/Deconv").
##' @param height_deconv Parameter in which a value is entered that determines the height of the bar chart, by default it is 10.
##' @param width_deconv Parameter in which a value is entered that determines the width of the bar chart, by default it is 9.
##' @param height_heatmap Parameter in which a value is entered that determines the height of the heatmap, by default it is 700.
##' @param width_heatmap Parameter in which a value is entered that determines the width of the heatmap, by default it is 700.
##' @param name Parameter to enter the name you want to be included in the generated graphics (Ex. "Data" -> "QuanTIseq_Deconv_Data_plot.png"), default would be ("QuanTIseq_Deconv_TIL10_plot.png").
##' @param number_format Parameter that allows to change the visualization of the numbers in the heatmap. ("\%.2f") for two decimals and ("\%.1e") for exponential notation.
##' @param byCond Parameter in which TRUE is introduced if we introduce a vector so that the function generates the graphs dividing the samples by its condition, by default it is FALSE.
##' @param cond Vector that assigns a condition to each sample of the data frame.
##' @param x.size This parameter allows you to select the font size of the X-axis in the bar chart. The default is 10.
##' @param y.size This parameter allows you to select the font size of the Y-axis in the bar chart. The default is 10.
##' @param l.size This parameter allows you to select the font size of the legend in the bar chart. The default is 10.
##' @return Returns the generated graphs and the df of the deconvolution.
##' @author Nidia Barco Armengol
##' @export
##'
##' @examples
##' c<- QuanTIseq_Deconv(matrix = matrix, results_dir = results_dir, byCond = TRUE, cond = fractions)

QuanTIseq_Deconv<- function(matrix, arrays= FALSE, tumor=FALSE, results_dir, height_deconv= 10, width_deconv= 9, height_heatmap= 700, width_heatmap= 700, name= NULL, number_format= "%.2f", byCond= FALSE, cond, data4Tyers= NULL, x.size= 10, y.size=10, l.size=10){
  require(EPIC)
  require(usethis)
  require(devtools)
  require(gplots)
  require(data.table)
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(immunedeconv)
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

  # The preprocessing for QuanTIseq in this case is only of the mixed matrix, because the reference matrix is internal TIL10, it is not allowed to use an external reference matrix.
  cpm<- Preprocesment_EPIC_FARDEEP(matrix = matrix)

  Q<- deconvolute_quantiseq.default(mix.mat = cpm, arrays = arrays, signame = "TIL10", tumor = tumor)

  if(byCond==FALSE){
    # Without condition vector
    if (is.null(name)){
      file <- c("QuanTIseq_fractions_", "TIL10" , "_","plot.png")
    }else{
      title <- name
      file <- c("QuanTIseq_fractions_", "TIL10_" , title,"_plot.png")
    }
    file <- paste(file, collapse = "")
    Q<-Q[,-1]
    Q<-t(Q)
    Q = as.data.frame(Q)
    Q$cell_type=rownames(Q)
    Deconvolution_graph(df= Q, file = file, results_dir = results_dir, height = height_deconv, width = width_deconv, x.size= x.size, y.size= y.size, l.size= l.size)
    if (is.null(name)){
      file <- c("QuanTIseq_heatmap_", "TIL10", "_","plot.png")
    }else{
      # With condition vector
      title <- name
      file <- c("QuanTIseq_heatmap_", "TIL10_", title,"_plot.png")
    }
    file <- paste(file, collapse = "")
    Heatmap_graph(df= Q, file = file, results_dir = results_dir, number_format = number_format, height = height_heatmap, width = width_heatmap)
  }else{
    Q<-Q[,-1]
    Q<-t(Q)
    Q = as.data.frame(Q)
    Q$cell_type=rownames(Q)
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
      df <- Q$cell_type
      df <- as.data.frame(df)
      df_names <- c("cell_type")
      for (i in e){
        c <- Q[,i]
        df <- cbind(df, c)
        df_names <- c(df_names, colnames(Q)[i])
      }
      colnames(df) <- df_names

      if (is.null(name)){
        file <- c("QuanTIseq_fractions_", "TIL10_", unics[x], "_","plot.png")
      }else{
        title <- name
        file <- c("QuanTIseq_fractions_", "TIL10_", title , "_", unics[x], "_", name,"_plot.png")
      }
      file <- paste(file, collapse = "")
      c$cell_type=rownames(c)
      Deconvolution_graph(df= c, file = file, results_dir = results_dir, height = height_deconv, width = width_deconv, x.size= x.size, y.size= y.size, l.size= l.size)
      if (is.null(name)){
        file <- c("QuanTIseq_heatmap_", "TIL10_" , unics[x], "_","plot.png")
      }else{
        title <- name
        file <- c("QuanTIseq_heatmap_", "TIL10_", title , "_", unics[x], "_","plot.png")
      }
      file <- paste(file, collapse = "")
      Heatmap_graph(df= c, file = file, results_dir = results_dir, number_format = number_format, height = height_heatmap, width = width_heatmap)
    }
  }
  return(Q)
}
