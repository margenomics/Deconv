
##' @name CIBERSORT_Deconv
##' @alias CIBERSORT_Deconv
##' @title CIBERSORT_Deconv
##'
##' @usage Function that generates the relative or absolute deconvolution plots with the CIBERSORT method, with the non-fractional plot, returns a df resulting from the deconvolution and two graphs of deconvolution data.
##' @param matrix Parameter in which a matrix (in .txt, .csv or .tsv format) is introduced, which in the first column has the names of the genes and in the rest of the columns the expressions for each gene of the different samples.
##' @param sig.matrix Parameter in which a matrix (in .txt, .csv or .tsv format) is introduced, with the first column containing the names of the genes and the rest of the columns containing the expression signatures of different cell types.
##' @param method Parameter where you enter "abs" if you want to do the absolute deconvolution of the matrix, by default it is "rel".
##' @param results_dir Parameter to enter the directory where the graphics will be generated (e.g. "/bicoh/nidia/Deconv").
##' @param height_deconv Parameter in which a value is entered that determines the height of the bar chart, by default it is 10.
##' @param width_deconv Parameter in which a value is entered that determines the width of the bar chart, by default it is 9.
##' @param height_heatmap Parameter in which a value is entered that determines the height of the heatmap, by default it is 700.
##' @param width_heatmap Parameter in which a value is entered that determines the width of the heatmap, by default it is 700.
##' @param name Parameter to enter the name you want to be included in the generated graphics (Ex. "LM22" -> "CIBERSORT_Deconv_LM22_rel_plot.png"), default would be ("CIBERSORT_Deconv_(sig.matrix)_rel_plot.png").
##' @param number_format Parameter that allows to change the visualization of the numbers in the heatmap. ("\%.2f") for two decimals and ("\%.1e") for exponential notation.
##' @param byCond Parameter in which TRUE is introduced if we introduce a vector so that the function generates the graphs dividing the samples by its condition, by default it is FALSE.
##' @param cond Vector that assigns a condition to each sample of the data frame.
##' @return Returns the generated graphs and the df of the deconvolution.
##' @author Nidia Barco Armengol
##' @export
##'
##' @examples
##' c<- CIBERSORT_Deconv(matrix = matrix, sig.matrix = sig.matrix, results_dir = results_dir, cibersortpath = path_ciber, by_Cond = TRUE, cond = fractions, method = "abs")

CIBERSORT_Deconv<- function(matrix, sig.matrix, method= "rel", results_dir, height_deconv= 10, width_deconv= 9, height_heatmap= 700, width_heatmap= 700, name= FALSE, number_format= "%.2f", by_Cond= FALSE, cond, cibersortpath){
  require(usethis)
  require(devtools)
  require(gplots)
  require(data.table)
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(immunedeconv)
  require(RColorBrewer)
  source(file.path(cibersortpath, "CIBERSORT.R"))
  if (method=="rel"){
    # The deconvolution is performed and gives results over 1.
    res_ciber = CIBERSORT(sig_matrix=sig.matrix,
                          mixture_file = matrix,
                          QN = F,perm=100)
    res_ciber = as.data.frame(t(res_ciber))
    res_ciber$cell_type=rownames(res_ciber)
    C <- res_ciber[1:(length(rownames(res_ciber))-3),]
  }else{
    # The deconvolution is performed and gives quantitative results.
    res_ciber = CIBERSORT(sig_matrix=sig.matrix,
                          mixture_file = matrix,
                          QN = F,perm=100,absolute = T,abs_method='no.sumto1')
    res_ciber = as.data.frame(t(res_ciber))
    res_ciber$cell_type=rownames(res_ciber)
    C <- res_ciber[1:(length(rownames(res_ciber))-4),]
  }
  if (by_Cond==FALSE){
    # If there is no condition it does not loop to generate several graphs.
    if (name==FALSE){
      title <- gsub(".txt","_",sig.matrix)
      file <- c("CIBERSORT_fractions_", title , "_" , method , "_","plot.png")
    }else{
      title <- name
      file <- c("CIBERSORT_fractions_", title , "_" , method , "_","plot.png")
    }
    file <- paste(file, collapse = "")
    C$cell_type=rownames(C)
    Deconvolution_graph(df= C, file = file, results_dir = results_dir, height = height_deconv, width = width_deconv)
    if (name==FALSE){
      title <- gsub(".txt","_",sig.matrix)
      file <- c("CIBERSORT_heatmap_", title , "_" , method , "_","plot.png")
    }else{
      title <- name
      file <- c("CIBERSORT_heatmap", title , "_" , method , "_","plot.png")
    }
    file <- paste(file, collapse = "")
    Heatmap_graph(df= C, file = file, results_dir = results_dir, number_format = number_format, width = width_heatmap, height = height_heatmap)
  }else{
    # If there is a condition vector, the data enters a llop that divides the df and generates as many graphs of each type as there are conditions in the vector.
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
      df <- rownames(C)
      df <- as.data.frame(df)
      df_names <- c("cell_type")
      for (i in e){
        c <- C[,i]
        df <- cbind(df, c)
        df_names <- c(df_names, colnames(C)[i])
      }
      colnames(df) <- df_names

      if (name==FALSE){
        title <- gsub(".txt","_",sig.matrix)
        file <- c("CIBERSORT_fractions_", title ,"_", method, "_", x, "_","plot.png")
      }else{
        title <- name
        file <- c("CIBERSORT_fractions_", title ,"_", method, "_", x, "_","plot.png")
      }
      file <- paste(file, collapse = "")
      c$cell_type=rownames(c)
      Deconvolution_graph(df= c, file = file, results_dir = results_dir, height = height_deconv, width = width_deconv)
      if (name==FALSE){
        title <- gsub(".txt","_",sig.matrix)
        file <- c("CIBERSORT_heatmap_", title ,"_", method, "_", x, "_","plot.png")
      }else{
        title <- name
        file <- c("CIBERSORT_heatmap_", title ,"_", method, "_", x, "_","plot.png")
      }
      file <- paste(file, collapse = "")
      Heatmap_graph(df= c, file = file, results_dir = results_dir, number_format = number_format, width = width_heatmap, height = height_heatmap)
    }
  }
  return(C)
}
