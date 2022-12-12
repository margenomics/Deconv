
##' @name FARDEEP_Deconv
##' @alias FARDEEP_Deconv
##' @title FARDEEP_Deconv
##'
##' @usage Function that generates the relative or absolute deconvolution plots with the FARDEEP method, with the non-fractional plot, and returns a df resulting from the deconvolution.
##' @param matrix Parameter in which a matrix (dataframe) is introduced, which in the first column has the names of the genes and in the rest of the columns the expressions for each gene of the different samples.
##' @param sig.matrix Parameter in which a matrix (dataframe) is introduced, with the first column containing the names of the genes and the rest of the columns containing the expression signatures of different cell types.
##' @param method Parameter where you enter "abs" if you want to do the absolute deconvolution of the matrix, by default it is "rel".
##' @param results_dir Parameter to enter the directory where the graphics will be generated (e.g. "/bicoh/nidia/Deconv").
##' @param height_deconv Parameter in which a value is entered that determines the height of the bar chart, by default it is 10.
##' @param width_deconv Parameter in which a value is entered that determines the width of the bar chart, by default it is 9.
##' @param height_heatmap Parameter in which a value is entered that determines the height of the heatmap, by default it is 700.
##' @param width_heatmap Parameter in which a value is entered that determines the width of the heatmap, by default it is 700.
##' @param name Parameter to enter the name you want to be included in the generated graphics (Ex. "LM22" -> "FARDEEP_Deconv_LM22_rel_plot.png"), default would be ("FARDEEP_Deconv_rel_plot.png").
##' @param number_format Parameter that allows to change the visualization of the numbers in the heatmap. ("\%.2f") for two decimals and ("\%.1e") for exponential notation.
##' @param byCond Parameter in which TRUE is introduced if we introduce a vector so that the function generates the graphs dividing the samples by its condition, by default it is FALSE.
##' @param cond Vector that assigns a condition to each sample of the data frame.
##' @param intercept Parameter that defines that the samples may contain percentages of other cell types not contemplated in the reference matrix, by default it is TRUE.
##' @param permn_number Parameter in which the number of permutations needed to calculate the PValue is noted, default is 10.
##' @param x.size This parameter allows you to select the font size of the X-axis in the bar chart. The default is 10.
##' @param y.size This parameter allows you to select the font size of the Y-axis in the bar chart. The default is 10.
##' @param l.size This parameter allows you to select the font size of the legend in the bar chart. The default is 10.
##' @return Returns the generated graphs and the df of the deconvolution.
##' @author Nidia Barco Armengol
##' @export
##'
##' @examples
##' c<- FARDEEP_Deconv(matrix = matrix, sig.matrix = sig.matrix, results_dir = results_dir, byCond = TRUE, cond = fractions, method = "abs")

FARDEEP_Deconv<- function(matrix, sig.matrix, method="rel", results_dir, height_deconv= 10, width_deconv= 9, height_heatmap= 700, width_heatmap= 700, name= NULL, number_format= "%.2f", byCond= FALSE, cond, data4Tyers= NULL, intercept= TRUE, permn_number= 10, x.size= 10, y.size=10, l.size=10){
  require(usethis)
  require(devtools)
  require(gplots)
  require(data.table)
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(immunedeconv)
  require(FARDEEP)
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
  # For FARDEEP the preprocessing returns a list with the two required df's, so they have to be selected as list items in order to be used.
  M<- Preprocesment_EPIC_FARDEEP(matrix = matrix, sig.matrix = sig.matrix)
  cpm<- M[[1]]
  sig.mtrx<- M[[2]]
  if (method=="rel"){

    # FARDEEP makes the deconvolution over 1.
    tryCatch({
      message("If it appears:")
      message("Error in dimnames(x) <- dn :")
      message("length of 'dimnames' [2] not equal to array extent")
      message("There are samples that do not identify with any cell type in the reference matrix, try intercept= FALSE to disregard the possible presence of other cell types or permn_number= 0, so that samples with no content are still displayed. ")
      RESULTS = t(FARDEEP::fardeep(sig.mtrx, cpm, nn = TRUE, intercept = intercept, permn = permn_number, QN = FALSE)$relative.beta)
    })

    RESULTS <- as.data.frame(RESULTS)
    RESULTS$cell_type=rownames(RESULTS)
    FA <- RESULTS
  }else{

    # FARDEEP makes the deconvolution over 1.
    tryCatch({
      message("If it appears:")
      message("Error in dimnames(x) <- dn :")
      message("length of 'dimnames' [2] not equal to array extent")
      message("There are samples that do not identify with any cell type in the reference matrix, try intercept= FALSE to disregard the possible presence of other cell types or permn_number= 0, so that samples with no content are still displayed. ")
      RESULTS = t(FARDEEP::fardeep(sig.mtrx, cpm, nn = TRUE, intercept = intercept, permn = permn_number, QN = FALSE)$abs.beta)
    })

    RESULTS <- as.data.frame(RESULTS)
    RESULTS$cell_type=rownames(RESULTS)
    FA <- RESULTS
  }
  if (byCond==FALSE){
    # If there is no condition vector.
    if (is.null(name)){
      file <- c("FARDEEP_fractions_",  method , "_","plot.png")
    }else{
      title <- name
      file <- c("FARDEEP_fractions_", title , "_" , method , "_","plot.png")
    }
    file <- paste(file, collapse = "")
    Deconvolution_graph(df= FA, file = file, results_dir = results_dir, height = height_deconv, width = width_deconv, x.size= x.size, y.size= y.size, l.size= l.size)
    if (is.null(name)){
      file <- c("FARDEEP_heatmap_", method , "_","plot.png")
    }else{
      # If there is condition vector.
      title <- name
      file <- c("FARDEEP_heatmap_", title , "_" , method , "_","plot.png")
    }
    file <- paste(file, collapse = "")
    Heatmap_graph(df= FA, file = file, results_dir = results_dir, number_format = number_format, height = height_heatmap, width = width_heatmap)
  }else{
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
      df <- FA$cell_type
      df <- as.data.frame(df)
      df_names <- c("cell_type")
      for (i in e){
        c <- FA[,i]
        df <- cbind(df, c)
        df_names <- c(df_names, colnames(FA)[i])
      }
      colnames(df) <- df_names

      if (is.null(name)){
        file <- c("FARDEEP_fractions_", method,"_", unics[x], "_","plot.png")
      }else{
        title <- name
        file <- c("FARDEEP_fractions_", title , "_", method,"_", unics[x], "_","plot.png")
      }
      file <- paste(file, collapse = "")
      c$cell_type=rownames(c)
      Deconvolution_graph(df= c, file = file, results_dir = results_dir, height = height_deconv, width = width_deconv, x.size= x.size, y.size= y.size, l.size= l.size)
      if (is.null(name)){
        file <- c("FARDEEP_heatmap_", method,"_", unics[x], "_","plot.png")
      }else{
        title <- name
        file <- c("FARDEEP_heatmap_", title , "_", method,"_", unics[x], "_","plot.png")
      }
      file <- paste(file, collapse = "")
      Heatmap_graph(df= c, file = file, results_dir = results_dir, number_format = number_format, height = height_heatmap, width = width_heatmap)
    }
  }
  return(FA)
}

