
##' @name Preprocesment_EPIC_FARDEEP
##' @alias Preprocesment_EPIC_FARDEEP
##' @title Preprocesment_EPIC_FARDEEP
##'
##' @usage Function for realize the pre-procesment for EPIC and FARDEEP methods
##' @param matrix Parameter in which a matrix (dataframe) is introduced, which in the first column has the names of the genes and in the rest of the columns the expressions for each gene of the different samples.
##' @param sig.matrix Parameter in which a matrix (dataframe) is introduced, with the first column containing the names of the genes and the rest of the columns containing the expression signatures of different cell types.
##' @return List with processed matrices.
##' @author Nidia Barco Armengol
##' @export
##'
##' @examples
##' M<- Preprocesment_EPIC_FARDEEP(matrix = matrix, sig.matrix = sig.matrix)

Preprocesment_EPIC_FARDEEP<- function(matrix, sig.matrix= NULL){
  # Removes duplicate genes, and shape the mixed matrix to be used for deconvolution.
  # In case of finding duplicate genes it mean the rows.
  colnames(matrix)[1]<- "Gene.names"
  matrix$Gene.names<- as.character(matrix$Gene.names)
  matrix <- matrix %>%
    group_by(Gene.names) %>%
    summarise_if(is.numeric,mean) %>%
    ungroup()
  names <- matrix$Gene.names
  matrix<- as.matrix(matrix[,-1])
  rownames(matrix)=names

  if (is.null(sig.matrix)){
    # Sig.matrix is null only when used for QuanTIseq, because this method does not allow external sig.matrix.
    return(matrix)
  }else{
    # Check that there are no duplicate genes in the reference matrix, and shape the matrix for deconvolution.
    colnames(sig.matrix)[1]<- "Gene.names"
    sig.matrix$Gene.names<- as.character(sig.matrix$Gene.names)
    sig.matrix <- sig.matrix %>%
      group_by(Gene.names) %>%
      summarise_if(is.numeric,mean) %>%
      ungroup()
    names <- sig.matrix$Gene.names
    sig.matrix=as.matrix(sig.matrix[,-1])
    rownames(sig.matrix)=names

    l<- list(matrix, sig.matrix)
    return(l)
  }
}
