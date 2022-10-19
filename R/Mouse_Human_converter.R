
##' @name Mouse_Human_converter
##' @alias Mouse_Human_converter
##' @title Mouse_Human_converter
##'
##' @usage This function convert mouse expression data into human expression data, changing the mouse genes for her analog in human.
##' @param matrix Parameter in which a matrix (dataframe) is introduced, which in the first column has the names of the genes and in the rest of the columns the expressions for each gene of the different samples.
##' @param resultant_file Parameter to enter the name we want the transformed matrix to have (Ex. "counts_mf_human.txt").
##' @param results_dir Parameter to enter the directory where the resultant_matrix.txt will be generated (e.g. "/bicoh/nidia/Deconv").
##' @return This function returns the expression data in human into new file.
##' @author Nidia Barco Armengol
##' @export
##'
##' @examples
##' M<- Mouse_Human_converter(matrix = matrix,resultant_file = "counts_mf_human.txt", results_dir = results_dir)

Mouse_Human_converter<- function(matrix, resultant_file, results_dir){

  require(biomaRt)

  # The gene names are assigned to the rownames.
  matrix<- as.matrix(matrix)
  colnames(matrix)[1]<- "Gene.names"
  rownames(matrix)<- matrix[,1]
  matrix<- matrix[,-1]

  gene.names.mouse <- rownames(matrix)
  matrix <- cbind(matrix,gene.names.mouse)

  # Human and mouse gene databases are taken and compared.
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = gene.names.mouse , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

  # The data is filtered to take only the genes that we have in our data and the genes in our original matrix are changed to their corresponding human genes.
  genesV2.f <- genesV2[genesV2$MGI.symbol %in% gene.names.mouse,]
  gene.names.mouse.f <- gene.names.mouse[gene.names.mouse %in% genesV2.f$MGI.symbol]

  genesV2.f.f <- genesV2.f[genesV2.f$MGI.symbol %in% gene.names.mouse.f,]
  genesV2.f.f <- genesV2.f.f[!duplicated(genesV2.f.f$MGI.symbol),]
  genesV2.f.f<- genesV2.f.f[match(gene.names.mouse.f,genesV2.f.f$MGI.symbol),]

  matrix_human <- matrix
  matrix_human <- matrix_human[rownames(matrix_human) %in% gene.names.mouse.f,]
  rownames(matrix_human) <- genesV2.f.f$HGNC.symbol

  last<- length(colnames(matrix_human))
  matrix_human<- matrix_human[,-last]
  gene.names<- rownames(matrix_human)
  matrix_human<- data.frame(apply(matrix_human, 2, function(x) as.numeric(as.character(x))))
  matrix_human<- cbind.data.frame(gene.names, matrix_human)
  # The resulting matrix is returned as an element and a .txt document is also generated in a selected directory with the transformed matrix.
  file_name<- c(results_dir, "/", resultant_file)
  file_name <- paste(file_name, collapse = "")
  write.table(matrix_human, file= file_name, row.names=FALSE, col.names=TRUE, sep = "\t")
  matrix_human<- as.data.frame(matrix_human)
  return(matrix_human)
}
