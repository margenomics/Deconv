
##' @name GEDIT_functions
##' @alias GEDIT_functions
##' @title GEDIT_functions
##'
##' @usage Document with the necessary functions for deconvolution with the GEDIT method.
##' @author Nidia Barco Armengol
##' @export

# These are the functions needed to make GEDIT, they have been taken from github (they have license that allows it).
A_GEDITDecon = function(MixMat, RefMat){
  args = commandArgs(TRUE)
  outFile = args[3]
  alpha = 0 #as.numeric(args[5])
  lambda = 0 #as.numeric(args[6])
  intercept = FALSE #as.logical(args[7])
  Predictions = C_GLMRegressionLoop(MixMat, RefMat, alpha, lambda, intercept)
  if (dim(Predictions)[2] == 1){
    Predictions[,2] = Predictions[,1]
  }
  Predictions = Predictions[2:dim(Predictions)[1],]
  Adjusted = DD_Adjust(Predictions)
  Adjusted = round(Adjusted[,],4)
  Adjusted = t(Adjusted)
  return(Adjusted)}

C_GLMRegressionLoop = function(Samples, Refmat, Alph, Lamb, Inter)
{
  Predictions = as.data.frame(matrix(0, ncol = dim(Samples)[2], nrow = 1+dim(Refmat)[2]))
  for (i in 1:dim(Samples)[2]){
    OutMat = as.matrix(coef(glmnet(as.matrix(Refmat), as.matrix(Samples[,i]),
                                   lower.limits = 0.0, alpha = Alph, lambda = Lamb, intercept = Inter)))
    LastColumn = dim(OutMat)[2]
    Predictions[,i] = OutMat[,LastColumn]
  }
  names(Predictions) = names(Samples)
  row.names(Predictions) = c("intercept",names(Refmat))
  return(Predictions)
}

D_AdjustLoop = function(Predictions){
  Adjusted = Predictions
  for (i in 1:length(Predictions)){
    Adjusted[[i]] = DD_Adjust(Predictions[[i]])}
  return(Adjusted)}

DD_Adjust = function(Predictions){
  Adjusted = Predictions
  numPreds = dim(Predictions)[2]
  numCTs   = dim(Predictions)[1]
  for (i in 1:numPreds){
    Adjusted[i] = round(Predictions[i]/sum(Predictions[1:numCTs,i]),5)
    if (sum(Predictions[1:numCTs,i]) == 0.0){
      Adjusted[i] = round(1.0/numCTs,5)
    }
  }
  return(Adjusted)
}

