#' Density plot after normalization
#' 
#' This function generates a density plot of beta-values distribution after normalization.
#' 
#' @param object A MethylSet object with normalized data
#' @param metadata Sample sheet information
#' @param sampGroups Variable of the metadata object for coloring groups
#' @param path path to the folder where plots will be saved. (Default ./)
#' 
#' @return plots
#' 
#' @import Minfi
#' @export


denplot<-function(object, metadata, sampGroups, path="./"){
  library(minfi)
  grDevices::png(file = paste0(path, "density_plot_after_norm.png"))
  Sample_Group <- factor(metadata[[sampGroups]])
  minfi::densityPlot(object, metadata$Sample_Group,main = 'Beta values distribution after normalization')
  dev.off()
}





