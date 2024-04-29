#' Density plot after normalization
#' @return dplot
#' @import Minfi
#' @export


denplot<-function(object, metadata, sampGroups, path="./"){
  library(minfi)
  grDevices::png(file = paste0(path, "density_plot_after_norm.png"))
  Sample_Group <- factor(metadata[[sampGroups]])
  minfi::densityPlot(object, metadata$Sample_Group,main = 'Beta values distribution after normalization')
  dev.off()
}





