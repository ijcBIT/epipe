#' Perform Principal Component Analysis (PCA) on Beta Values
#'
#' This function applies Principal Component Analysis (PCA) to a matrix of beta values.
#'
#' @param beta_top100 A matrix of beta values.
#' @param scale Logical, indicating whether to scale the variables (default: TRUE).
#' @param center Logical, indicating whether to center the variables (default: TRUE).
#' @return A prcomp object representing the results of PCA.
#'
#' @importFrom stats prcomp ggplot2
#'
#' @examples
#' # Example usage:
#' # pca_res(beta_top100, scale = TRUE, center = TRUE)
#'
#' @export
pca_res <- function(beta_top100,ss,sampGroups,path="./", scale = TRUE, center = TRUE,filename = "") {
  library(ggplot2)
  library(ggfortify)
   # Perform PCA on the transposed matrix
  pca_result <- prcomp(t(beta_top100), scale = scale, center = center)
  Sample_Group <- factor(ss[[sampGroups]])
  
  #Get proportion of variance explained bye each pca
  proportion<-summary(pca_result)$importance['Proportion of Variance',]
  
  #PC1 and PC2
  plt<-ggplot(pca_result,aes(x=PC1,y=PC2,color=ss$Sample_Group))+
    geom_point()+
    ggtitle('PC1 and PC2')+
    theme_bw()+ labs(color = 'Sample Group', 
                          x = paste('PC1 (', proportion[1]*100, '%)', sep = ''), 
                          y = paste('PC2 (', proportion[2]*100, '%)', sep = ''))
  save_plot(plt, path = path, filename = paste0(filename,'PC1PC2')) #"pcaplot_PC1PC2")
  
  #Plot PC3 and PC4
  plt34<-ggplot(pca_result,aes(x=PC3,y=PC4,color=ss$Sample_Group))+
    geom_point()+
    ggtitle('PC3 and PC4')+
    theme_bw()+ labs(color = 'Sample Group',
                     x = paste('PC3 (', proportion[3]*100, '%)', sep = ''),
                     y = paste('PC4 (', proportion[4]*100, '%)', sep = ''))
  save_plot(plt34, path = path,filename = paste0(filename,'PC3PC4'))
  
  
  return(pca_result)
}

