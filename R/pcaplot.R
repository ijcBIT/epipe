#' PCA plot
#'
#' 
#' @return pcaplot
#' @import ggplot

pcaplot<-function(pca,ss,sampGroups, path = "./"){
  Sample_Group <- factor(ss[[sampGroups]])
  plt<-ggplot(pca,aes(x=PC1,y=PC2,color=ss$Sample_Group))+geom_point()+ggtitle('PCA1 and PC2')+theme_classic()+labs(color='Sample Groups')
  save_plot(plt, path = path, filename = "pcaplot")
  return(Sample_Group)
}


#library(ggfortify)
#plt<-autoplot(pca,data=ss,colour='Sample_Group',size=2)+theme_classic()+ggtitle('PCA1 and PC2')





#pplot<-pcaplot(pca,
 #             ss,
  #            sampGroups = names(ss)[attributes(ss)$category == "mgroups"][1])
              #path=custom_paths[[bplots_folder]])
              



