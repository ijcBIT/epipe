##############################
# Volcano and manhattan plots
###############################


#' Volcano plot
#'
#' @param object Data frame containing DMPs or DMRs.
#' @path path to the folder where the plots will be saved.

volcanoplot<-function(object,path="./"){

  library(ggplot2)
  library(ggrepel)
  # split the data frame by contrasts
  df_contrasts<-split(object,object$Contrast)

  list_element<-1

  for (df in df_contrasts){

    contrast_name<-names(df_contrasts)[list_element]


    # To plot just the DMP that fall in a GENE region:
    # Subset the DMPs that have a gene associated to
    #df<-df[which(df$UCSC_RefGene_Name!=''),]
    
    
    # Add a column to the data frame to specify if they are UP or DOWN regulated
    df$diffexpression<-'NS' # No significative
    
    
    # IF there are DMPS with an adjusted.pvalue less than 0.05. Use the adj.P.Val as a treshold.
    
    
    if (!is.null(nrow(df[df$adj.P.Val < 0.05, ]))){
      threshold<-'adj.P.Val'
    }else{
      threshold<-'P.Value'
    }
    
    df$diffexpression[df$logFC>0.5 & df[[threshold]]<0.05]<-'UP'
    df$diffexpression[df$logFC< -0.5 & df[[threshold]]<0.05]<-'DOWN'
  
     # Select TOP genes DOWN and UP regulated.
    df_upregulated <- df[df$diffexpression == "UP", ]
    df_upregulated <- df_upregulated[order(df_upregulated[[threshold]], -df_upregulated$logFC), ]
  
    df_downregulated <- df[df$diffexpression == "DOWN", ]
    df_downregulated <- df_downregulated[order(df_downregulated[[threshold]], df_downregulated$logFC), ]
  
    top_upregulated <- head(df_upregulated, 10)
    top_downregulated <- head(df_downregulated, 10)
  
    top_genes <- rbind(top_upregulated, top_downregulated)
  
    colors <- c('DOWN' = '#FF4455', 'NS' = 'black', 'UP' = '#00AFBB')
    colors <- colors[unique(df$diffexpression)]
  
    volcano<-ggplot(data=df, aes(x=logFC, y=-log10(!!rlang::sym(threshold)),col = diffexpression))+
      geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
      geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
      geom_point(size=1)+
      labs(x = "log Fold Change", y = paste0("-log10(",threshold,")"),
            title = "Volcano Plot",color='Expression')+
      ggtitle(paste('Volcano Plot: ',contrast_name))+
      scale_color_manual(values= colors,
                           labels = c('DOWN' = 'Downregulated', 'NS' = 'Not significant', 'UP' = 'Upregulated'))+
      geom_text_repel(data=top_genes,aes(label=ProbeID),max.overlaps = 30,
                       box.padding = 0.35,
                       point.padding = 0.3,show.legend = F)
       

    # Save the plot
    ggplot2::ggsave(filename = paste0(path,"volacanoplot_",contrast_name,".png"),
                    device = "png",plot = volcano,
                    width = 10, height = 12)
    print(paste0(path,"volacanoplot_",contrast_name,".png"))

    list_element <- list_element + 1
  }

}


#volcanoplot(dmps_ex_EPICv2)

# df<-df_contrasts[[2]]

#Add a column to the data frame to specify if they are UP or DOWN regulated
# df_1$diffexpression<-'NS' # No significative
# df_1$diffexpression[df_1$logFC>0.4 & df_1$adj.P.Val<0.01]<-'UP'
# df_1$diffexpression[df_1$logFC< -0.4 & df_1$adj.P.Val<0.01]<-'DOWN'

# Select TOP genes DOWN and UP regulated.
# df_upregulated <- df_1[df_1$diffexpression == "UP", ]
# df_upregulated <- df_upregulated[order(df_upregulated$adj.P.Val, -df_upregulated$logFC), ]
#
# df_downregulated <- df_1[df_1$diffexpression == "DOWN", ]
# df_downregulated <- df_downregulated[order(df_downregulated$adj.P.Val, df_downregulated$logFC), ]
#
# top_upregulated <- head(df_upregulated, 10)
# top_downregulated <- head(df_downregulated, 10)
#
# top_genes <- rbind(top_upregulated, top_downregulated)

#top_genes$UCSC_RefGene_Name <- sub(";.*", "", top_genes$UCSC_RefGene_Name)

# Volcano plot:
# ggplot(data=df_1, aes(x=logFC, y=-log10(P.Value),col = diffexpression))+
#   geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
#   geom_vline(xintercept = c(-0.4, 0.4), col = "gray", linetype = 'dashed') +
#   geom_point(size=0.3)+
#   labs(x = "log Fold Change", y = "-log10(adj.p.value)",
#        title = "Volcano Plot",color='expression')+
#   scale_color_manual(values= c('#FF4455','black','#00AFBB'),
#                      labels= c('Downregulated','Not significant','Upregulated'))+
#   geom_text_repel(data=top_genes,aes(label=ProbeID),
#                   max.overlaps = 30,
#                   box.padding = 0.35,
#                   point.padding = 0.3,show.legend = F)




########################################################################################
########################################################################################
#
# ## MANHATTAN PLOT


#' Manhtattan plot
#'
#' @param object Data frame containing DMPs or DMRs.
#' @path path to the folder where the plots will be saved.

manhattanplot<-function(object,path="./"){

  library(karyoploteR)

  df_contrasts<-split(object,object$Contrast)

  list_element<-1

  for (df in df_contrasts){

    contrast_name<-names(df_contrasts)[list_element]
    
    if (!is.null(nrow(df[df$adj.P.Val < 0.05, ]))){
      threshold<-'adj.P.Val'
    }else{
      threshold<-'P.Value'
    }
    
    png(filename = paste(path,'manhattanplot_',contrast_name,'.png',sep=''),width = 2000,height = 1600,res = 300,units = 'px')
    # Create a Granges object
    gr<-makeGRangesFromDataFrame(df,start.field = 'pos',end.field = 'pos')
    mcols(gr) <- df[, c(threshold, 'B')]
    pval <- mcols(gr)[[threshold]]

    kp <- plotKaryotype(chromosomes = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7","chr8",
                                        "chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                                        "chr17","chr18","chr19","chr20","chr21","chr22"), plot.type = 4,ideogram.plotter = NULL,labels.plotter = NULL)


    kp <- kpPlotManhattan(kp, data=gr,pval=pval,points.col = c('blue4','orange3'),points.cex = 0.5,
                          genomewideline = 0, ymax=8)

    kpAxis(kp, ymin=0, ymax=kp$latest.plot$computed.values$ymax)
    kpAddChromosomeNames(kp, srt=90, cex=1.2, yoffset=20)
    kpAddLabels(kp,labels = paste('-log10(',threshold,')'),srt=90,label.margin = 0.038,pos=2)

    cpgs <- kp$latest.plot$computed.values$data
    suggestive.thr <- kp$latest.plot$computed.values$suggestiveline
    top.cpgs <- tapply(seq_along(cpgs), seqnames(cpgs), function(x) {
      in.chr <- cpgs[x]
      top.snp <- in.chr[which.max(in.chr$y)]
      return(names(top.snp))
    })
    top.cpgs <- top.cpgs[cpgs[top.cpgs]$y>suggestive.thr]
    top.cpgs <- cpgs[top.cpgs]

    if (length(seqnames(top.cpgs))!=0) {

      kpPlotMarkers(kp, data=top.cpgs, labels=names(top.cpgs), srt=45, y=0.9,
                    ymax=8, r0=0.8, line.color="red")
      kpSegments(kp, data=top.cpgs, y0=top.cpgs$y, y1=6.5 ,ymax=8, col="red")

    }
    par(mar = c(5, 5, 4, 2) + 0.1)
    title=paste('Manhattan plot: ',contrast_name)
    kpAddMainTitle(kp, main=title)

    print(paste0(path,"volacanoplot_",contrast_name,".png"))

    dev.off()

    list_element <- list_element + 1
  }

}



#manhattanplot(dmps_ex_EPICv2,path)





###############

# library(qqman)
#
# df_1$chr<-as.numeric(sub('chr','',df_1$chr))
#
# threshold1 <- -log10(0.000005)
#
# manhattan(df_1, chr = "chr", bp = "pos", p = "adj.P.Val",snp='ProbeID',
#           col=c('blue4','orange3'),genomewideline = F,suggestiveline = F,cex=0.5,cex.axis=0.6)
#
# abline(h = threshold1, col = "red", lty = 1)






