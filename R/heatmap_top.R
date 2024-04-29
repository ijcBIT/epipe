#' Generate nice HEATMAP for top beta values
#' 
#' This function generates a heatmap of top 1000 variable CpG probes. 
#' 
#' @title Heatmap
#' @param top_beta_1000 A matrix of the top 1000 beta values. 
#' @param metadata A data frame containing sample metadata
#' @param pal Color palette for sample condition(default: Brewer palette "Dark2")
#' @param path Path to save the plots (default: "./")
#' 
#' @import ComplexHeatmap
#' 
#' @return Heatmap
#' 
#' @export
#' 
#' 
#' 
heatmap_top<-function(top_beta_100,ss,sampGroups,pal = RColorBrewer::brewer.pal(8, "Dark2"),
                      path = "./",filename = ""){
  
  library(ComplexHeatmap)
  
  # Add top annotations
  
  Condition <- factor(ss[[sampGroups]])
  unique_conditions <- unique(Condition) 
  color_palette <- pal[1:length(unique_conditions)]
  
  condition_colors<-setNames(color_palette,unique_conditions)
  
  ha <- HeatmapAnnotation(
    df = data.frame(Sample_Group = Condition),
    col = list(Sample_Group = condition_colors)
  )
  
  
  
  # Generate the heatmap
  plot<-Heatmap(top_beta_100, 
          name='methylation',
          row_title = 'CpG probes',
          column_title='Samples',
          show_row_names = F,top_annotation = ha)
  
  # Save the plot
  save_plot(plot, path = path, filename = filename)
  
  return(plot)
  
}



#' 
#' #'
#' 
#' 
#' heatmap_dmps<-()
#' 
#' dmps_f
#' 
#' 
#' 
#' result <- list()  # Create an empty list to store the results
#' 
#' for (i in levels(dmps_f$Contrast)) {
#'   c <- head(dmps_f[dmps_f$Contrast == i, 'ProbeID'], 4)
#'   probe_names<-as.character(c$ProbeID)
#'   top4_betas_DMPs<-betas[which(rownames(betas) %in% probe_names),]
#'   
#' }



