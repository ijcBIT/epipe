#' Plot Differentially Methylated Probes (DMPs)
#'
#' This function generates various plots to visualize different aspects of differentially methylated probes (DMPs).
#'
#' @title Plot DMPs
#' @param DMPann A data.table containing DMP annotations.
#' @param names A character vector specifying the names for saving the plots.
#' @param path Path to save the plots. If NULL, plots won't be saved.
#'
#' @return A list of ggplot2 objects representing the generated plots.
#'
#' @author Original version by unknown author, modifications and documentation by ChatGPT.
#' @export
#' @import ggplot2
#' @import data.table
#'
#' @examples
#' # Example usage:
#' DMP_annotations <- data.table(
#'   Contrast = c("Group1", "Group1", "Group2", "Group2"),
#'   Type = c("Hyper", "Hyper", "Hypo", "Hypo"),
#'   Relation_to_Island = c("Island", "Shore", "Island", "Shore"),
#'   UCSC_RefGene_Group = c("Gene1;Gene2", "Gene3", "Gene4", ""),
#'   UCSC_RefGene_Group_short = c("Gene1", "Gene3", "Gene4", ".")
#' )
#' plot_list <- plotDMP(DMPann = DMP_annotations, names = "DMP_plots", path = "analysis/DNA_methylation/")
#'
plotDMP <- function(DMPann, names= c( "DMP_count.png","DMP_count_facet.png","DMP_annCGI.png", "DMP_annGenomic.png"), path = NULL) {
  if (nrow(DMPann) > 0) {
    library(ggplot2)
    library(data.table)
    data.table::setDT(DMPann)

    # Plot DMPs (hypo/hyper)
    g1 <- ggplot2::ggplot(DMPann, aes(Contrast, fill = Type)) +
      geom_bar(position = "dodge", stat = "count") +
      theme_bw() +
      scale_fill_manual(values = c("red", "skyblue")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "", y = "Count", fill = 'Methylation') +
      ggtitle('Differentially Methylated Probes')

    # Plot with facets
    g2 <- ggplot2::ggplot(DMPann, aes(Contrast, fill = Type)) +
      geom_bar(position = "dodge", stat = "count") +
      facet_wrap(. ~ Type, scales = "free_x") +
      theme_bw() +
      scale_fill_manual(values = c("red", "skyblue")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "", y = "Count", fill = 'Methylation') +
      ggtitle('Differentially Methylated Probes')

    # Plot proportion of DMPs in CGI
    DMP_annCGI <- data.frame(DMPann[, .(Contrast, Type, Relation_to_Island)])
    g3 <- ggplot2::ggplot(DMP_annCGI, aes(Contrast, fill = Relation_to_Island)) +
      facet_wrap(. ~ Type, scales = "free_x") +
      geom_bar(position = "fill", width = 0.8) +
      theme_bw() +
      scale_fill_brewer(palette = "Set1") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
      ylab("DMPs") +
      xlab("")

    # Plot proportion of DMPs in genomic elements
    DMPann$UCSC_RefGene_Group[which(DMPann$UCSC_RefGene_Group == "")] <- "."
    DMP_annGenomic <- DMPann[, .(UCSC_RefGene_Group_short = unlist(lapply(strsplit(UCSC_RefGene_Group, ";"), '['))),
                             by = c("Contrast", "Type")]
    
    DMP_annGenomic$USCS_refgene_Group_summary<-ifelse(startsWith(DMP_annGenomic$UCSC_RefGene_Group_short,'exon'),'EXON',DMP_annGenomic$UCSC_RefGene_Group_short)
    #DMP_annGenomic$USCS_refgene_Group_summary<-ifelse(DMP_annGenomic$USCS_refgene_Group_summary=='.','NO INFORMATION',DMP_annGenomic$USCS_refgene_Group_summary)
    
    g4 <- ggplot2::ggplot(DMP_annGenomic, aes(Contrast, fill = USCS_refgene_Group_summary)) +
      facet_wrap(. ~ Type, scales = "free_x") +
      geom_bar(position = "fill", width = 0.8) +
      theme_bw() +
      scale_fill_brewer(palette = "Set1") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
      labs(fill = "RefGene") +
      ylab("DMPs") +
      xlab("")

    plt_list <- list(g1, g2, g3, g4)

    if (!is.null(path)) {
      sapply(1:length(plt_list), function(x) {
        save_plot(
          filename = paste0(names, "_", x, ".png"),
          object = plt_list[[x]],
          path = path
        )
      })
    }
    return(plt_list)
  } else {
    warning("Empty data.frame. Please try again modifying filter params.")
  }
}
