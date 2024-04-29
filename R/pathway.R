
#' Perform Pathway Enrichment Analysis for Differentially Methylated Regions (DMRs)
#'
#' This function conducts pathway analysis using the gprofiler2 package for each contrast defined in the "Contrast" column of the input DMRs table.
#' The results for each contrast are combined into a single data.table.
#'
#' @param dmrs Data table containing information about DMRs.
#'             It is assumed to have columns named "Contrast" and "Gene."
#' @param path Path to save the combined pathway analysis results (default: "results/pathways.csv").
#' @param cols Columns to include in the pathway results (default: c("term_size", "query_size", "intersection_size")).
#' @param pval P-value threshold for significance in pathway analysis (default: 0.05).
#' @param topN Minimum number of terms for each group (default: 50).
#' @param savefile Logical indicating whether to save the results to a file (default: FALSE).
#'
#' @return A data.table containing pathway analysis results for all contrasts.
#'
#' @author izar de Villasante and Marina Vilardell
#' @export
#' @import gprofiler2
#' @import data.table
#'
#' @examples
#' # Example usage:
#' dmrs <- data.table(Contrast = c("Group1", "Group1", "Group2", "Group2"),
#'                    Gene = c("Gene1", "Gene2", "Gene3", "Gene4"))
#' result <- pathway(dmrs)
#'
#' @seealso
#' \code{\link{path_results}} for additional details on result formatting.
#'
#' @references
#' For gprofiler2 package documentation, see: https://cran.r-project.org/package=gprofiler2
#'
#' @keywords pathway enrichment DMRs methylation gprofiler2 data.table
#'
#' @examples
#' # Example usage:
#' dmrs <- data.table(Contrast = c("Group1", "Group1", "Group2", "Group2"),
#'                    Gene = c("Gene1", "Gene2", "Gene3", "Gene4"))
#' result <- pathway(dmrs)
#'
#' @seealso
#' \code{\link{path_results}} for additional details on result formatting.
#'
#' @references
#' For gprofiler2 package documentation, see: https://cran.r-project.org/package=gprofiler2
#'
#' @keywords pathway enrichment DMRs methylation gprofiler2 data.table
#'
pathway <- function(dmrs, path = "results/pathways.csv", cols = c("term_size", "query_size", "intersection_size"), pval = 0.05, topN = 50, savefile = FALSE,re.folder ='./',enrich_plots=FALSE) {
  require(gprofiler2)
  require(data.table)
  require(ggplot2)

  # Loop through each unique contrast in the "Contrast" column of dmrs
  pathways <- lapply(unique(dmrs$Contrast), function(cont) {

    # Perform pathway analysis using gprofiler2::gost
    p2 <- gprofiler2::gost(signif = TRUE, unique(dmrs[Contrast == cont, Gene], user_threshold = pval),sources = c('GO:MF','GO:CC','GO:BP','KEGG','REAC','TF','HP'))
    
    p<-p2[[1]]
    # Convert the result to a data.table
    dth <- data.table::as.data.table(p)

    # Prepare the result data.table
    pat <- data.table()

    # Adapt output naming to convention similar to missmethyl
    if (length(dth) > 0) {
      dth[, FDR := p_value]
      dth[, TERM := term_name]
      dth[, source := factor(source)]
      dth[, Contrast := cont]

      # Save results using path_results function
      pat <- path_results(
        pathway = dth,
        topN = topN,
        method = "source",
        path = path,
        cols = cols,
        pval = pval,
        savefile = savefile
      )
    }
    
    if (!is.null(p2) & enrich_plots==TRUE){  
      p2$result$query<-cont
      #Make interactive plots 
      enrichplot_interactive<-gostplot(p2,capped=FALSE,interactive = TRUE)
      htmlwidgets::saveWidget(enrichplot_interactive,file = file.path(re.folder,paste0('enrichplot_interactive_',cont,'.html')))
      
      #Make static plots 
      enrichplot_static<-gostplot(p2,capped=FALSE,interactive = FALSE)
      save_plot(enrichplot_static, path = re.folder, filename = paste0('enrichplot_static_',cont))
      
      #Select top 3 pathways per each category
      library(dplyr)
      top_results<- p2$result %>%
        group_by(source)  %>%
        slice_min(order_by = p_value, n = 3)
      
      #highlight the terms and create a table
      
      
      png(file=file.path(re.folder,paste0('table_plot_',cont,'.png')),width=1000,height = 1000)
      publish_gostplot(enrichplot_static, highlight_terms = top_results$term_id)
      dev.off() 
      
      
      #save_plot(publish_gos, path = re.folder, filename = paste0('highlight_table_plot_',cont))
     
    }
   
    return(pat)
  })

  # Combine results for all contrasts into a single data.table
  result <- do.call("rbind", pathways)
  
  return(result)
}

#' Generate Pathway Results in an Excel Sheet
#'
#' This function filters pathway analysis results for at least the topN pathways for each group and includes all terms with FDR <= FDR. The filtered results are then saved to a CSV file and returned as a data.table.
#'
#' @title Generate Results Excel Sheet for Pathway Analysis
#' @param pathway A data.table with pathway analysis results compatible with gopath function output values.
#' @param topN Numeric value. The minimum number of terms to include for each group.
#' @param method A character or a vector of characters specifying the columns to group by. Default is "method".
#' @param pval FDR threshold to filter out terms (default: 0.05).
#' @param path The path to save the results CSV file (default: "results/pathways.csv").
#' @param cols Columns to include in the pathway results (default: NULL).
#' @param savefile Logical indicating whether to save the results to a file (default: FALSE).
#'
#' @return A data.table containing filtered pathway analysis results.
#'
#' @author izar de Villasante
#' @export
#' @import data.table
#'
#' @examples
#' # Example usage:
#' pathway_results <- data.table(
#'   Contrast = c("Group1", "Group1", "Group2", "Group2"),
#'   FDR = c(0.01, 0.02, 0.03, 0.04),
#'   term_size = c(20, 25, 18, 22),
#'   query_size = c(15, 18, 12, 16),
#'   intersection_size = c(10, 12, 8, 11),
#'   TERM = c("Pathway1", "Pathway2", "Pathway3", "Pathway4"),
#'   method = c("Method1", "Method1", "Method2", "Method2")
#' )
#' result <- path_results(pathway_results)
#'
#' @seealso
#' \code{\link{pathway}} for pathway analysis.
#'
#' @keywords pathway analysis results data.table filtering
#'
#' @return A data.table containing filtered pathway analysis results.
#'
#' @examples
#' # Example usage:
#' pathway_results <- data.table(
#'   Contrast = c("Group1", "Group1", "Group2", "Group2"),
#'   FDR = c(0.01, 0.02, 0.03, 0.04),
#'   term_size = c(20, 25, 18, 22),
#'   query_size = c(15, 18, 12, 16),
#'   intersection_size = c(10, 12, 8, 11),
#'   TERM = c("Pathway1", "Pathway2", "Pathway3", "Pathway4"),
#'   method = c("Method1", "Method1", "Method2", "Method2")
#' )
#' result <- path_results(pathway_results)
#'
#' @seealso
#' \code{\link{pathway}} for pathway analysis.
#'
#' @keywords pathway analysis results data.table filtering
#'
path_results <- function(pathway, topN = 50, method = "method", pval = 0.05, path = "results/pathways.csv", cols = NULL, savefile = FALSE) {
  require(data.table)

  # Add a new column 'method' based on the specified grouping columns
  pathway$method <- pathway[[method]]

  # Order the pathway data.table by grouping columns and FDR
  data.table::setorder(pathway, method, FDR)

  # Identify indices of significant terms (FDR < pval) for each group
  sig_idx <- pathway[,.I[FDR < pval], by = method]$V1

  # Identify indices of topN terms for each group
  head_idx <- pathway[,.I[1:min(..topN, .N)], by = c(method, "Contrast")]$V1

  # Combine significant indices and topN indices using union
  selected_indices <- base::union(sig_idx, head_idx)

  # Subset the pathway data.table using the selected indices
  res <- pathway[selected_indices,]

  # Extract relevant columns for the results
  results <- res[,.SD, .SDcols = c("Contrast", "FDR", "p_value",cols, "TERM", "method")]

  # Order the results by Contrast, method, and FDR
  data.table::setorder(results, Contrast, method, FDR)

  # Save the results to a CSV file if savefile is TRUE
  if (savefile == TRUE) data.table::fwrite(results, path)

  return(results)
}
