#' Get Pathways function
#' @title Get Pathways
#' @description Get pathway results based on DMRS
#' @param dmrs Data table containing DMRs.
#' @param res.folder Path to save the plots and intermediate files.
#' @param cols Columns to include in the results.
#' @param pval P-value threshold.
#' @param topN Number of pathways to show on the summary.
#' @param savefile Logical, whether to save the results to a file.
#' @return Pathway results as a data.table.
get_pathways <- function(dmrs, res.folder = "results/pathways/", cols = c("term_size", "query_size", "intersection_size"), pval = 0.05, topN = 50, savefile = FALSE) {
  require(data.table)
  # Convert dmrs to data.table
  data.table::setDT(dmrs)
  varnames = c("Contrast", "FDR", "p_value",cols, "TERM", "method")
  empty_dt <- setNames(data.table(matrix(nrow = 0, ncol = length(varnames))), varnames)
  # Check if there are genes in the input data
  if (nrow(dmrs) < 1) {
    warning("No genes supplied")

    # Create an empty results table with column names

    results <-empty_dt
  } else {
    # Create the results folder if it doesn't exist
    if (savefile == TRUE) suppressWarnings(dir.create(res.folder))

    # Obtain pathway results for all, hyper, and hypo DMRs
    full_pathways <- pathway(dmrs[, .SD, .SDcols = c("Gene", "Contrast")], cols = cols, pval = pval, topN = topN,re.folder=res.folder,enrich_plots=TRUE)
    if(nrow(full_pathways) == 0) full_pathways <- empty_dt
    if (savefile == TRUE) data.table::fwrite(full_pathways, file.path(res.folder, "full_pathway.csv"))

    hyper_pathways <- pathway(dmrs[meandiff > 0, .SD, .SDcols = c("Gene", "Contrast")], cols = cols, pval = pval, topN = topN)
    if(nrow(hyper_pathways) == 0) hyper_pathways <- empty_dt
    if (savefile == TRUE) data.table::fwrite(hyper_pathways, file.path(res.folder, "hyper_pathway.csv"))

    hypo_pathways <- pathway(dmrs[meandiff < 0, .SD, .SDcols = c("Gene", "Contrast")], cols = cols, pval = pval, topN = topN)
    if (savefile == TRUE) data.table::fwrite(hypo_pathways, file.path(res.folder, "hypo_pathway.csv"))
    if(nrow(hypo_pathways) == 0) hypo_pathways <- empty_dt


    # Add a status column to each subset of pathway results
    hypo_pathways$status = "hypo"
    hyper_pathways$status = "hyper"
    full_pathways$status = "both"

    # Combine results for hypo, hyper, and full pathways
    results <- rbind(hypo_pathways, hyper_pathways, full_pathways)
  }

  return(results)
}
