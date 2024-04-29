#' Create Directories for Data Storage in the Pipeline
#'
#' This function creates directories for storing various types of results and plots
#' generated during the pipeline analysis.
#'
#' @param results_folder The main folder to store all results. Default is "./results/".
#' @param analysis_folder The folder to store analysis-specific results. Default is "./analysis/".
#' @param subf A subfolder name to organize results further.
#'
#' @return A list of folder paths created for different purposes.
#'
#' @examples
#' # Create directories with default paths
#' make_results_dirs(subf = "sample_analysis")
#'
#' # Create directories with custom paths
#' make_results_dirs(results_folder = "./my_results/", analysis_folder = "./my_analysis/", subf = "custom_analysis")
#'
#' @export
make_results_dirs <- function(subf, results_folder = "./results/", analysis_folder = "./analysis/") {
  params <- list(
    results_folder = results_folder,
    qc_folder = paste(results_folder, subf, "/QC/", sep = .Platform$file.sep),
    ss_clean_path = paste(analysis_folder, subf, sep = .Platform$file.sep),
    sexplot_folder = paste(results_folder, subf, "plots/sex/", sep = .Platform$file.sep),
    bplots_folder = paste(results_folder, subf, "plots/pca/bplots/", sep = .Platform$file.sep),
    corrplot_folder = paste(results_folder, subf, "plots/pca/corrplot/", sep = .Platform$file.sep),
    dmp_folder = paste(results_folder, subf, "dmps/", sep = .Platform$file.sep),
    dmpplots_folder = paste(results_folder, subf, "dmps/", sep = .Platform$file.sep),
    dmrs_folder = paste(results_folder, subf, "dmrs/", sep = .Platform$file.sep),
    pathway_folder = paste(results_folder, subf, "gopath/", sep = .Platform$file.sep),
    report_folder = paste(results_folder, subf, "reports/", sep = .Platform$file.sep),
    temp = paste(analysis_folder, subf,"temp", sep = .Platform$file.sep),
    heatmap_folder=paste(results_folder, subf, "plots/heatmap/", sep = .Platform$file.sep),
    betas_folder=paste(results_folder, subf, "betas/", sep = .Platform$file.sep)
  )
  sapply(params, function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE))

  return(params)
}
