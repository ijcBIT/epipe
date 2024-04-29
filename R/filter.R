#' Filter RGChannelSet based on Quality Control Criteria
#'
#' This function filters an RGChannelSet based on quality control criteria, including mean detection p-values.
#'
#' @param rgSet RGChannelSet object to be filtered.
#' @param sg Column name indicating sample groups.
#' @param sampNames Column name in 'rgSet@colData' containing sample names.
#' @param save_barplot Logical, whether to save the barplot of mean detection p-values.
#' @param frac Fraction of probes not passing the p-value filter for sample exclusion.
#' @param pval P-value threshold for sample and probe filtering.
#' @param remove_sex Logical, whether to remove sex chromosomes.
#' @param arraytype Array type, used to update annotation for Illumina arrays.
#' @param cols Vector of colors for sample groups. If NULL, colors are generated using get_cols.
#' @param qc_folder Folder path for saving quality control plots.
#'
#' @return Filtered RGChannelSet.
#'
#' @importFrom data.table setdiff
#' @importFrom minfi detectionP
#' @importFrom grDevices jpeg
#'
#' @examples
#' # Filter RGChannelSet based on quality control criteria
#' filter(rgSet, sg = "Sample_Group", sampNames = "Sample_Name")
#'
#' @export
filter <- function(rgSet, sg = NULL, sampNames = "Sample_Name", save_barplot = TRUE,
                   frac = 0.1, pval = 0.01, remove_sex = FALSE, arraytype = NULL, cols = NULL,
                   qc_folder = "analysis/intermediate/QC") {
  metadata(rgSet)$bad_samples <- character(0)
  metadata(rgSet)$bad_probes <- character(0)
  targets <- data.table::as.data.table(rgSet@colData, keep.rownames = "rn")
  # Check if 'sampNames' is in 'rgSet@colData', otherwise use colnames(rgSet)
  if (!(sampNames %in% names(targets))){
    sampNames <- colnames(rgSet)
  }else{
    sampNames <- targets[[sampNames]]
  }

  n <- ncol(rgSet)

  # Create 'qc_folder' if it doesn't exist
  dir.create(qc_folder, recursive = TRUE, showWarnings = FALSE)

  if (length(sg) > 1) sg <- sg[1]
  if (is.null(sg) | !(sg %in% names(targets))) {
    SampGroups <- rep(1, ncol(rgSet))
  } else {sampGroups <- as.factor(targets[[sg]])}

  names(sampGroups) <- sampNames

  # Generate colors for sample groups
  if (is.null(cols)) cols <- get_cols(sampGroups)

  # 1. Check quality of the combined probe signal (testing against negative controls)
  detP <- minfi::detectionP(rgSet, type = "m+u")
  if(save_barplot == TRUE){
    grDevices::jpeg(file = file.path(qc_folder, "mean_detection_pvalues.jpeg"), width = 960, height = 1240)
    ylabels <- colnames(detP)
    par(mar = c(max(4.1, max(nchar(ylabels)) / 2.2), 4.1, 4.1, 2.1))
    barplot(colMeans(detP), col = cols, las = 2,
            cex.names = 0.8, ylim = c(0, max(0.002, max(colMeans(detP)) * 2)), main = "Mean detection p-values")
    if(max(0.002, max(colMeans(detP)) * 2)>0.002){
      graphics::abline(h = 0.01, col = "red")
    }else{
      graphics::abline(h = 0.001, col = "red")
    }
    graphics::legend("topleft", legend = levels(sampGroups), fill = levels(factor(cols)),
                     bg = "white")
    grDevices::dev.off()
  }

  # 2. Remove low-quality samples (with fraction of probes not passing pval)
  bad_samples <- colnames(detP)[colSums(detP >= pval) / nrow(detP) > frac]
  if (length(bad_samples) > 0) {
    warning("The following samples will be discarded since they fail to pass the p-value filter ( ",
            frac * 100, "% of the probes with p-val >", pval, "): \n ", paste(bad_samples, collapse = ", "))
    rgSet <- rgSet[, setdiff(colnames(detP), bad_samples)]
  } else {
    cat("All samples passed detection P-value filter")
  }

  # 3. Remove low-quality probes (with p-value below pval)
  bad_probes <- which(rowSums(detP < pval) < ncol(rgSet) * (1 - frac))
  rgSet <- rgSet[-c(bad_probes), ]
  if (length(bad_probes) > 0) {
    warning("The following probes will be discarded since more than", frac * 100,
            "% of the samples have detection p-values > ", pval, "): \n ", paste(bad_probes, collapse = ", "))
  } else {
    cat("All samples passed detection P-value filter")
  }
  metadata(rgSet)$bad_samples <- bad_samples
  metadata(rgSet)$bad_probes <- bad_probes

  return(rgSet)
}

