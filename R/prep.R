#' Preprocess Methylation Data
#'
#' This function performs preprocessing steps on a MethylSet object.
#'
#' @param mSetSqn MethylSet object to be preprocessed.
#' @param remove_sex Logical, whether to remove sex chromosomes. Default is TRUE.
#' @param sexchr Vector of chromosome names representing sex chromosomes (default is c("chrX", "chrY")).
#' @param arraytype Type of array used for preprocessing. If provided, array-specific annotation is used.
#' @return Preprocessed MethylSet object.
#'
#' @import minfi
#' @import maxprobes
#' @export
prep <- function(mSetSqn, remove_sex = TRUE, sexchr = c("chrX", "chrY"), arraytype = NULL,sexplot_folder= NULL) {
  # Save the initial set of probe IDs
  probeID_start <- rownames(mSetSqn)
  
  # Step 0: Set array-specific annotation if not provided
  if (length(annotation(mSetSqn)) < 2) {
    # Determine array type if not provided
    if (is.null(arraytype)) {
      arraytype <- guess_arraytype(nrow(mSetSqn))
    } else {
      arraytype <- switch(
        arraytype,
        "EPICv2" = "IlluminaHumanMethylationEPICv2",
        "EPIC" = "IlluminaHumanMethylationEPIC",
        "450K" = "IlluminaHumanMethylation450k",
        "27K" = "IlluminaHumanMethylation27k",
        "mammal" = "HorvathMammalMethylChip40",
        arraytype
      )
    }
    # Set array-specific annotation
    mSetSqn <- set_array_annotation(mSetSqn, arraytype)
  }

  # Step 1: Remove probes with known SNPs at CpG site
  mSetSqn <- minfi::mapToGenome(mSetSqn)
  mSetSqn <- minfi::dropLociWithSnps(mSetSqn)
  Snps_rm <- rownames(mSetSqn)
  metadata(mSetSqn)$removed_Snps <- setdiff(probeID_start, Snps_rm)

  # Step 2: Remove cross-reactive probes
  remove_cross_reactive_probes(mSetSqn, sexchr)

  # Step 3: Sex prediction & removal
  mSetSqn$predictedSex <- minfi::getSex(mSetSqn, cutoff = -2)$predictedSex
  grDevices::png(file = paste0(sexplot_folder,"sex_estimation.png"))
  plotSex(addSex(minfi::mapToGenome(mSetSqn)))
  grDevices::dev.off()
  
  
  
  # Check if sex chromosomes need to be removed
  if (remove_sex) {
    mSetSqn <- remove_sex_chromosomes(mSetSqn, sexchr)
  }

  # Save the final set of probe IDs
  probeID_end <- rownames(mSetSqn)

  # Capture the probes that were removed during preprocessing and add to metadata
  metadata(mSetSqn)$preprocess <- setdiff(probeID_start, probeID_end)

  # Return preprocessed MethylSet object
  return(mSetSqn)
}


# Helper function to remove cross-reactive probes
remove_cross_reactive_probes <- function(mSetSqn, sexchr) {
  probeID_start <- rownames(mSetSqn)
  if (annotation(mSetSqn)[1] == "IlluminaHumanMethylationEPICv2") {
    # Read the cross-reactive probe mask for EPICv2
    mask <- qs::qread("data/EPICv2_xreactive.qs")
    mSetSqn <- mSetSqn[!rownames(mSetSqn) %in% mask, ]
  } else {
    # For other arrays, use maxprobes::dropXreactiveLoci
    mSetSqn <- maxprobes::dropXreactiveLoci(mSetSqn)
  }
  metadata(mSetSqn)$removed_sex <- setdiff(probeID_start,rownames(mSetSqn))
  return(mSetSqn)
}

# Helper function to remove sex chromosomes
remove_sex_chromosomes <- function(mSetSqn, sexchr) {
  probeID_start <- rownames(mSetSqn)
  # Get Annotation
  anno <- minfi::getAnnotation(mSetSqn)

  # Remove sex chromosomes
  anno <- anno[!(anno$chr %in% sexchr),]
  mSetSqn <- mSetSqn[rownames(anno),]
  sex_rm <- rownames(mSetSqn)
  metadata(mSetSqn)$removed_sex <- setdiff(probeID_start,rownames(mSetSqn))
  return(mSetSqn)
}
