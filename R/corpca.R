#' Perform Correlation-based PCA and Generate Plot
#'
#' This function performs a correlation-based PCA on beta values and generates an eigencorplot.
#'
#' @param beta_top100 A matrix of beta values.
#' @param metadata A data frame containing sample metadata.
#' @param vars Variables in the metadata to include in the eigencorplot.
#' @param idcol Column name in metadata containing unique identifiers for samples (default: "barcode").
#' @param path Directory path for saving the plot (default: "./").
#' @param filename Name of the file to save the plot (default: "").
#' @param title Title for the eigencorplot (default: 'PC1-6 clinical correlations').
#' @return A PCAtools eigencorplot object.
#'
#' @importFrom PCAtools pca eigencorplot getComponents save_plot
#'
#' @export
corpca <- function(beta_top100, metadata, vars = NULL, idcol = "barcode",
                   path = "./", filename = "", title = 'PC1-6 clinical correlations') {
  requireNamespace("PCAtools")
  
  #Check if the columns are present in the data frame and remove them if they are.
  if ("Basename" %in% colnames(metadata)) {
    metadata <- subset(metadata, select = -Basename)
  }
  if ("barcode" %in% colnames(metadata)) {
    metadata <- subset(metadata, select = -barcode)
  }
  if ("Sample_Name" %in% colnames(metadata)) {
    metadata <- subset(metadata, select = -Sample_Name)
  }
  if ("Sentrix_Position" %in% colnames(metadata)) {
    metadata <- subset(metadata, select = -Sentrix_Position)
  }

  
  # Convert metadata to data frame
  metadata <- data.frame(metadata, stringsAsFactors = TRUE)
  rownames(metadata) <- metadata[[idcol]]
  colnames(beta_top100) <- rownames(metadata)

  # Perform PCA with PCAtools package
  p <- PCAtools::pca(beta_top100, metadata = metadata, removeVar = 0.1)

  # If vars parameter is not provided, select metadata variables with more than one unique value
  if (is.null(vars)) {
    vars <- names(p$metadata)[sapply(p$metadata, function(x) !any(is.na(x)) & length(unique(x)) > 1)]
  }

  # Filter vars based on unique values
  vars <- vars[sapply(p$metadata[, vars], function(x) length(unique(x)) > 1)]
  
  #par(cex.main = 0.8)  
  
  # Generate eigencorplot
  pcaplt <- PCAtools::eigencorplot(
    p,
    components = PCAtools::getComponents(p, 1:6),
    metavars = vars,
    col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
    cexCorval = 0.7,
    colCorval = 'white',
    fontCorval = 2,
    posLab = 'bottomleft',
    rotLabX = 45,
    posColKey = 'top',
    cexLabColKey = 1.5,
    scale = TRUE,
    main = title,
    colFrame = 'white',
    plotRsquared = FALSE,
    cexMain = 0.8
  )

  # Save the plot
  save_plot(pcaplt, path = path, filename = filename)

  return(pcaplt)
}
