#' Get Category from Input Sample Sheet
#'
#' @description Sets the category for each variable in the sample sheet so information in each column is treated properly ("ids", "covs", "batch", "mgroup").
#'
#' @param object data.table or data.frame, sample sheet.
#' @return A vector with each column's category.
#' @export
#' @import data.table
get_category <- function(object) {
  
  library(dplyr)
  # Check if SentrixID is in the dataframe
  
  if (!'Sentrix_ID' %in% names(object)) {
      # If 'Sentrix_ID' is not found and not specified by the user, extract from 'barcode'
      object <- object %>%
        mutate(Sentrix_ID = sub("_.*", "", barcode),
               Sentrix_Position = sub(".*_", "", barcode))
  }
  
  
  # Convert to data.table if not already in that format
  object <- data.table::setDT(object)

  # Define categories for columns
  ids <- c("Sample_Name", "barcode", "Basename","Sentrix_Position")
  batch <- c("Sentrix_ID", "batch")
  covs <- c("Type", "Condition",'Age','predictedSex')
  mgroups <- c("Sample_Group")

  
  category<-NULL
  
  # If category is not defined, assign categories based on column names
  if (is.null(category)) {
    category <- ifelse(colnames(object) %in% ids, "ids",
                       ifelse(colnames(object) %in% batch, "batch",
                              ifelse(colnames(object) %in% covs, "covs",
                                     ifelse(colnames(object) %in% mgroups, "mgroups",
                                            "unknown")
                              )
                       )
    )
  }

  # Perform sanity checks
  if (any(category == "unknown")) {
    warning("Failed to detect category for the following columns: ", paste(colnames(object)[category == "unknown"], collapse = ", "))
  }

  # Check if 'ids' category is present, if not, throw an error
  if (!"ids" %in% category) {
    stop("Sample sheet object must have the category attribute with at least one column as 'ids'.")
  }
  sapply(category,as.character)
  # Add 'category' attribute to the sample sheet
  data.table::setattr(object, "category", category)
  return(object)
}
