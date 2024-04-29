#' Apply Filtered DMPs and Generate Plots
#'
#' @param dmps Data table containing DMPs.
#' @param dev Device to use for plotting (e.g., "png").
#' @param p.value Vector of p-values to explore.
#' @param mDiff Vector of mean methylation differences to explore.
#' @param path Path to save the plots and intermediate files.
#'
#' @return None (Plots are saved to the specified path).
apply_filter_dmps <- function(dmps, dev = "png", p.value = seq(0.00, 0.1, .01),
                              mDiff = seq(0.15, 0.5, .05), path = "analysis/intermediate/dmps/") {
  if (length(dmps) < 1 | is.null(dmps)) {
    warning("No DMPs available...")
  } else {
    require(ggplot2)

    # Create the directory if it doesn't exist
    dir.create(path)

    # Generate combinations of p-values and mean methylation differences
    params <- expand.grid(p.value, mDiff, T)

    # Define a function to filter DMPs and add parameters
    res1 <- with(params, Map(function(a, b, c) {
      dt <- filter_dmps(dmps, a, b, c)
      dt$p.val <- a
      dt$mDiff <- b
      dt
    }, Var1, Var2, Var3))

    # Combine the results into a data.table
    pdata <- rbindlist(res1)
    pdata[, All := Hypo + Hyper]

    # Reshape the data for plotting
    pd <- melt(pdata, measure.vars = c("Hyper", "Hypo", "All"))

    pd$variable <- factor(pd$variable)

    # Generate plots for each variable
    if (NROW(pd) > 0 & length(levels(pd$variable)) > 1) {
      plt_list <- list()
      plt_list[["dmp_params"]] <- ggplot2::ggplot(data = pd, aes(x = mDiff, y = value, group = p.val, color = p.val)) +
        geom_line(aes(linetype = factor(p.val))) +
        facet_grid(variable ~ Contrast)

      # Save plots to the specified path using custom save_plot function
      lapply(1:length(plt_list), function(x)
        save_plot(object = plt_list[[x]],
                  filename = paste0(names(plt_list[x]), ".", dev),
                  path = path))
    }
  }
}

#' Filter DMPs Function
#'
#' @param dmps Data table containing DMPs.
#' @param p.value P-value threshold.
#' @param mDiff Mean methylation difference threshold.
#' @param s Logical indicating whether to include summary statistics.
#'
#' @return Filtered DMPs.
filter_dmps <- function(dmps, p.value = 0.01, mDiff = 0.3, s = F) {
  require(data.table)
  dmps <- data.table::as.data.table(dmps)
  dmps_f <- dmps[adj.P.Val <= p.value &
                   abs(diff_meanMeth) >= mDiff, ]
  dmps_s <- summary_dmps(dmps_f)
  if(nrow(dmps_f)==0)warning("no DMPs detected with p.value = ",p.value," & mDiff = ",mDiff)
  if (s == T) {
    out <- dmps_s #[, .SD, .SDcols = c("Hypo", "Hyper", "Contrast")]
  } else {
    out <- dmps_f
  }
  return(out)
}

#' Summarize DMPs and Write Results to File
#'
#' @param DMPextr_object Object containing DMPs.
#' @param dir Directory to save the results.
#' @param name Name for the summary (used in file naming).
#' @param write Logical indicating whether to write results to file (default is FALSE).
#'
#' @return Summary of DMPs.
summary_dmps <- function(DMPextr_object, dir = "./results/dmps/", name = "raw", write = FALSE) {
  require(data.table)

  # Convert DMPextr_object to a data.table
  dt <- data.table::as.data.table(DMPextr_object)
 
  
  # Calculate summary statistics for Hyper and Hypo DMPs by Contrast
  dt_summary <- dt[, .(Total_DMPs = .N, Hyper_DMPs=sum(Type == "Hyper")), by = c("Contrast")]
  
  
  dmrs.l <- dt[, .(Total_DMPS=.N,Hyper_DMPS = sum(Type == "Hyper"), Hypo_DMPS = sum(Type == "Hypo")), by = c("Contrast")]
  genes.l <- dt[, .(Total_Genes_DMPS=length(unique(unlist(strsplit(UCSC_RefGene_Name,';')))),Hyper_Genes_DMPS = length(unique(unlist(strsplit(UCSC_RefGene_Name[Type == "Hyper"], ";")))),
                      Hypo_Genes_DMPS = length(unique(unlist(strsplit(UCSC_RefGene_Name[Type == "Hypo"], ";"))))), by = c("Contrast")]
  summary <- merge(dmrs.l, genes.l)
  
  #Order the data frame:
  library(dplyr)

  dt_summary <- summary %>%
    select(Contrast, matches("^Total"),matches("Hyper"),matches("Hypo"))
  
  
  # Write raw DMPs and summary statistics to file if requested
  if (write) {
    data.table::fwrite(dt, paste0(dir, name, "_dmp_raw.csv.gz"))
    data.table::fwrite(dt_summary, paste0(dir, name, "_dmp_summary.txt"))
  }

  return(dt_summary)
}

