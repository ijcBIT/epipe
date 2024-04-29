#' Apply Filtering to DMRs
#'
#' This function applies filtering to Differentially Methylated Regions (DMRs) and generates plots.
#'
#' @param dmrs Data.table containing DMR information.
#' @param plots Logical, indicating whether to generate plots.
#' @param p.value Vector of p-values.
#' @param mDiff Vector of methylation differences.
#' @param min.cpg Vector of minimum CpG values.
#' @param path Path to save intermediate DMR files.
#'
#' @return A melted data.table containing the filtered DMRs.
#'
#' @author Izar de Villasante
#' @import ggplot2
#'
#' @examples
#' # Example usage:
#' filtered_dmrs <- apply_filter_dmrs(dmrs = dmrs_data, plots = TRUE)
#'

apply_filter_dmrs <- function(dmrs, plots = TRUE, p.value = seq(0.001, 0.11, .025),
                              mDiff = seq(0.15, 0.5, .05), min.cpg = seq(3, 6, 1), path = "analysis/intermediate/dmrs/") {
  if (nrow(dmrs) < 1 | is.null(dmrs)) {

    return(warning("No DMRs available..."))
  } else {
    require(ggplot2)
    dir.create(path)
    params <- expand.grid(p.value, mDiff, min.cpg, s = TRUE)

    res1 <- with(params, Map(function(a, b, c, d) {
      dt <- filter_dmrs(dmrs, a, b, c, d)
      dt$p.val <- a
      dt$mDiff <- b
      dt$min.cpg <- c
      dt
    }, Var1, Var2, Var3, s))
    pdata <- rbindlist(res1)
    colA <- names(pdata)[endsWith(names(pdata), "DMRS")]
    colB <- names(pdata)[endsWith(names(pdata), "Genes")]
    pd <- melt(pdata, measure.vars = list(colA, colB), value.name = c("DMRS", "Genes"),
               variable.name = "type", variable.factor = TRUE)
    levels(pd$type) <- c("Hyper", "Hypo")
    if (plots) {
      plt_list <- list()

      pd2 <- make_ribbon_dt(dt = pd, var = "min.cpg")

      merge(pd, pd2) -> pd3
      (plt_list[["p_Genes"]] <- ggplot2::ggplot(data = pd3, aes(x = mDiff, y = DMRS, group = factor(type), color = factor(type))) +
          geom_line(aes(y = Genes)) +
          theme_minimal() +
          facet_grid(Contrast ~ min.cpg)
      )
      (plt_list[["p_r_Genes"]] <- ggplot2::ggplot(data = pd3, aes(x = mDiff, y = DMRS, group = factor(type), color = factor(type))) +
          geom_line(aes(y = Genes)) +
          geom_ribbon(aes(ymin = min.min.cpg.Genes, ymax = max.min.cpg.Genes,)) +
          facet_grid(Contrast ~ ., margins = FALSE)
      )

      (plt_list[["p_DMRS"]] <- ggplot2::ggplot(data = pd3, aes(x = mDiff, y = DMRS, group = factor(type), color = factor(type))) +
          geom_line() +
          facet_grid(Contrast ~ min.cpg, margins = "min.cpg")
      )
      (plt_list[["p_r_DMRS"]] <- ggplot2::ggplot(data = pd3, aes(x = mDiff, y = DMRS, group = factor(type), color = factor(type))) +
          geom_ribbon(aes(ymin = min.min.cpg.Genes, ymax = max.min.cpg.Genes)) +
          facet_grid(Contrast ~ ., margins = FALSE)
      )

      lapply(1:length(plt_list), function(x)
        save_plot(object = plt_list[[x]],
               filename = paste0(names(plt_list[x]), ".png"),
               path = path
               ))
    }
    return(pd)
  }
}

#' Filter DMRs
#'
#' This function filters Differentially Methylated Regions (DMRs) based on specified criteria.
#'
#' @param dmrs Data.table containing DMR information.
#' @param p.value False Discovery Rate threshold.
#' @param mDiff Methylation difference threshold.
#' @param min.cpg Minimum number of CpGs in a DMR.
#' @param s Logical, indicating whether to generate a summary of filtered DMRs.
#'
#' @return A data.table containing the filtered DMRs.
#'
#' @author Izar de Villasante
#'

filter_dmrs <- function(dmrs, p.value = 0.01, mDiff = 0.2, min.cpg = 5, s = FALSE) {
  if (length(dmrs) < 1 | is.null(dmrs)) {
    warning("No DMRs available...")
    return(dmrs <- data.table::as.data.table(dmrs))
  } else {

    require(data.table)
    dmrs <- data.table::as.data.table(dmrs)
    out <- dmrs[HMFDR <= p.value & abs(meandiff) >= mDiff & no.cpgs >= min.cpg, ]
    if (s) out <- summary_dmrs(out, write = FALSE)
    return(out)
  }
}

#DMRs
#' Summarize DMRs
#'
#' This function summarizes Differentially Methylated Regions (DMRs) by contrast.
#'
#' @param dmrs Data.table containing DMR information.
#' @param path Path to save the summary, defaults to "/results/dmrs/".
#' @param write Logical, indicating whether to write the summary to a file.
#'
#' @return A data.table summarizing DMRs by contrast.
#'
#' @author Izar de Villasante
#'

summary_dmrs <- function(dmrs, path = "/results/dmrs/", write = FALSE) {
  dmrs[, Type := ifelse(meandiff > 0, "Hyper", "Hypo")]
  dmrs.l <- dmrs[, .(Total_DMRS=.N,Hyper_DMRS = sum(Type == "Hyper"), Hypo_DMRS = sum(Type == "Hypo")), by = c("Contrast")]
  genes.l <- dmrs[, .(Total_Genes_DMRS=length(unique(unlist(strsplit(overlapping.genes,',')))),Hyper_Genes_DMRS = length(unique(unlist(strsplit(overlapping.genes[Type == "Hyper"], ",")))),
                      Hypo_Genes_DMRS = length(unique(unlist(strsplit(overlapping.genes[Type == "Hypo"], ","))))), by = c("Contrast")]
  summary <- merge(dmrs.l, genes.l)
  
  # Order the data frame:
  library(dplyr)
  
  summary <- summary %>%
    select(Contrast, matches("^Total"),matches("Hyper"),matches("Hypo"))
  
  
  if (write) data.table::fwrite(summary, path)
  return(summary)
}

# summary_df <- dmrs %>%
#   group_by(Contrast) %>%
#   summarize(
#     Total_DMRs = n(),
#     Total_genes = length(unique(unlist(strsplit(overlapping.genes, ',')))),# Total number of rows per Contrast
#     Hyper_DMRs = sum(Type == 'Hyper'), # Total number of rows with Type 'Hyper'
#     Hyper_genes= length(unique(unlist(strsplit(overlapping.genes[Type == "Hyper"], ",")))),
#     Hypo_DMRs = sum(Type == 'Hypo'), # Total number of rows with Type 'Hypo'
#     Hypo_genes=length(unique(unlist(strsplit(overlapping.genes[Type == "Hypo"], ","))))
#   )
# 





#' Make Ribbon Names
#'
#' This function generates ribbon names based on a variable.
#'
#' @param var Variable for which ribbon names are generated.
#'
#' @return A vector of ribbon names.
#'
#' @author Izar de Villasante
#'
#' @examples
#' # Example usage:
#' ribbon_names <- make_ribbon_names(var = "min.cpg")
#'

make_ribbon_names <- function(var) {
  varnames <- apply(expand.grid(c("min", "max"), var, c("DMRS", "Genes")), 1, function(x) paste(x, collapse = "."))
  return(varnames)
}

#' Make Ribbon Data Table
#'
#' This function generates a data.table for ribbon plotting based on a variable and data.table.
#'
#' @param var Variable for which ribbon data table is generated.
#' @param dt Data.table containing DMR information.
#'
#' @return A data.table for ribbon plotting.
#'
#' @author Izar de Villasante
#'
#' @examples
#' # Example usage:
#' ribbon_dt <- make_ribbon_dt(var = "min.cpg", dt = dt_data)
#'

make_ribbon_dt <- function(var, dt) {
  n <- names(dt)
  sdcols <- setdiff(n, c("Genes", "DMRS", var))
  dt2 <- dt[, make_ribbon_names(var) := .(
    min(DMRS),
    max(DMRS),
    min(Genes),
    max(Genes)
  ), by = sdcols]
  return(dt2)
}
