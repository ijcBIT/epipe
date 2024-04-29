#' Generate qc plots
#' @title generate qc plots for signal distribution prior to filtering
#' @param rgSet rgSet object containing channel intenisty values
#' @param sampNames variable containing barcodes or ids matching colnames of the rgsetdata
#' @param pal Color palette to use
#' @param sampGroups variables to use for coloring groups
#' @param qc_folder path to the folder where plots will be saved
#' @return plots
#' @author izar de Villasante
#' @export
qc<-function(rgSet,sampGroups=NULL, sampNames= "Sample_Name", pal=NULL,
             qc_folder="analysis/intermediate/QC/",idcol="barcode"){
  require(S4Vectors)
  require(Biostrings)
  require(Biobase)
  require(minfi)
  path=paste0(qc_folder,"Report.pdf")
  colData = rgSet@colData
  ss <- data.table::as.data.table(colData)
  Sample_Group = factor(colData[[sampGroups]])
  if (!(sampNames %in% names(rgSet@colData))) sampNames <- colnames(rgSet)

  unlist(ss[order(Sample_Group),..idcol])->idx
  rgSet<-rgSet[,idx]
  colData = rgSet@colData
  Sample_Group = factor(colData[[sampGroups]])
  
  Sample_Group<-setNames(Sample_Group,idx)
  
  cols<-get_cols(factor_variable = Sample_Group, pal = pal)
  dir.create(qc_folder,recursive=T,showWarnings = F)
  minfi::qcReport(rgSet = rgSet,
                  pdf = paste0(qc_folder,"Report.pdf"),
                  sampGroups= rgSet@colData[[sampGroups]])
                  #sampNames = rgSet@colData[[sampNames]])
  
  #### Density plot
  
  grDevices::png(file = paste0(qc_folder, "density_plot_before_norm.png"))
  minfi::densityPlot(rgSet, 
                     sampGroups = rgSet@colData[[sampGroups]],main = 'Beta values distribution before normalization')
  grDevices::dev.off()
  
  # dp <- densPlot(
  #   rgSet, sg = Sample_Group,main = "Beta",
  #   pal =pal#[!duplicated(cols)] 
  # )
  # save_plot(dp, path=qc_folder, filename = "density_plot")

  # ggplot2::ggsave(filename = paste0(qc_folder,"density_plot.png"),
  #                 device = "png",plot = dp,
  #                 width = 10, height = 12)
  
  
  #### Bean plot
  
  # grDevices::png(file = paste0(qc_folder,"bean_plot.png"))
  # minfi::densityBeanPlot(
  #   rgSet, sampGroups = rgSet@colData[[sampGroups]],
  #   pal=cols[!duplicated(cols)]
  # ) 
  # grDevices::dev.off()

  #save_plot(beanplot, path=qc_folder, filename = "bean_plot")
  
  
  #### Mean QC plot
  
  mSet <- minfi::preprocessRaw(rgSet)
  qc   <- minfi::getQC(mSet)
  grDevices::png(file = paste0(qc_folder,"mean_qc.png"))
  minfi::plotQC(qc)
  grDevices::dev.off()
  #save_plot(mean_qc, path=qc_folder, filename = "mean_qc")
  
}

#' Density Plot for Beta Values (ggplot2)
#'
#' Generates a density plot for beta values from an 'RGChannelSet', 'MethylSet', or matrix using ggplot2.
#'
#' @param dat Data object. Must be an 'RGChannelSet', a 'MethylSet', or a matrix.
#' @param sg Sample groups. If NULL, defaults to one group for all samples.
#' @param main Main title for the plot.
#' @param xlab Label for the x-axis.
#' @param pal Color palette for sample groups (default: Brewer palette "Dark2").
#' @param cols Vector of colors for sample groups. If NULL, colors are generated using get_cols.
#' @param xlim Range for the x-axis.
#' @param ylim Range for the y-axis.
#' @param add Logical, whether to add to an existing plot.
#' @param legend Logical, whether to include a legend.
#' @return A ggplot object.
#'
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom methods is
#' @importFrom minfi getBeta
#' @import data.table
#' @export
densPlot <- function(dat, sg = NULL, main = "", xlab = "Beta",
                     pal = RColorBrewer::brewer.pal(8, "Dark2"), cols = NULL,
                     xlim, ylim, add = TRUE, legend = TRUE, ...) {
  # Load necessary libraries
  library(ggplot2)
  library(data.table)

  # Check the type of 'dat' and extract beta values accordingly
  if (is(dat, "RGChannelSet") || is(dat, "MethylSet")) {
    b <- getBeta(dat)
  } else if (is(dat, "matrix")) {
    b <- dat
  } else {
    stop("argument 'dat' must be an 'RGChannelSet', a 'MethylSet', or matrix.")
  }

  # Handle sample groups
  if (is.null(sg)) {
    sg <- rep(1, ncol(b))
  } else if (length(sg) == 1) {
    sg <- rep(sg, ncol(b))
  }

  # Convert sample groups to factor with proper names
  sampGroups <- as.factor(sg)
  names(sampGroups) <- colnames(b)

  # Generate colors for sample groups
  if (is.null(cols)) cols <- get_cols(sampGroups, pal = pal)
  cols <- factor(cols)

  # Use melt from data.table package to reshape data for ggplot2
  ggplot_data <- melt(as.data.table(b, keep.rownames = "cg"), id.vars = "cg", variable.name = "sample", value.name = xlab)
  ggplot_data[, grp := sampGroups[ggplot_data$sample]]

  # Create ggplot object
  p <- ggplot(ggplot_data, aes(x = .data[[xlab]], color = grp, group = sample)) +
    geom_density() +
    labs(title = main, x = xlab, y = "Density") +
    scale_color_manual(values = levels(cols)) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12))

  # Return the ggplot object
  return(p)
}



# densPlot <- function(dat, sampGroups = NULL, main = "", xlab = "Beta",
#                      pal = RColorBrewer::brewer.pal(8, "Dark2"),cols=NULL,
#                      xlim, ylim, add = TRUE, legend = TRUE, ...) {
#
#   if (is(dat, "RGChannelSet") || is(dat, "MethylSet")) {
#     b <- getBeta(dat)
#   } else if (is(dat, "matrix")) {
#     b <- dat
#   } else {
#     stop("argument 'dat' must be an 'RGChannelSet', a 'MethylSet' or ",
#          "matrix.")
#   }
#   # NOTE: Have to ensure data are in memory as an ordinary vector before
#   #       computing density
#   d <- apply(b, 2, function(x) density(as.vector(x), na.rm = TRUE))
#   if (missing(ylim)) ylim <- range(sapply(d, function(i) range(i$y)))
#   if (missing(xlim)) xlim <- range(sapply(d, function(i) range(i$x)))
#   if (is.null(sampGroups)) {
#     sampGroups <- rep(1, ncol(b))
#   } else if (length(sampGroups) == 1) {
#     sampGroups <- rep(sampGroups, ncol(b))
#   }
#   sampGroups <- as.factor(sampGroups)
#   if(is.null(cols))cols<-get_cols(sampGroups)
#
#   # Save current graphical parameters
#   opar <- par(no.readonly = TRUE)
#
#   # Make sure d and sampGroups are in the same order:
#
#   # Change the margins of the plot (the first is the bottom margin)
#   par(mar = c(8, 4.1, 4.1, 2.1))
#
#   if (add) {
#     plot(x = 0,
#          type = "n",
#          ylim = ylim,
#          xlim = xlim,
#          ylab = "Density",
#          xlab = xlab,
#          main = main, ...)
#     abline(h = 0, col = "grey80")
#   }
#   for (i in seq_along(d)) {
#     lines(d[[i]], col = cols[i])
#   }
#
#   if (legend & length(levels(sampGroups)) > 1) {
#     legend(x = "bottom",
#            inset = c(0, -0.25), # You will need to fine-tune the second
#            # value depending on the windows size
#            legend = levels(sampGroups),
#            text.col = cols[!duplicated(cols)],
#            # text.width=0.001,
#            cex=0.7, pch=1, pt.cex = 1,
#            lwd = 2,
#            xpd = TRUE, # You need to specify this graphical parameter to add
#            # the legend outside the plot area
#            horiz = TRUE) # Horizontal legend. You can also set the number
#     # of columns with the argument ncol
#     # if horiz = FALSE
#
#     # Back to the default graphical parameters
#
#   }
#   on.exit(par(opar))
# }
