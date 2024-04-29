
#' Apply Noob normalization to RGChannelSet.
#'
#' @param rgSet RGChannelSet object to be normalized.
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessNoob
#' @export
noob <- function(rgSet){
  require(minfi)
  minfi::preprocessNoob(rgSet)
}

#' Apply SWAN normalization to RGChannelSet.
#'
#' @param rgSet RGChannelSet object to be normalized.
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessSWAN
#' @export
swan <- function(rgSet){
  require(minfi)
  minfi::preprocessSWAN(rgSet)
}

#' Apply FunNorm normalization to RGChannelSet.
#'
#' @param rgSet RGChannelSet object to be normalized.
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessFunnorm
#' @export
funn <- function(rgSet){
  require(minfi)
  minfi::preprocessFunnorm(rgSet)
}

#' Apply Noob followed by Quantile normalization to RGChannelSet.
#'
#' @param rgSet RGChannelSet object to be normalized.
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessNoob preprocessQuantile
#' @export
noob_pq <- function(rgSet){
  require(minfi)
  minfi::preprocessNoob(rgSet) |> minfi::preprocessQuantile()
}

#' Apply Quantile normalization to RGChannelSet.
#'
#' @param rgSet RGChannelSet object to be normalized.
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessQuantile
#' @export
pq <- function(rgSet){
  require(minfi)
  minfi::preprocessQuantile(rgSet)
}

#' Apply ENmix normalization to RGChannelSet.
#'
#' @param rgSet RGChannelSet object to be normalized.
#' @param arraytype Type of array used for preprocessing. If provided, array-specific annotation is used.
#' @return Normalized GenomicRatioSet object.
#'
#' @importFrom ENmix preprocessENmix QCinfo mpreprocess
#' @importFrom minfi makeGenomicRatioSetFromMatrix
#' @export
Em2 <- function(rgSet, arraytype = NULL){
  pd <- rgSet@colData
  qc <- ENmix::QCinfo(rgSet, distplot = FALSE)
  mdat <- ENmix::preprocessENmix(rgSet, QCinfo = qc)
  mSetSqn <- ENmix::mpreprocess(rgSet = rgSet, impute = TRUE)
  if(is.null(arraytype)){
    arraytype <- ifelse(500000 < nrow(mSetSqn), "EPIC", "450K" )
  }

  if(arraytype == "EPIC"){
    mSetSqn <- minfi::makeGenomicRatioSetFromMatrix(
      mat = mSetSqn,
      array = "IlluminaHumanMethylationEPIC",
      annotation = "ilm10b4.hg19"
    )
  } else {
    mSetSqn <- minfi::makeGenomicRatioSetFromMatrix(mSetSqn)
  }

  mSetSqn@colData <- pd
  return(mSetSqn)
}

#' Apply ENmix normalization to RGChannelSet and map to genome.
#'
#' @param rgSet RGChannelSet object to be normalized.
#' @param arraytype Type of array used for preprocessing. If provided, array-specific annotation is used.
#' @return Normalized GenomicRatioSet object.
#'
#' @importFrom ENmix preprocessENmix QCinfo
#' @importFrom minfi mapToGenome
#' @import SummarizedExperiment
#' @export
Em <- function(rgSet, arraytype = NULL){
  require(SummarizedExperiment)
  pd <- rgSet@colData
  qc <- ENmix::QCinfo(rgSet, distplot = FALSE)
  mdat <- ENmix::preprocessENmix(rgSet, QCinfo = qc)
  mSetSqn <- minfi::mapToGenome(mdat)
  mSetSqn@colData <- pd
  return(mSetSqn)
}




makelist<-function(vals){
  vals$norm<-as.character(vals$norm)
  return(vals)
  
}








