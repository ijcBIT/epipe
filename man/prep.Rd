% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep.R
\name{prep}
\alias{prep}
\title{Preprocess Methylation Data}
\usage{
prep(mSetSqn, remove_sex = TRUE, sexchr = c("chrX", "chrY"), arraytype = NULL)
}
\arguments{
\item{mSetSqn}{MethylSet object to be preprocessed.}

\item{remove_sex}{Logical, whether to remove sex chromosomes. Default is TRUE.}

\item{sexchr}{Vector of chromosome names representing sex chromosomes (default is c("chrX", "chrY")).}

\item{arraytype}{Type of array used for preprocessing. If provided, array-specific annotation is used.}
}
\value{
Preprocessed MethylSet object.
}
\description{
This function performs preprocessing steps on a MethylSet object.
}
