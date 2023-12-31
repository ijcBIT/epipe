% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corpca.R
\name{corpca}
\alias{corpca}
\title{Perform Correlation-based PCA and Generate Plot}
\usage{
corpca(
  beta_top100,
  metadata,
  vars = NULL,
  idcol = "barcode",
  path = "./",
  filename = "",
  title = "PC1-6 clinical correlations"
)
}
\arguments{
\item{beta_top100}{A matrix of beta values.}

\item{metadata}{A data frame containing sample metadata.}

\item{vars}{Variables in the metadata to include in the eigencorplot.}

\item{idcol}{Column name in metadata containing unique identifiers for samples (default: "barcode").}

\item{path}{Directory path for saving the plot (default: "./").}

\item{filename}{Name of the file to save the plot (default: "").}

\item{title}{Title for the eigencorplot (default: 'PC1-6 clinical correlations').}
}
\value{
A PCAtools eigencorplot object.
}
\description{
This function performs a correlation-based PCA on beta values and generates an eigencorplot.
}
