% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_dmrs.R
\name{find_dmrs}
\alias{find_dmrs}
\title{Find DMRs}
\usage{
find_dmrs(
  object = NULL,
  betas = NULL,
  model,
  fdr = 0.05,
  p.value = "fdr",
  bcutoff = 0.3,
  min.cpg = 5,
  ncores = NULL
)
}
\arguments{
\item{object}{DMPextr object or NULL. If NULL, it will be constructed using the model and beta values.}

\item{betas}{Beta values matrix. Ignored if object is provided.}

\item{model}{Linear model used for differential methylation analysis.}

\item{fdr}{False Discovery Rate threshold.}

\item{p.value}{P-value threshold (not recommended, use fdr).}

\item{bcutoff}{Beta value cutoff.}

\item{min.cpg}{Minimum number of CpGs in a DMR.}

\item{ncores}{Number of CPU cores to use.}
}
\value{
A data.table containing the identified DMRs.
}
\description{
Find Differentially Methylated Regions (DMRs)
}
\details{
This function identifies differentially methylated regions (DMRs) from beta values and a linear model.
}
\examples{
# Example usage:
model <- lm(betas ~ condition, data = beta_data)
dmrs <- find_dmrs(object = NULL, betas = beta_data, model = model)

}
\author{
Izar de Villasante
}
