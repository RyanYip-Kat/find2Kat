% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ArchRUtils.R
\name{DOCRPointPlot}
\alias{DOCRPointPlot}
\title{DOCR Score Rank Point Plot}
\usage{
DOCRPointPlot(
  p2g = NULL,
  genes = NULL,
  nShowGenes = 10,
  save = FALSE,
  outDir = NULL
)
}
\arguments{
\item{p2g}{object from ArchR::getPeak2GeneLinks}

\item{genes}{names,from ArchR::getFeatures(proj,useMatrix="GeneScoreMatrix")}

\item{nShowGenes}{number of gene name to label in point plot}
}
\value{
ggplot object
}
\description{
DOCR Score Rank Point Plot
}
