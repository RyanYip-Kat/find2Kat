% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreateGeneAnnotation.R
\name{createGenomeAnnotation}
\alias{createGenomeAnnotation}
\title{create Genome Annotation}
\usage{
createGenomeAnnotation(
  genome = NULL,
  chromSizes = NULL,
  blacklist = NULL,
  filter = TRUE,
  filterChr = c("chrM")
)
}
\arguments{
\item{genome}{genome object}

\item{chromSizes}{chromSizes value}

\item{blacklist}{blacklist Grange  object}

\item{filter}{whether filter chroms}

\item{filterChr}{the chroms to be filtered}
}
\description{
create Genome Annotation
}
