% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chipSeeker.R
\name{summitsChipSeekerAnnote}
\alias{summitsChipSeekerAnnote}
\title{apply ChipSeeker annote for summits result from MACS2}
\usage{
summitsChipSeekerAnnote(
  summitFiles = NULL,
  summitNames = NULL,
  filterChr = TRUE,
  genome = "hg38",
  outdir = "chipSeeker",
  plot = FALSE
)
}
\arguments{
\item{summitFiles}{summitFiles rds file from CallSummitsMACS}

\item{summitNames}{relativte name for summitFiles}

\item{filterChr}{whether file GRanges object}

\item{genome}{genome name,hg38 or mm10}

\item{outdir}{output path to save result}

\item{plot}{whether plot}
}
\value{
peakAnnoList object
}
\description{
apply ChipSeeker annote for summits result from MACS2
}
