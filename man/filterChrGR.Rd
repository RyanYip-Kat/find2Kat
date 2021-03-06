% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreateGeneAnnotation.R, R/annoteBED.R
\name{filterChrGR}
\alias{filterChrGR}
\title{filter Chrom GRanges}
\usage{
filterChrGR(
  gr = NULL,
  remove = NULL,
  underscore = TRUE,
  standard = TRUE,
  pruningMode = "coarse"
)

filterChrGR(
  gr = NULL,
  remove = NULL,
  underscore = TRUE,
  standard = TRUE,
  pruningMode = "coarse"
)
}
\arguments{
\item{gr}{A `GRanges` object or another object containing `seqlevels`.}

\item{remove}{A character vector indicating the seqlevels that should be removed if manual removal is desired for certain seqlevels.
If no manual removal is desired, `remove` should be set to `NULL`.}

\item{underscore}{A boolean value indicating whether to remove all seqlevels whose names contain an underscore (for example "chr11_KI270721v1_random").}

\item{standard}{A boolean value indicating whether only standard chromosomes should be kept. Standard chromosomes are defined by
`GenomeInfoDb::keepStandardChromosomes()`.}

\item{pruningMode}{The name of the pruning method to use (from`GenomeInfoDb::seqinfo()`) when seqlevels must be removed from a `GRanges` object.
When some of the seqlevels to drop from the given `GRanges` object are in use (i.e. have ranges on them), the ranges on these sequences need
to be removed before the seqlevels can be dropped. Four pruning modes are currently defined: "error", "coarse", "fine", and "tidy".}
}
\description{
This function allows for removal of manually designated or more broadly undesirable seqlevels from a Genomic Ranges object or similar object
}
