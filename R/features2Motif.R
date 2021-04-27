#' @importFrom Rcpp evalCpp
#' @importFrom methods setClass setMethod is slot slot<- new as
#' slotNames
#' @importClassesFrom Matrix dgCMatrix
#' @useDynLib find2Kat
NULL


#' Create motif object
#'
#' Create a \code{\link{Motif-class}} object.
#'
#' @param data A motif x region matrix
#' @param pwm A named list of position weight matrices or position frequency
#' matrices matching the motif names in \code{data}.
#' Can be of class PFMatrixList.
#' @param motif.names A named list of motif names. List element names
#' must match the names given in \code{pwm}. If NULL, use the names from the
#' list of position weight or position frequency matrices. This can be used to
#' set a alternative common name for the motif. If a PFMatrixList is passed to
#' \code{pwm}, it will pull the motif name from the PFMatrixList.
#' @param positions A \code{\link[GenomicRanges]{GRangesList}} object containing
#' exact positions of each motif.
#' @param meta.data A data.frame containing metadata
#' @export
#' @return Returns a \code{\link{Motif}} object
#' @concept motifs
#' @examples
#' motif.matrix <- matrix(
#'   data = sample(c(0,1),
#'     size = 100,
#'     replace = TRUE),
#'   ncol = 5
#' )
#' motif <- CreateMotifObject(data = motif.matrix)
CreateMotifObject <- function(
  data = NULL,
  pwm = NULL,
  motif.names = NULL,
  positions = NULL,
  meta.data = NULL
) {

  data <- SetIfNull(x = data, y = new(Class = "dgCMatrix"))
  meta.data <- SetIfNull(x = meta.data, y = data.frame())
  if (
    !(inherits(x = data, what = "matrix") |
      inherits(x = data, what = "dgCMatrix"))
  ) {
    stop("Data must be matrix or sparse matrix class. Supplied ",
         class(x = data))
  }
  if (inherits(x = data, what = "matrix")) {
    data <- as(Class = "dgCMatrix", object = data)
  }
  if ((nrow(x = data) > 0) & (length(x = pwm) > 0)) {
    if (!all(names(x = pwm) == colnames(x = data))) {
      stop("Motif names in data matrix and PWM list are inconsistent")
    }
  }
  if ((nrow(x = data) > 0) & (nrow(x = meta.data) > 0)) {
    if (!all(rownames(x = meta.data) == rownames(x = data))) {
      stop("Motif names in data matrix and metadata are inconsistent")
    }
  }
  if (inherits(x = pwm, what = "list")) {
    if (is.null(names(x = pwm))) {
      stop("PWM must be a named list")
    }
  }
  if (!is.null(x = motif.names)) {
    if (length(x = motif.names) != length(x = pwm)) {
      stop("Number of motif names supplied does not match the number of motifs")
    }
  }
  if (
    inherits(x = pwm, what = "PFMatrixList") |
    inherits(x = pwm, what = "PWMatrixList")
  ) {
    pwm.converted <- lapply(X = as.list(x = pwm), FUN = PFMatrixToList)
    pwm <- lapply(X = pwm.converted, FUN = "[[", 1)
    motif.names <- lapply(X = pwm.converted, FUN = "[[", 2)
  }
  pwm <- SetIfNull(x = pwm, y = list())
  if (is.null(x = motif.names)) {
    motif.names <- as.list(x = names(x = pwm))
    names(motif.names) <- names(x = pwm)
  }
  motif.obj <- new(
    Class = "Motif",
    data = data,
    pwm = pwm,
    motif.names = motif.names,
    positions = positions,
    meta.data = meta.data
  )
  return(motif.obj)
}


#' Create motif matrix
#'
#' Create a motif x feature matrix from a set of genomic ranges,
#' the genome, and a set of position weight matrices.
#'
#' Requires that motifmatchr is installed
#' \url{https://www.bioconductor.org/packages/motifmatchr/}.
#'
#' @param features A GRanges object containing a set of genomic features
#' @param pwm A \code{\link[TFBSTools]{PFMatrixList}} or
#' \code{\link[TFBSTools]{PWMatrixList}}
#' object containing position weight/frequency matrices to use
#' @param genome Any object compatible with the \code{genome} argument
#' in \code{\link[motifmatchr]{matchMotifs}}
#' @param score Record the motif match score, rather than presence/absence
#' (default FALSE)
#' @param use.counts Record motif counts per region. If FALSE (default),
#' record presence/absence of motif. Only applicable if \code{score=FALSE}.
#' @param sep A length-2 character vector containing the separators to be used
#' when constructing matrix rownames from the GRanges
#' @param ... Additional arguments passed to
#' \code{\link[motifmatchr]{matchMotifs}}
#'
#' @return Returns a sparse matrix
#' @export
#' @concept motifs
#' @examples
#' \dontrun{
#' library(JASPAR2018)
#' library(TFBSTools)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#'
#' pwm <- getMatrixSet(
#'   x = JASPAR2018,
#'   opts = list(species = 9606, all_versions = FALSE)
#' )
#' motif.matrix <- CreateMotifMatrix(
#'   features = granges(atac_small),
#'   pwm = pwm,
#'   genome = BSgenome.Hsapiens.UCSC.hg38
#' )
#' }
CreateMotifMatrix <- function(
  features,
  pwm,
  genome,
  score = FALSE,
  use.counts = FALSE,
  sep = c("-", "-"),
  ...
) {
  if (!requireNamespace("motifmatchr", quietly = TRUE)) {
    stop("Please install motifmatchr.
         https://www.bioconductor.org/packages/motifmatchr/")
  }
  motif_ix <- motifmatchr::matchMotifs(
    pwms = pwm,
    subject = features,
    genome = genome,
    out = "scores",
    ...
  )
  if (score) {
    motif.matrix <- motifmatchr::motifScores(object = motif_ix)
  } else {
    if (use.counts) {
      motif.matrix <- motifmatchr::motifCounts(object = motif_ix)
    } else {
      motif.matrix <- motifmatchr::motifMatches(object = motif_ix)
      motif.matrix <- as(Class = "dgCMatrix", object = motif.matrix)
    }
  }
  rownames(motif.matrix) <- GRangesToString(grange = features, sep = sep)
  return(motif.matrix)
}

#' GRanges to String
#'
#' Convert GRanges object to a vector of strings
#'
#' @param grange A GRanges object
#' @param sep Vector of separators to use for genomic string. First element is
#' used to separate chromosome and coordinates, second separator is used to
#' separate start and end coordinates.
#' @importMethodsFrom GenomicRanges start end seqnames
#' @examples
#' GRangesToString(grange = blacklist_hg19)
#' @return Returns a character vector
#' @export
#' @concept utilities
GRangesToString <- function(grange, sep = c("-", "-")) {
  regions <- paste0(
    as.character(x = seqnames(x = grange)),
    sep[[1]],
    start(x = grange),
    sep[[2]],
    end(x = grange)
  )
  return(regions)
}

##############################################
## Create Motif Object and Plot for PeakSet ##
##############################################

#' get Motifs DataBase
#' @param species which species to use(eg,"Mus musculus" or "Homo sapiens"),default "Mus musculus
#' @param collection which collection to use,deafault CORE
#' @return PWMatrixList Object
#' @export
getMotifsDB=function (species = "Mus musculus", collection = "CORE", ...)
{
  opts <- list()
  opts["species"] <- species
  opts["collection"] <- collection
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out))))
    names(out) <- paste(names(out), TFBSTools::name(out),
                        sep = "_")
  return(out)
}

#' Make Signac Motif Object with Granges Object
#' @param object Granges object,like peakset from ArchR::getPeakSet
#' @param genome BSgenome Object ,eg,BSgenome.Mmusculus.UCSC.mm10
#' @param pfm PWMatrixList object,eg,chromVARmotifs::mouse_pwms_v2
#' @return motif.matrix,and you can plot with like MotifPlot(motif.matrix,motifs =c("Fosb","Mga")
#' @export
GRangesAddMotifs <- function(
  object,
  genome,
  pfm,
  verbose = TRUE,
  ...
) {
  if (!requireNamespace("motifmatchr", quietly = TRUE)) {
    stop("Please install motifmatchr.\n",
         "https://www.bioconductor.org/packages/motifmatchr/")
  }
  if (verbose) {
    message("Building motif matrix")
  }
  motif.matrix <- CreateMotifMatrix(
    features = object,
    pwm = pfm,
    genome = genome,
    use.counts = FALSE
  )
  if (verbose) {
    message("Finding motif positions")
  }
  motif.positions <- motifmatchr::matchMotifs(
    pwms = pfm,
    subject = object,
    out = 'positions',
    genome = genome
  )
  if (verbose) {
    message("Creating Motif object")
  }
  motif <- CreateMotifObject(
    data = motif.matrix,
    positions = motif.positions,
    pwm = pfm
  )
  return(motif)
}


# Set a default value if an object is null
#
# @param x An object to set if it's null
# @param y The value to provide if x is null
# @return Returns y if x is null, otherwise returns x.
SetIfNull <- function(x, y) {
  if (is.null(x = x)) {
    return(y)
  } else {
    return(x)
  }
}

#' Get Motif Data
#' @param slot Information to pull from object (data, pwm, meta.data)
#' @rdname GetMotifData
#' @method GetMotifData Motif
#' @export
#' @concept motifs
#' @examples
#' motif.obj <- GetMotifData(
#'   object = atac_small[['peaks']], slot = "motifs"
#' )
#' GetMotifData(object = motif.obj)
GetMotifData <- function(object, slot = "data", ...) {
  return(slot(object = object, name = slot))
}


#' Plot DNA sequence motif
#'
#' Plot position weight matrix or position frequency matrix for different DNA
#' sequence motifs.
#'
#' @param object A Seurat object
#' @param motifs A list of motifs to plot
#' @param assay Name of the assay to use
#' @param use.names Use motif names stored in the motif object
#' @param ... Additional parameters passed to \code{\link[ggseqlogo]{ggseqlogo}}
#'
#' @importFrom ggseqlogo ggseqlogo
#' @export
#' @concept visualization
#' @concept motifs
#' @return Returns a \code{\link[ggplot2]{ggplot}} object
#' @examples
#' \donttest{
#' motif.obj <- Seurat::GetAssayData(atac_small, slot = "motifs")
#' MotifPlot(atac_small, motifs = head(colnames(motif.obj)))
#' }
MotifPlot <- function(
  object,
  motifs,
  assay = NULL,
  use.names = TRUE,
  ...
) {
  require("ggseqlogo")
  require("stringr")
  motif.names <- GetMotifData(object=object,slot = "motif.names")
  all.motifs <- as.character(unlist(motif.names))
  motifs<-names(motif.names)[str_to_lower(all.motifs)%in%str_to_lower(motifs)]

  data.use <- GetMotifData(object = object, assay = assay, slot = "pwm")
  if (length(x = data.use) == 0) {
    stop("Position weight matrix list for the requested assay is empty")
  }
  data.use <- data.use[motifs]
  if (use.names) {
    names(x = data.use) <- GetMotifData(
      object = object, assay = assay, slot = "motif.names"
    )[motifs]
  }
  p <- ggseqlogo(data = data.use, ...)
  return(p)
}

# Convert PFMMatrix to
# @param x A PFMatrix
PFMatrixToList <- function(x) {
  if (!requireNamespace("TFBSTools", quietly = TRUE)) {
    stop("Please install TFBSTools.
         https://www.bioconductor.org/packages/TFBSTools/")
  }
  position.matrix <- TFBSTools::Matrix(x = x)
  name.use <- TFBSTools::name(x = x)
  return(list("matrix" = position.matrix, "name" = name.use))
}

#' The Motif class
#'
#' The Motif class is designed to store DNA sequence motif information,
#' including motif PWMs or PFMs, motif positions, and metadata.
#'
#' @slot data A sparse, binary, feature x motif matrix. Columns
#' correspond to motif IDs, rows correspond to genomic features
#' (peaks or bins). Entries in the matrix should be 1 if the
#' genomic feature contains the motif, and 0 otherwise.
#' @slot pwm A named list of position weight matrices
#' @slot motif.names A list containing the name of each motif
#' @slot positions A \code{\link[GenomicRanges]{GRangesList}} object containing
#' exact positions of each motif.
#' @slot meta.data A dataframe for storage of additional
#' information related to each motif. This could include the
#' names of proteins that bind the motif.
#'
#' @name Motif-class
#' @rdname Motif-class
#' @exportClass Motif
#' @concept motifs
Motif <- setClass(
  Class = "Motif",
  slots = list(
    data = "dgCMatrix",
    pwm = "list",
    motif.names = "list",
    positions = "ANY",
    meta.data = "data.frame"
  )
)

