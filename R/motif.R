.summarizeJASPARMotifs <- function(motifs = NULL){

  motifNames <- lapply(seq_along(motifs), function(x){
    namex <- make.names(motifs[[x]]@name)
    if(substr(namex,nchar(namex),nchar(namex))=="."){
      namex <- substr(namex,1,nchar(namex)-1)
    }
    namex <- paste0(namex, "_", x)
    namex
  }) %>% unlist(.)

  motifDF <- lapply(seq_along(motifs), function(x){
    data.frame(
      row.names = motifNames[x],
      name = motifs[[x]]@name[[1]],
      ID = motifs[[x]]@ID,
      strand = motifs[[x]]@strand,
      symbol = ifelse(!is.null(motifs[[x]]@tags$symbol[1]), motifs[[x]]@tags$symbol[1], NA) ,
      family = ifelse(!is.null(motifs[[x]]@tags$family[1]), motifs[[x]]@tags$family[1], NA),
      alias = ifelse(!is.null(motifs[[x]]@tags$alias[1]), motifs[[x]]@tags$alias[1], NA),
      stringsAsFactors = FALSE
    )
  }) %>% Reduce("rbind", .) %>% DataFrame

  names(motifs) <- motifNames

  out <- list(motifs = motifs, motifSummary = motifDF)

  return(out)

}

.summarizeChromVARMotifs <- function(motifs = NULL){

  motifNames <- lapply(seq_along(motifs), function(x){
    namex <- make.names(motifs[[x]]@name)
    if(grepl("LINE", namex)){
      splitNamex <- stringr::str_split(motifs[[x]]@ID, pattern="\\_", simplify = TRUE)
      namex <- splitNamex[1, grep("LINE",splitNamex[1,]) + 1]
    }
    if(substr(namex,nchar(namex),nchar(namex))=="."){
      namex <- substr(namex,1,nchar(namex)-1)
    }
    namex <- paste0(namex, "_", x)
    namex
  }) %>% unlist(.)

  motifNames2 <- lapply(seq_along(motifs), function(x){
    namex <- make.names(motifs[[x]]@name)
    if(grepl("LINE", namex)){
      splitNamex <- stringr::str_split(motifs[[x]]@ID, pattern="\\_", simplify = TRUE)
      namex <- splitNamex[1, grep("LINE",splitNamex[1,]) + 1]
    }
    if(substr(namex,nchar(namex),nchar(namex))=="."){
      namex <- substr(namex,1,nchar(namex)-1)
    }
    namex
  }) %>% unlist(.)

  motifDF <- lapply(seq_along(motifs), function(x){
    df <- data.frame(
      row.names = motifNames[x],
      name = motifNames2[[x]],
      ID = motifs[[x]]@ID,
      strand = motifs[[x]]@strand,
      stringsAsFactors = FALSE
    )
  }) %>% Reduce("rbind", .) %>% DataFrame

  names(motifs) <- motifNames

  out <- list(motifs = motifs, motifSummary = motifDF)

  return(out)

}

.safeSaveRDS<-function (object = NULL, file = "", ascii = FALSE, version = NULL,
                        compress = TRUE, refhook = NULL)
{
  testDF <- data.frame(a = 1, b = 2)
  canSave <- suppressWarnings(tryCatch({
    saveRDS(object = testDF, file = file, ascii = ascii,
            version = version, compress = compress, refhook = refhook)
    TRUE
  }, error = function(x) {
    FALSE
  }))
  if (!canSave) {
    dirExists <- dir.exists(dirname(file))
    if (dirExists) {
      stop("Cannot saveRDS. File Path : ", file)
    }
    else {
      stop("Cannot saveRDS because directory does not exist (",
           dirname(file), "). File Path : ", file)
    }
  }
  else {
    saveRDS(object = object, file = file, ascii = ascii,
            version = version, compress = compress, refhook = refhook)
  }
}


#' Add motif annotations to an ArchRProject
#'
#' This function adds information about which peaks contain motifs to a given peakSet Granges Object. For each peak, a binary value
#' is stored indicating whether each motif is observed within the peak region.
#'
#' @param peakSer An `Granges` object.
#' @param motifSet The motif set to be used for annotation. Options include: (i) "JASPAR2016", "JASPAR2018", "JASPAR2020"
#' which gives the 2016, 2018 or 2020 version of JASPAR motifs or (ii) one of "cisbp", "encode", or "homer" which gives the
#' corresponding motif sets from the `chromVAR` package.
#' @param name The name of the `peakAnnotation` object to be stored in the provided `ArchRProject`
#' @param species The name of the species relevant to the supplied `ArchRProject`. This is used for identifying which motif to be
#' used from CisBP/JASPAR. By default, this function will attempt to guess the species based on the value from `getGenome()`.
#' @param collection If one of the JASPAR motif sets is used via `motifSet`, this parameter allows you to indicate the JASPAR
#' collection to be used. See `getMatrixSet()` from `TFBSTools` for all options to supply for collection.
#' @param motifPWMs A custom set of motif PWMs as a PWMList for adding motif annotations.
#' @param cutOff The p-value cutoff to be used for motif search. The p-value is determined vs a background set of sequences
#' (see `MOODS` for more details on this determination).
#' @param width The width in basepairs to consider for motif matches. See the `motimatchr` package for more information.
#' @param version An integer specifying version 1 or version 2 of chromVARmotifs see github for more info GreenleafLab/chromVARmotifs.#' it already exists in the given `ArchRProject`.
#' @param ... Additional parameters to be passed to `TFBSTools::getMatrixSet` for getting a PWM object.
#' @export
makeMotifAnnotations <- function(
  peakSet=NULL,
  motifSet = "cisbp",
  name = "Motif",
  species = NULL,
  genome="BSgenome.Mmusculus.UCSC.mm10",
  collection = "CORE",
  motifPWMs = NULL,
  cutOff = 5e-05,
  width = 7,
  version = 2,
  save=TRUE,
  ...
){

  if(!is.null(motifPWMs)){
    if(!is(motifPWMs, "PWMatrixList")){
      stop("User Supplied motifPWMS must be a PWMatrixList!")
    }
    motifSet <- "Custom"
  }

  if(is.null(motifSet)){
    stop("Must provide motifSet or motifPWMs!")
  }

  require("motifmatchr")

  if(grepl("JASPAR|CISBP", motifSet, ignore.case = TRUE) & is.null(species)){
    if(grepl("hg19",genome, ignore.case = TRUE)){
      species <- "Homo sapiens"
    }
    if(grepl("hg38",genome, ignore.case = TRUE)){
      species <- "Homo sapiens"
    }
    if(grepl("mm9",genome, ignore.case = TRUE)){
      species <- "Mus musculus"
    }
    if(grepl("mm10",genome, ignore.case = TRUE)){
      species <- "Mus musculus"
    }
  }

  #############################################################
  # Get PWM List adapted from chromVAR!
  #############################################################


  if(tolower(motifSet)=="jaspar2020"){

    require("JASPAR2020")
    args <- list(species = species, collection = collection, ...)
    motifs <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, args)
    obj <- .summarizeJASPARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary

  }else if(tolower(motifSet)=="jaspar2018"){

    require("JASPAR2018")
    args <- list(species = species, collection = collection, ...)
    motifs <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, args)
    obj <- .summarizeJASPARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary

  }else if(tolower(motifSet)=="jaspar2016"){

    require("JASPAR2016")
    args <- list(species = species, collection = collection, ...)
    motifs <- TFBSTools::getMatrixSet(JASPAR2016::JASPAR2016, args)
    obj <- .summarizeJASPARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary

  }else if(tolower(motifSet)=="cisbp"){

    require("chromVARmotifs")
    if(tolower(species) == "mus musculus"){
      if(version == 1){
        message("Using version 1 motifs!")
        data("mouse_pwms_v1")
        motifs <- mouse_pwms_v1
      }else if(version == 2){
        message("Using version 2 motifs!")
        data("mouse_pwms_v2")
        motifs <- mouse_pwms_v2
      }else{
        stop("Only versions 1 and 2 exist!")
      }
      obj <- .summarizeChromVARMotifs(motifs)
      motifs <- obj$motifs
      motifSummary <- obj$motifSummary
    }else if(tolower(species) == "homo sapiens"){
      if(version == 1){
        message("Using version 1 motifs!")
        data("human_pwms_v1")
        motifs <- human_pwms_v1
      }else if(version == 2){
        message("Using version 2 motifs!")
        data("human_pwms_v2")
        motifs <- human_pwms_v2
      }else{
        stop("Only versions 1 and 2 exist!")
      }
      obj <- .summarizeChromVARMotifs(motifs)
      motifs <- obj$motifs
      motifSummary <- obj$motifSummary
    }else{
      stop("Species not recognized homo sapiens, mus musculus supported by CisBP!")
    }

  }else if(tolower(motifSet)=="encode"){

    require("chromVARmotifs")
    data("encode_pwms")
    motifs <- encode_pwms
    obj <- .summarizeChromVARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary

  }else if(tolower(motifSet)=="homer"){

    requirePackage("chromVARmotifs")
    data("homer_pwms")
    motifs <- homer_pwms
    obj <- .summarizeChromVARMotifs(motifs)
    motifs <- obj$motifs
    motifSummary <- obj$motifSummary

  }else if(tolower(motifSet)=="custom"){

    obj <- NULL
    motifs <- motifPWMs
    motifSummary <- NULL

  }else{

    stop("Error MotifSet Not Recognized!")

  }

  #############################################################
  # Get BSgenome Information!
  #############################################################
  require(genome, character.only = TRUE)
  BSgenome <- eval(parse(text = genome))
  BSgenome <- validBSgenome(BSgenome)

  #############################################################
  # Calculate Motif Positions
  #############################################################
  if(is.null(peakSet)){
    stop("peakSet is NULL. You need a peakset to run addMotifAnnotations! See addReproduciblePeakSet!")
  }
  motifPositions <- motifmatchr::matchMotifs(
    pwms = motifs,
    subject = peakSet,
    genome = BSgenome,
    out = "positions",
    p.cutoff = cutOff,
    w = width
  )

  #############################################################
  # Motif Overlap Matrix
  #############################################################
  allPositions <- unlist(motifPositions)
  overlapMotifs <- findOverlaps(peakSet, allPositions, ignore.strand=TRUE)
  motifMat <- Matrix::sparseMatrix(
    i = queryHits(overlapMotifs),
    j = match(names(allPositions),names(motifPositions))[subjectHits(overlapMotifs)],
    x = rep(TRUE, length(overlapMotifs)),
    dims = c(length(peakSet), length(motifPositions))
  )
  colnames(motifMat) <- names(motifPositions)
  motifMat <- SummarizedExperiment::SummarizedExperiment(assays=SimpleList(matches = motifMat), rowRanges = peakSet)

  out <- SimpleList(
    motifSummary = motifSummary,
    motifMatches = motifMat,   #  can use for computeDeviations's annotation
    motifPositions = motifPositions,
    motifList = motifs,
    date = Sys.Date()
  )
  if(save){
    outdir<-file.path("addMotif","Annotations")
    makedir(outdir)
    savePositions <- file.path(outdir, paste0(name,"-Positions-In-Peaks.rds"))
    saveMatches <- file.path(outdir,  paste0(name,"-Matches-In-Peaks.rds"))

    .safeSaveRDS(out,file.path(outdir, paste0(name,"-In-Peaks-Summary.rds")), compress = FALSE)
    .safeSaveRDS(out$motifPositions, savePositions, compress = FALSE)
    .safeSaveRDS(out$motifMatches, saveMatches, compress = FALSE)
  }
  return(out)
}

#' get valid BSgenome
#' @param genome genome name,eg.BSgenome.Hsapiens.UCSC.hg19
#' export
validBSgenome<-function (genome = NULL, masked = FALSE)
{
  stopifnot(!is.null(genome))
  if (inherits(genome, "BSgenome")) {
    return(genome)
  }
  else if (is.character(genome)) {
    genome <- tryCatch({
      require(genome, character.only = TRUE)
      bsg <- eval(parse(text = genome))
      if (inherits(bsg, "BSgenome")) {
        return(bsg)
      }
      else {
        stop("genome is not a BSgenome valid class!")
      }
    }, error = function(x) {
      BSgenome::getBSgenome(genome, masked = masked)
    })
    return(genome)
  }
  else {
    stop("Cannot validate BSgenome options are a valid BSgenome or character for getBSgenome")
  }
}

#' get RangedSummarizedExperiment Object Deviation Variability
#' @param seCounts RangedSummarizedExperiment Object
#' @param matchesAnno motif matches object from function `makeMotifAnnotations`
#' @param genome genome name,default:BSgenome.Hsapiens.UCSC.hg19
#' @param n number of top rank motifs to show
#' @param plot whether plot default true
#' @examples
#' data(example_counts, package = "chromVAR")
#' peakset=rowRanges(example_counts)
#' out=makeMotifAnnotations(peakSet=peakset,genome="BSgenome.Hsapiens.UCSC.hg19")
#' p=plotDeVariability(seCounts=example_counts,matchesAnno=out[["motifMatches"]])
#' @export
DeVariability<-function(seCounts=NULL,
                            matchesAnno=NULL,
                            genome="BSgenome.Hsapiens.UCSC.hg19",
                            n=25,
                            plot=TRUE){
  require("chromVAR")
  require(genome, character.only = TRUE)
  BSgenome <- eval(parse(text = genome))
  BSgenome <- validBSgenome(BSgenome)

  message("addGC Bias..")
  seCounts<- addGCBias(seCounts, genome = BSgenome)
  message("computeDeviations..")
  if(is.null(matchesAnno)){
    peakset<-rowRanges(seCounts)
    out<-makeMotifAnnotations(peakSet = peakSet,genome = genome)
    matchesAnno<-out[["motifMatches"]]
  }
  dev <- computeDeviations(seCounts,annotations=matchesAnno)
  message("computeVariability")
  variability <- computeVariability(dev)


  rowV <- variability[order(variability$variability, decreasing = TRUE),]
  rowV$rank <- seq_len(nrow(rowV))

  res<-list("variability"=rowV,"dev"=dev)
  if(plot){
    message("make plot..")
    options(ggrepel.max.overlaps = Inf)
    p<-ggplot(rowV, aes(x = rank, y = variability, color = variability)) +
      geom_point(size = 1) +
      scale_color_gradientn(colors = yipCat::paletteContinuous(set = "comet")) +
      ggrepel::geom_label_repel(
        data = rowV[rev(seq_len(n)), ], aes(x = rank, y = variability, label =name),
        size = 1.5,
        color = "black",
        nudge_x = 2
      ) + yipCat::theme_ArchR() + ylab("Variability") + xlab("Rank Sorted Annotations")
    res$plot<-p
  }
  res
}


#' get row Sums for RangedSummarizedExperiment object
#' @param rseObj RangedSummarizedExperiment object
#' @param seqnames seqnames use
.getRowSums<-function (rseObj=NULL, Seqnames = NULL,
                       verbose = TRUE, filter0 = FALSE, threads = 1)
{
  if (is.null(Seqnames)) {
    Seqnames <- seqlevels(rseObj)
  }
  peaks<-rowRanges(rseObj)
  metadata<-as.data.frame(colData(rseObj))
  if(!"Sample"%in%colnames(metadata)){
    stop("Please add Sample message in colData of Object!!!")
  }
  Samples<-unique(metadata$Sample)
  Counts<-assay(rseObj)
  SeqnamesList<-as.character(seqnames(peaks))
  summaryDF <- .safelapply(seq_along(Seqnames), function(x) {
    chr<- Seqnames[x]
    idx<- which(SeqnamesList==chr)
    for (y in seq_along(Samples)) {
      sam<-Samples[y]
      idy<-which(metadata$Sample==sam)
      count<-Counts[idx,idy,drop=FALSE]

      if (y == 1) {
        sumy<-Matrix::rowSums(count)
      }
      else {

        sumy1<-Matrix::rowSums(count)
        if (length(sumy1) != length(sumy)) {
          stop("rowSums lengths do not match in RSE Object for a seqname!")
        }
        else {
          sumy <- sumy + sumy1
        }
      }
    }
    DataFrame(seqnames = Rle(Seqnames[x], lengths = length(sumy)),
              idx = seq_along(sumy), rowSums = as.vector(sumy),start=start(peaks)[idx],end=end(peaks)[idx])
  }, threads = threads) %>% Reduce("rbind", .)

  if (filter0) {
    summaryDF <- summaryDF[which(summaryDF$rowSums > 0),
                           , drop = FALSE]
  }
  return(summaryDF)
}

#' chromVar method to caculate Backgroup Peaks
#' @param rseObj RangedSummarizedExperiment Object
#' @param Seqnames seqnames use
#' @examples
#' data(example_counts, package = "chromVAR")
#' example_counts=addGCBias(example_counts,genome="BSgenome.Hsapiens.UCSC.hg19")
#' .computeBgdPeaks(example_counts)
#' @return RangedSummarizedExperiment Object
#' @export
computeBgdPeaks<-function (rseObj = NULL, Seqnames=NULL,nIterations = 50, w = 0.1, binSize = 50,
                            seed = 1, outFile = "Background-Peaks.rds")
{
  set.seed(1)
  require("chromVAR")

  rS <- suppressMessages(.getRowSums(rseObj = rseObj,
                                     Seqnames = Seqnames, filter0 = FALSE))
  peakset<-rowRanges(rseObj)
  rS$start <- start(peakset)
  rS$end <- end(peakset)
  GCcol<-c("bias","GC")[c("bias","GC")%in%colnames(mcols(peakset))]
  if(is.null(GCcol)){
    stop("Please compute GC first")
  }
  rS$GC <- mcols(peakset)[[GCcol]]

  uniqueDist <- unique(rS$end - rS$start)

  se <- SummarizedExperiment::SummarizedExperiment(assays = SimpleList(counts = as.matrix(data.frame(rS$rowSums,
                                                                                                     1))), rowData = DataFrame(bias = rS$GC, start = rS$start,
                                                                                                                               end = rS$end))

  bgdPeaks <- tryCatch({
      chromVAR::getBackgroundPeaks(object = se, bias = rowData(se)$bias,
                                   niterations = nIterations, w = w, bs = binSize)
    }, error = function(e) {
      message("Error with chromVAR::getBackgroundPeaks!")
      .ArchRBdgPeaks(
        object = se,
        bias = rowData(se)$bias,
        nIterations = nIterations
      )

    })

  bgdPeaks <- SummarizedExperiment(assays = SimpleList(bgdPeaks = bgdPeaks),
                                   rowRanges = GRanges(rS$seqnames, IRanges(rS$start, rS$end),
                                                       value = rS$rowSums, GC = rS$GC))
  biasDF <- data.frame(rowSums = Matrix::rowSums(assay(se)),
                       bias = rowData(se)$bias, length = rowData(se)$end - rowData(se)$start)
  rowData(bgdPeaks)$bgdSumMean <- round(rowMeans(matrix(biasDF[assay(bgdPeaks),
                                                               1], nrow = nrow(bgdPeaks))), 3)
  rowData(bgdPeaks)$bgdSumSd <- round(matrixStats::rowSds(matrix(biasDF[assay(bgdPeaks),
                                                                        1], nrow = nrow(bgdPeaks))), 3)
  rowData(bgdPeaks)$bgdGCMean <- round(rowMeans(matrix(biasDF[assay(bgdPeaks),
                                                              2], nrow = nrow(bgdPeaks))), 3)
  rowData(bgdPeaks)$bgdGCSd <- round(matrixStats::rowSds(matrix(biasDF[assay(bgdPeaks),
                                                                       2], nrow = nrow(bgdPeaks))), 3)
  rowData(bgdPeaks)$bgdLengthMean <- round(rowMeans(matrix(biasDF[assay(bgdPeaks),
                                                                  3], nrow = nrow(bgdPeaks))), 3)
  rowData(bgdPeaks)$bgdLengthSd <- round(matrixStats::rowSds(matrix(biasDF[assay(bgdPeaks),
                                                                           3], nrow = nrow(bgdPeaks))), 3)
  if (!is.null(outFile)) {
    saveRDS(bgdPeaks, outFile, compress = FALSE)
  }
  return(bgdPeaks)
}

.ArchRBdgPeaks <- function(object = NULL, bias = NULL, nIterations = 50){

  .cleanSelf <- function(x){
    xn <- matrix(0, nrow = nrow(x), ncol = ncol(x))
    for(i in seq_len(nrow(x))){
      xi <- x[i, ]
      idx <- which(xi != i)
      xn[i, seq_along(idx)] <- xi[idx]
    }
    idx <- which(colSums(xn == 0) > 0)
    if(length(idx) > 0){
      xn <- xn[,-idx]
    }
    xn
  }

  #Bias Dataframe
  biasDF <- data.frame(
    rowSums = Matrix::rowSums(assay(object)),
    bias = bias,
    length = rowData(object)$end - rowData(object)$start
  )

  #Quantile Normalize
  biasDFN <- apply(biasDF, 2, .getQuantiles)

  #Get KNN
  knnObj <- nabor::knn(
    data =  biasDFN,
    k = nIterations + 1
  )[[1]]

  #Filter Self
  knnObj <- .cleanSelf(knnObj)

  #Shuffle
  idx <- seq_len(ncol(knnObj))
  knnObj2 <- matrix(0, nrow = nrow(knnObj), ncol = ncol(knnObj))
  for(x in seq_len(nrow(knnObj2))){
    knnObj2[x,] <- knnObj[x, sample(idx, length(idx))]
  }

  knnObj2

}


.customDeviationsSingle <- function(
  annotationsVector = NULL,
  countsMatrix = NULL,
  countsPerSample = NULL,
  backgroudPeaks = NULL,
  out = c("deviations", "z"),
  expectation = NULL,
  threshold = 1,
  prefix = ""
){

  .binarizeMat <- function(mat = NULL){
    mat@x[mat@x > 0] <- 1
    mat
  }

  if (length(annotationsVector@x) == 0) {
    out <- list(
      z = rep(NA, ncol(countsMatrix)),
      dev = rep(NA, ncol(countsMatrix)),
      expFG = NA,
      expBG = NA,
      matches = 0,
      overlap = NA
    )
    return(out)
  }

  outList <- tryCatch({

    ################################
    # Fore Ground Deviations
    ################################
    observed <- as.vector(Matrix::t(annotationsVector) %*% countsMatrix)
    expected <- as.vector(Matrix::t(annotationsVector) %*% expectation %*% countsPerSample)
    observed_deviation <- (observed - expected)/expected

    #Filter those with no matches at all
    fail_filter <- which(expected == 0)

    ################################
    # Back Ground Deviations
    ################################
    if("z" %in% tolower(out)){

      #Compute Background Null Per Iteration
      niterations <- ncol(backgroudPeaks)
      sampleMat <- Matrix::sparseMatrix(
        j = as.vector(backgroudPeaks[annotationsVector@i + 1, seq_len(niterations)]),
        i = rep(seq_len(niterations), each = length(annotationsVector@x)),
        x = rep(annotationsVector@x, niterations),
        dims = c(niterations, nrow(countsMatrix))
      )
      sampled <- as.matrix(sampleMat %*% countsMatrix)
      sampledExpected <- sampleMat %*% expectation %*% countsPerSample
      sampledDeviation <- (sampled - sampledExpected)/sampledExpected
      bgOverlap <- Matrix::mean(.binarizeMat(sampleMat) %*% .binarizeMat(annotationsVector)) / length(annotationsVector@x)

      #Summary
      meanSampledDeviation <- Matrix::colMeans(sampledDeviation)
      sdSampledDeviation <- apply(as.matrix(sampledDeviation), 2, sd)

      #Norm Deviation
      normdev <- (observed_deviation - meanSampledDeviation)
      z <- normdev/sdSampledDeviation
      if (length(fail_filter) > 0) {
        z[fail_filter] <- NA
        normdev[fail_filter] <- NA
      }

    }else{

      #Compute Background Null Per Iteration
      niterations <- ncol(backgroudPeaks)
      sampleMat2 <- Matrix::sparseMatrix(
        j = as.vector(backgroudPeaks[annotationsVector@i + 1, seq_len(niterations)]),
        i = rep(1, niterations * length(annotationsVector@x)),
        x = rep(annotationsVector@x, niterations),
        dims = c(1, nrow(countsMatrix))
      )
      sampled2 <- (sampleMat2 %*% countsMatrix)[1,]
      sampledExpected2 <- (sampleMat2 %*% expectation %*% countsPerSample)[1,]
      ######################
      # Equivalent to above
      # colMeans(sampled) - colMeans(sampledExpected))/colMeans(sampledExpected)
      ######################
      sampledDeviation2 <- (sampled2 - sampledExpected2)/sampledExpected2
      bgOverlap <- NA

      #Norm Deviation
      normdev <- (observed_deviation - sampledDeviation2)
      z <- NULL
      if (length(fail_filter) > 0) {
        normdev[fail_filter] <- NA
      }

    }

    outList <- list(
      z = z,
      dev = normdev,
      matches = length(annotationsVector@x) / nrow(countsMatrix),
      overlap = bgOverlap
    )

    outList

  }, error = function(e){
    message(".customDeviationsSingle error!")
    NULL
  })

  return(outList)

}



############################################################################
# Adapted from chromVAR, Approved by Alicia Schep for Modification.
############################################################################
#' add Deviation Matrix
#' @param rseObj RangedSummarizedExperiment Object
#' @param matches object from makeMotifAnnotations motifMatches
#' @param threads  number of threads
#' @return SummarizedExperiment object
#' @export
getDeviationsMatrix <- function(
  rseObj = NULL,
  matches  = NULL,
  prefix = "",
  out = c("deviations", "z"),
  threads = 1,
  verbose = TRUE
){
  message("get expectation..")
  peakset<-rowRanges(rseObj)
  rS<-.getRowSums(rseObj = rseObj,filter0 = FALSE,threads = threads)
  rownames(rS) <- paste(rS$seqnames,rS$start,rS$end,sep="_")
  #rS <- rS[paste0(seqnames(peakset), "_", mcols(peakset)$idx),]
  rS<-rS[paste(seqnames(peakset),start(peakset),end(peakset),sep="_"),]
  #rS$start <- start(peakset)
  #rS$end <- end(peakset)
  GCcol<-c("bias","GC")[c("bias","GC")%in%colnames(mcols(peakset))]
  if(length(GCcol)==0){
    stop("Please compute GC first")
  }
  rS$GC <- mcols(peakset)[[GCcol]]
  #rownames(rS) <- paste0(rS$seqnames, "_", rS$start, "_", rS$end)
  expectation<-rS$rowSums/sum(rS$rowSums)

  message("get backgroudPeaks.. ")
  bgdPeaks<-computeBgdPeaks(rseObj)
  backgroudPeaks<-SummarizedExperiment::assay(bgdPeaks)


  message("get matches annotation matrix..")

  rownames(matches)<-paste0(seqnames(matches), "_", start(matches), "_", end(matches))
  macthes<-matches[rownames(rS),]
  annotationsMatrix <- SummarizedExperiment::assay(matches)
  annotationsMatrix <- as(annotationsMatrix, "dgCMatrix")
  gc()


  countsMatrix<-assay(rseObj)
  colData <- DataFrame(seq_len(ncol(countsMatrix)), row.names = colnames(countsMatrix))[,FALSE]
  norm_expectation <- expectation / sum(expectation) #Double check this sums to 1!
  countsPerSample <- Matrix::colSums(countsMatrix)


  d <- max(floor(ncol(annotationsMatrix)/20), 1)
  m <- 0
  results <- .safelapply(seq_len(ncol(annotationsMatrix)), function(x){
    if(x %% d == 0){
      cat(sprintf("%s : Deviations for Annotation %s of %s\n", prefix, x, ncol(annotationsMatrix)))
      m <- 1 #Print to console
    }
    if(x %% max(floor(d/5), 2) == 0){
      if(m != 1){
        cat(sprintf("%s : Deviations for Annotation %s of %s\n", prefix, x, ncol(annotationsMatrix)))
      }else{
        m <- 0 #Reset
      }
    }
    if(x %% max(c(d, 10)) == 0){
      gc()
    }
    .customDeviationsSingle(
      annotationsVector = annotationsMatrix[, x, drop=FALSE],
      countsMatrix = countsMatrix,
      backgroudPeaks = backgroudPeaks,
      countsPerSample = countsPerSample,
      expectation = norm_expectation,
      out = out,
      prefix = prefix
    )
  }, threads = threads)
  cn <- colnames(countsMatrix)
  rm(countsMatrix)
  gc()

  if("z" %in% tolower(out)){
    z <- t(vapply(results, function(x) x[["z"]], rep(0, length(cn))))
    if(length(cn)==1){
      z <- matrix(z, ncol=length(cn))
    }
  }else{
    z <- matrix(0, nrow = ncol(annotationsMatrix), ncol = length(cn))
  }
  if("deviations" %in% tolower(out)){
    dev <- t(vapply(results, function(x) x[["dev"]], rep(0, length(cn))))
    if(length(cn)==1){
      dev <- matrix(dev, ncol=length(cn))
    }
  }else{
    dev <- matrix(0, nrow = ncol(annotationsMatrix), ncol = length(cn))
  }
  colnames(z) <- cn
  colnames(dev) <- cn

  #Check First
  nullOverlap <- is.null(results[[1]]$overlap)
  rowData <- lapply(seq_along(results), function(x){
    resx <- results[[x]]
    if(nullOverlap){
      data.frame(fractionMatches = resx$matches)
    }else{
      data.frame(fractionMatches = resx$matches, fractionBackgroundOverlap = resx$overlap)
    }
  }) %>% Reduce("rbind",.)
  rownames(rowData) <- colnames(annotationsMatrix)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      deviations = dev,
      z = z
    ),
    colData = colData,
    rowData = rowData
  )
  SummarizedExperiment::assays(se) <- SummarizedExperiment::assays(se)[tolower(out)]

  return(se)

}
