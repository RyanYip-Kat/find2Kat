#' filter chroms
#' @gr Granges object
#' @param genome genome name
#' @export
keepGR<-function(gr,genome,pruningMode="coarse"){
  seqNames <- GenomeInfoDb::seqlevels(gr)
  chrRemove<- c()
  if(tolower(genome%in%c("hg38","hg19"))){
    keepChr=paste0("chr",c(1:22,"X"))
  }else{
    keepChr=paste0("chr",c(1:19,"X"))
  }
  chrRemove <- c(chrRemove, which(!seqNames %in% keepChr))
  chrKeep <- seqNames[-chrRemove]
  GenomeInfoDb::seqlevels(gr, pruning.mode=pruningMode) <- chrKeep
  gr
}

#' Filters unwanted seqlevels from a Genomic Ranges object or similar object
#'
#' This function allows for removal of manually designated or more broadly undesirable seqlevels from a Genomic Ranges object or similar object
#'
#' @param gr A `GRanges` object or another object containing `seqlevels`.
#' @param remove A character vector indicating the seqlevels that should be removed if manual removal is desired for certain seqlevels.
#' If no manual removal is desired, `remove` should be set to `NULL`.
#' @param underscore A boolean value indicating whether to remove all seqlevels whose names contain an underscore (for example "chr11_KI270721v1_random").
#' @param standard A boolean value indicating whether only standard chromosomes should be kept. Standard chromosomes are defined by
#' `GenomeInfoDb::keepStandardChromosomes()`.
#' @param pruningMode The name of the pruning method to use (from`GenomeInfoDb::seqinfo()`) when seqlevels must be removed from a `GRanges` object.
#' When some of the seqlevels to drop from the given `GRanges` object are in use (i.e. have ranges on them), the ranges on these sequences need
#' to be removed before the seqlevels can be dropped. Four pruning modes are currently defined: "error", "coarse", "fine", and "tidy".
#' @export
filterChrGR <- function(
  gr = NULL,
  remove = NULL,
  underscore = TRUE,
  standard = TRUE,
  pruningMode="coarse"
){

  #first we remove all non standard chromosomes
  if(standard){
    gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = pruningMode)
  }
  #Then check for underscores or specified remove
  seqNames <- GenomeInfoDb::seqlevels(gr)
  chrRemove <- c()
  #first we remove all chr with an underscore
  if(underscore){
    chrRemove <- c(chrRemove, which(grepl("_", seqNames)))
  }
  #next we remove all chr specified in remove
  if(!is.null(remove)){
    chrRemove <- c(chrRemove, which(seqNames %in% remove))
  }
  if(length(chrRemove) > 0){
    chrKeep <- seqNames[-chrRemove]
  }else{
    chrKeep <- seqNames
  }
  #this function restores seqlevels
  GenomeInfoDb::seqlevels(gr, pruning.mode=pruningMode) <- chrKeep

  return(gr)

}
#' call macs2 for bam file or bed file(pass)
#' @export
callSummitsMACS2 <- function(
  File = NULL,
  format="BED",
  pathToMacs2 = "macs2",
  genomeSize = 2.7e9,
  shift = -75,
  extsize = 150,
  cutOff = 0.05,
  method = "q",
  additionalParams = "--nomodel --nolambda",
  outdir="result",
  prefix="macs2",
  save=TRUE,
  cleanUp=FALSE
){

  stopifnot(tolower(method) %in% c("p","q"))
  stopifnot(!is.null(genomeSize))
  utility <- .checkPath(pathToMacs2)

  #Output Files
  #bedName <- gsub("\\.insertions.bed", "", bedFile)
  bedName=prefix

  #Create MACS2 Command
  cmd <- sprintf("callpeak -g %s --name %s --treatment %s --outdir %s --format %s --call-summits --keep-dup all %s",
                 genomeSize, bedName,File,outdir,format, additionalParams)

  if(!is.null(shift) & !is.null(extsize)){
    cmd <- sprintf("%s --shift %s --extsize %s", cmd , shift, extsize)
  }

  if(tolower(method) == "p"){
    cmd <- sprintf("%s -p %s", cmd , cutOff)
  }else{
    cmd <- sprintf("%s -q %s", cmd , cutOff)
  }


  run <- system2(pathToMacs2, cmd, wait=TRUE, stdout=NULL, stderr=NULL)
  summitsFile <- file.path(outdir,paste0(bedName, "_summits.bed"))
  narrowPeaksFile <- file.path(outdir,paste0(bedName, "_peaks.narrowPeak"))
  xlsFile <- file.path(outdir,paste0(bedName, "_peaks.xls"))

  #Read Summits!
  out <- data.table::fread(summitsFile, select = c(1,2,3,5))
  out <- GenomicRanges::GRanges(out$V1, IRanges::IRanges(out$V2 + 1, out$V3), score = out$V5)
  if(cleanUp){
    #Remove Files
    r2 <- suppressWarnings(file.remove(summitsFile, narrowPeaksFile, xlsFile))
  }
  if(save){
    outFile=file.path(outdir,paste0(prefix,"_summits.rds"))
    saveRDS(out,outFile)
  }
  return(out)

}

.getQuantiles<-function (v = NULL, len = length(v))
{
  if (length(v) < len) {
    v2 <- rep(0, len)
    v2[seq_along(v)] <- v
  }
  else {
    v2 <- v
  }
  p <- trunc(rank(v2))/length(v2)
  if (length(v) < len) {
    p <- p[seq_along(v)]
  }
  return(p)
}

nonOverlappingGR<-function (gr = NULL, by = "score", decreasing = TRUE, verbose = FALSE)
{
  require("GenomicRanges")
  require("BSgenome")
  stopifnot(by %in% colnames(mcols(gr)))
  .clusterGRanges <- function(gr = NULL, filter = TRUE, by = "score",
                              decreasing = TRUE) {
    gr <- sort(sortSeqlevels(gr))
    r <- GenomicRanges::reduce(gr, min.gapwidth = 0L, ignore.strand = TRUE)
    o <- findOverlaps(gr, r, ignore.strand = TRUE)
    mcols(gr)$cluster <- subjectHits(o)
    gr <- gr[order(mcols(gr)[, by], decreasing = decreasing),
    ]
    gr <- gr[!duplicated(mcols(gr)$cluster), ]
    gr <- sort(sortSeqlevels(gr))
    mcols(gr)$cluster <- NULL
    return(gr)
  }
  if (verbose) {
    message("Converging", appendLF = FALSE)
  }
  i <- 0
  grConverge <- gr
  while (length(grConverge) > 0) {
    if (verbose) {
      message(".", appendLF = FALSE)
    }
    i <- i + 1
    grSelect <- .clusterGRanges(gr = grConverge, filter = TRUE,
                                by = by, decreasing = decreasing)
    grConverge <- subsetByOverlaps(grConverge, grSelect,
                                   invert = TRUE, ignore.strand = TRUE)
    if (i == 1) {
      grAll <- grSelect
    }
    else {
      grAll <- c(grAll, grSelect)
    }
  }
  message(sprintf("Converged after %s iterations!", i))
  if (verbose) {
    message("\nSelected ", length(grAll), " from ", length(gr))
  }
  grAll <- sort(sortSeqlevels(grAll))
  return(grAll)
}

#' identify Reproducible Peaks(pass)
#' @param summitFiles summitFiles after macs2 callpeak
#' @param summitNames relative names for summitFiles
#' @param reproducibility reproducibility value
#' @param extendSummits extendSummits number
#' @param blacklist blacklist Granges object
#' @export
identifyReproduciblePeaks <- function(
  summitFiles = NULL,
  summitNames = NULL,
  reproducibility = 0.51,
  extendSummits = 250,
  blacklist = NULL,
  prefix = NULL
){
  require("GenomicRanges")
  require("BSgenome")
  require("dplyr")

  nonOverlapPassES <- tryCatch({

    summits <- lapply(seq_along(summitFiles), function(x){
      grx <- readRDS(summitFiles[x])
      grx <- subsetByOverlaps(grx, blacklist, invert = TRUE) #Not Overlapping Blacklist!
      grx$GroupReplicate <- paste0(summitNames[x])
      grx
    })
    summits <- Reduce("c", as(summits, "GRangesList"))

    extendedSummits <- resize(summits, extendSummits * 2 + 1, "center")
    extendedSummits <- lapply(split(extendedSummits, extendedSummits$GroupReplicate), function(x){
      nonES <- nonOverlappingGR(x, by = "score", decreasing = TRUE)
      nonES$replicateScoreQuantile <- round(.getQuantiles(nonES$score),3)
      nonES
    })
    extendedSummits <- Reduce("c", as(extendedSummits, "GRangesList"))

    nonOverlapES <- nonOverlappingGR(extendedSummits, by = "replicateScoreQuantile", decreasing = TRUE)

    overlapMat <- lapply(split(extendedSummits, extendedSummits$GroupReplicate), function(x){
      overlapsAny(nonOverlapES, x)
    }) %>% Reduce("cbind", .)

    if(length(summitFiles) > 1){
      nonOverlapES$Reproducibility <- rowSums(overlapMat)
      nonOverlapES$ReproducibilityPercent <- round(rowSums(overlapMat) / ncol(overlapMat) , 3)
      n <- length(summitFiles)
      minRep <- eval(parse(text=reproducibility))
      if(!is.numeric(minRep)){
        stop("Error reproducibility not numeric when evaluated!")
      }
      idxPass <- which(nonOverlapES$Reproducibility >= minRep)
      nonOverlapPassES <- nonOverlapES[idxPass]
    }else{
      nonOverlapES$Reproducibility <- rep(NA, length(nonOverlapES))
      nonOverlapPassES <- nonOverlapES
    }

    nonOverlapPassES$groupScoreQuantile <- round(.getQuantiles(nonOverlapPassES$replicateScoreQuantile),3)
    mcols(nonOverlapPassES) <- mcols(nonOverlapPassES)[,c("score","replicateScoreQuantile", "groupScoreQuantile", "Reproducibility", "GroupReplicate")]

    nonOverlapPassES

  }, error = function(e){
    message("INFO : identifyReproduciblePeaks Error!")

  })

  return(nonOverlapPassES)

}

#' check valid Granges
#' @export
.validGRanges<-function (gr = NULL)
{
  stopifnot(!is.null(gr))
  if (inherits(gr, "GRanges")) {
    return(gr)
  }
  else {
    stop("Error cannot validate genomic range!")
  }
}

#' check path
#' @export
.checkPath<-function (u = NULL, path = NULL, throwError = TRUE)
{
  require("dplyr")
  if (is.null(u)) {
    out <- TRUE
  }
  out <- lapply(u, function(x, error = TRUE) {
    if (Sys.which(x) == "") {
      if (!is.null(path) && file.exists(file.path(path,
                                                  x))) {
        o <- TRUE
      }
      else {
        if (throwError) {
          stop(x, " not found in path, please add ",
               x, " to path!")
        }
        else {
          o <- FALSE
        }
      }
    }
    else {
      o <- TRUE
    }
    return(o)
  }) %>% unlist %>% all
  return(out)
}

#' check valid BSgenome
#' @export
.validBSgenome<-function (genome = NULL, masked = FALSE)
{
  stopifnot(!is.null(genome))
  if (inherits(genome, "BSgenome")) {
    return(genome)
  }
  else if (is.character(genome)) {
    genome <- tryCatch({
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

#' extend Granges object
#' @param gr Granges object
#' @param upstream upstream value
#' @param downstream downstream value
extendGR<-function (gr = NULL, upstream = 50000, downstream = 50000)
{
  require("GenomicRanges")
  isMinus <- BiocGenerics::which(strand(gr) == "-")
  isOther <- BiocGenerics::which(strand(gr) != "-")
  start(gr)[isOther] <- start(gr)[isOther] - upstream
  end(gr)[isOther] <- end(gr)[isOther] + downstream
  end(gr)[isMinus] <- end(gr)[isMinus] + upstream
  start(gr)[isMinus] <- start(gr)[isMinus] - downstream
  return(gr)
}

#' annotate peak Grange Object
#' @param peaks peak Grange object
#' @param BSgenome BSgenome object
#' @param geneAnnotation geneAnnotation object
#' @export
fastAnnoPeaks <- function(
  peaks = NULL,
  BSgenome = NULL,
  geneAnnotation = NULL,
  promoterRegion = c(2000, 100)
){

  require("GenomicRanges")
  #Validate
  peaks <- .validGRanges(peaks)
  peakSummits <- resize(peaks,1,"center")
  geneAnnotation$genes <- .validGRanges(geneAnnotation$genes)
  geneAnnotation$exons <- .validGRanges(geneAnnotation$exons)
  geneAnnotation$TSS <- .validGRanges(geneAnnotation$TSS)
  BSgenome <- .validBSgenome(BSgenome)

  #First Lets Get Distance to Nearest Gene Start
  distPeaks <- distanceToNearest(peakSummits, resize(geneAnnotation$genes, 1, "start"), ignore.strand = TRUE)
  mcols(peaks)$distToGeneStart <- mcols(distPeaks)$distance
  mcols(peaks)$nearestGene <- mcols(geneAnnotation$genes)$symbol[subjectHits(distPeaks)]
  promoters <- extendGR(resize(geneAnnotation$genes, 1, "start"), upstream = promoterRegion[1], downstream = promoterRegion[2])
  op <- overlapsAny(peakSummits, promoters, ignore.strand = TRUE)
  og <- overlapsAny(peakSummits, geneAnnotation$genes, ignore.strand = TRUE)
  oe <- overlapsAny(peakSummits, geneAnnotation$exons, ignore.strand = TRUE)
  type <- rep("Distal", length(peaks))
  type[which(og & oe)] <- "Exonic"
  type[which(og & !oe)] <- "Intronic"
  type[which(op)] <- "Promoter"
  mcols(peaks)$peakType <- type

  #First Lets Get Distance to Nearest TSS's
  distTSS <- distanceToNearest(peakSummits, resize(geneAnnotation$TSS, 1, "start"), ignore.strand = TRUE)
  mcols(peaks)$distToTSS <- mcols(distTSS)$distance
  if("symbol" %in% colnames(mcols(geneAnnotation$TSS))){
    mcols(peaks)$nearestTSS <- mcols(geneAnnotation$TSS)$symbol[subjectHits(distPeaks)]
  }else if("tx_name" %in% colnames(mcols(geneAnnotation$TSS))){
    mcols(peaks)$nearestTSS <- mcols(geneAnnotation$TSS)$tx_name[subjectHits(distPeaks)]
  }

  #Get NucleoTide Content
  nucFreq <- BSgenome::alphabetFrequency(getSeq(BSgenome, peaks))
  mcols(peaks)$GC <- round(rowSums(nucFreq[,c("G","C")]) / rowSums(nucFreq),4)
  mcols(peaks)$N <- round(nucFreq[,c("N")] / rowSums(nucFreq),4)
  peaks

}

#' annotate bed or bam files(pass)
#' @param Files vector of bedfile
#' @param Names relative bedfile names
#' @param blacklist blacklist GRanges object
#' @param format input data format
#' @param genome genome name
#' @param outdir outdir to save result
#' @export
annoteFile<-function(Files=NULL,
                    Names=NULL,
                    blacklist=NULL,
                    format="BED",
                    genome="hg38",
                    outdir="annotation"
){

  stopifnot(inherits(blacklist,"GRanges"))
  stopifnot(tolower(genome)%in%c("hg38","mm10"))
  if(genome=="hg38"){
    require(BSgenome.Hsapiens.UCSC.hg38)
    bsgenome=BSgenome.Hsapiens.UCSC.hg38
    genomeSize=2.7e9
    keepChr=paste0("chr",c(1:22,"X"))
  }else{
    bsgenome=BSgenome.Mmusculus.UCSC.mm10
    genomeSize=1.87e9
    keepChr=paste0("chr",c(1:19,"X"))
  }
  if(is.null(blacklist)){
    blacklist<-.getBlacklist(genome)
  }
  message("INFO : create Gene Annotation")
  geneAnnotation=createGeneAnnotation(genome)

  message("INFO :CallPeaks ..")
  summitFiles<-c()
  summitNames<-c()

  if(!dir.create(outdir)){
    dir.create(outdir,recursive = TRUE,showWarnings = TRUE)
  }

  for(i in seq_along(Files)){
      name=Names[i]
      bedfile=Files[i]
      cat(sprintf("INFO :   %d of %d  ---  %s \n",i,length(Files),name))
      peak=callSummitsMACS2(File=bedfile,
                            format=format,
                            genomeSize=genomeSize,
                            outdir=outdir,
                            prefix=name,
                            cleanUp=FALSE,
                            save=TRUE)
      peak = keepGR(peak,genome)
      keepChrPeak=file.path(outdir,paste0(name,"_final_summits.rds"))
      cat(sprintf("INFO : save into : %s\n",keepChrPeak))
      saveRDS(peak,keepChrPeak)
      summitFiles=c(summitFiles,keepChrPeak)
      summitNames=c(summitNames,name)
  }

  gr<-tryCatch({
      message("INFO : identify ReproduciblePeaks...")
      identifyReproduciblePeaks(summitFiles=summitFiles,
                               summitNames=summitNames,
                               blacklist=blacklist,
                               reproducibility = 0.51,
                               extendSummits = 250)
    },error=function(e){
      message("IdentifyReproduciblePeaks Error!")
      NULL
    })
  gr<-tryCatch({
    message("INFO : anno Peaks ...")
    fastAnnoPeaks(gr,bsgenome,geneAnnotation)
  },error=function(e){
    message("Annotate Peaks Error!")
    NULL
  })
  if(!is.null(gr)){
    saveRDS(gr,file.path(outdir,"AnnoPeaks_summits.rds"))
  }
  message("INFO : Done!")
}
