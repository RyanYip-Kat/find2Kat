library(data.table)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)

#' filter chroms
keepFilteredChromosomes <- function(x, remove = c("chrM"), underscore = TRUE, standard = TRUE, pruning.mode="coarse"){
  #first we remove all non standard chromosomes
  if(standard){
    x <- GenomeInfoDb::keepStandardChromosomes(x, pruning.mode = pruning.mode)
  }
  #Then check for underscores or specified remove
  seqNames <- GenomeInfoDb::seqlevels(x)
  chrRemove <- c()
  #first we remove all chr with an underscore
  if(underscore){
    chrRemove <- c(chrRemove, which(grepl("_", seqNames)))
  }
  #next we remove all chr specified in remove
  chrRemove <- c(chrRemove, which(seqNames %in% remove))
  if(length(chrRemove) > 0){
    chrKeep <- seqNames[-chrRemove]
  }else{
    chrKeep <- seqNames
  }
  #this function restores seqlevels
  GenomeInfoDb::seqlevels(x, pruning.mode=pruning.mode) <- chrKeep
  return(x)
}


#' reading GTF file
#' @param file gtf file
#' @export
getGeneGTF <- function(file){
  #Import
  message("Reading in GTF...")
  importGTF <- rtracklayer::import(file)
  #Exon Info
  message("Computing Effective Exon Lengths...")
  exonGTF <- importGTF[importGTF$type=="exon",]
  exonList <- reduce(split(exonGTF, mcols(exonGTF)$gene_id))
  exonReduced <- unlist(exonList, use.names=TRUE)

  mcols(exonReduced)$gene_id <- names(exonReduced)
  mcols(exonReduced)$widths <- width(exonReduced)
  exonSplit <- split(exonReduced$widths, mcols(exonReduced)$gene_id)
  exonLengths <- lapply(seq_along(exonSplit), function(x) sum(exonSplit[[x]])) %>%
    unlist %>% data.frame(row.names=names(exonSplit), effLength=.)
  #Gene Info
  message("Constructing gene GTF...")
  geneGTF1 <- importGTF[importGTF$type=="gene",]
  geneGTF2 <- GRanges(
    seqnames=paste0("chr",seqnames(geneGTF1)),
    ranges=ranges(geneGTF1),
    strand=strand(geneGTF1),
    gene_name=geneGTF1$gene_name,
    gene_id=geneGTF1$gene_id
  ) %>% keepFilteredChromosomes %>% sortSeqlevels %>% sort(.,ignore.strand=TRUE)
  mcols(geneGTF2)$exonLength <- exonLengths[geneGTF2$gene_id,]
  return(geneGTF2)
}
