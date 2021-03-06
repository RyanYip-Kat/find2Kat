#' link peak to TF
#' @param diffRNA object from ArchR::getMarkerFeatures(ArchRProject,useMatrix="GeneScoreMatrix" or "GeneIntegraMatrix")
#' @param diffATAC object from ArchR::getMarkerFeatures(ArchRProject,useMatrix="PeakMatrix")
#' @param matches object from ArchR::getMatches
#' @param p2gLinks object from `getP2GLinks` function
#' @export
link2TF<-function(diffRNA=NULL,
                  diffATAC=NULL,
                  p2gLinks=NULL,
                  matches=NULL,
                  log2fc=0.1,
                  pval=0.05
){
  require("SummarizedExperiment")
  require("Matrix")

  #stopifnot(inherits(diffRNA,"SummarizedExperiment"))
  #stopifnot(class(diffATAC)=="RangedSummarizedExperiment")
  #stopifnot(class(p2gLinks)=="GRanges")
  #stopifnot(class(matches)=="RangedSummarizedExperiment")

  rownames(diffRNA)<-rowData(diffRNA)$name
  rownames(matches) <- paste(seqnames(matches),start(matches),end(matches), sep = "_")
  rownames(diffATAC) <- paste(seqnames(diffATAC),start(diffATAC),end(diffATAC), sep = "_")

  #Make P2G Mats
  message("Make P2G Mats..")
  names(p2gLinks) <- paste0("l",seq_along(p2gLinks))
  dAp2g <- diffATAC[featureName(p2gLinks),colnames(diffRNA)]
  dRp2g <- diffRNA[mcols(p2gLinks)$gene_name[mcols(p2gLinks)$gene_name%in%rownames(diffRNA)],]
  rownames(dAp2g) <- names(p2gLinks)
  rownames(dRp2g) <- names(p2gLinks)[mcols(p2gLinks)$gene_name%in%rownames(diffRNA)]
  dAp2g=dAp2g[rownames(dRp2g),]

  pdR=paste(seqnames(dAp2g),start(dAp2g),end(dAp2g), sep = "_")
  matches=matches[pdR,]

  #Identify Significant Peaks and Genes per MPAL Subpopulation
  message("Identify Significant Peaks and Genes per MPAL Subpopulation..")
  sigMatP2G <- abs(assays(dAp2g)[["Log2FC"]]) >= log2fc & assays(dAp2g)[["Pval"]] <= pval & abs(assays(dRp2g)[["Log2FC"]]) >= log2fc & assays(dRp2g)[["Pval"]] <= pval


  #Which have at least 1 MPAL subpopulation with a linked (within p2glinks) diff peak and diff gene
  message("Which have at least 1 MPAL subpopulation with a linked (within p2glinks) diff peak and diff gene..")
  sigP2G <- p2gLinks[which(Matrix::rowSums(sigMatP2G) > 0)]

  #List of TFs to identify Targerts
  tfs <-  colnames(matches)
  #Identify Targets
  message("Identify Targets..")
  t2gDF <- lapply(seq_along(tfs), function(x){
    tryCatch({
      #Subset Linked Peaks that contain TF Motif
      peaksx <- names(which(assay(matches[,tfs[x]])[,1]))
      #Subset sigMat by linked peaks with TF Motif
      sigMatP2Gx <- sigMatP2G[rownames(sigMatP2G) %in% peaksx,]

      #Subset links by linked peaks with TF Motif
      linksx <- sigP2G[featureName(sigP2G) %in% peaksx]

      #Figure out which samples are true for each link
      #sigMatP2Gx <- sigMatP2G[names(linksx),]

      #Compute sum of MPAL subpopulations that are diff peak and diff peak for every peak connected to the gene
      sigDFListx <- split(data.frame(sigMatP2Gx), mcols(linksx)$gene_name) %>% lapply(., Matrix::colSums)

      #Convert to Data Frame
      sigDFx <- data.frame(Reduce("rbind",sigDFListx))
      rownames(sigDFx) <- names(sigDFListx)

      #Set to maximum of 1 for each MPAL subpop for all diff peak gene combo
      sigDFx[sigDFx > 1] <- 1

      #Compute number of positive MPAL subpopulations
      nSubpop <- Matrix::rowSums(sigDFx)
      cat(sprintf("INFO : [ N = %d , P = %f ]\n",nSubpop,nSubpop / ncol(sigMatP2Gx)))
      #Summed R2 linkage for each linked positive diff peak gene combo (Max of Cor Healthy Cor Cancer then Squared for each Gene)
      maxCor <- pmax(mcols(linksx)$CorrelationHealthy,mcols(linksx)$CorrelationDisease,na.rm=TRUE)
      linkageScore <- split(maxCor^2, mcols(linksx)$gene_name) %>% lapply(., sum) %>% unlist(use.names=TRUE)
      #Return Summary Data Frame
      data.frame(TF = tfs[x], Gene = names(nSubpop), N = nSubpop, P = nSubpop / ncol(sigMatP2Gx), linkageScore = linkageScore[names(nSubpop)])
    },error=function(e){cat(sprintf("INFO : Invalid TF : [ %s ]\n",tfs[x]))})
  })
  t2gDF=do.call(rbind,t2gDF)
  return(t2gDF)
}

#' caculate P2G Link ,please reference to (https://github.com/GreenleafLab/MPAL-Single-Cell-2019)(pass)
#' @param seA seATAC RangedSummarizedExperiment Object,can get from ArchR::getMatrxFromProject (after use function gr2Feature)
#' @param seB seRNA  RangedSummarizedExperiment Object(after use function seRNA2Rse)
#' @param heathP2G can get from `metadata(ArchRProj@peakSet)$Peak2GeneLinks`
#' @param diseaseP2G same as heathP2G
#' @param totalP2G same as heathP2G
#' @param corCutOff Pearson Correlation Cutoff,default 0.2
#' @param fdrCutOff FDR Cutoff,default 0.1
#' @param distCutOff Min Dist to TSS,default 2500
#' @param gtf gene gtf file or gene Granges object(find2Kat::GeneGTF)
#' @export
getP2GLinks<-function(seA=NULL,
                      seB=NULL,
                      corCutOff=0.2, #Pearson Correlation Cutoff
                      fdrCutOff=0.1, #FDR Cutoff
                      distCutOff=2500, #Min Dist to TSS
                      heathP2G=NULL,
                      diseaseP2G=NULL,
                      totalP2G=NULL,
                      GTF="gene.gtf",
                      ...){

  #require("GenomicRanges")
  #require("SummarizedExperiment")
  if(class(seB)!="RangedSummarizedExperiment"){
    seB<-seRNA2Rse(seRNA=seB,gtf = GTF)
  }
  seA<-gr2Feature(seA)

  colnames(totalP2G)=c('idxTotalA','idxTotalB','CorrelationTotal','FDRTotal','VarTotalA','VarTotalB')
  colnames(diseaseP2G)=c('idxDiseaseA','idxDiseaseB','CorrelationDisease','FDRDisease','VarDiseaseA','VarDiseaseB')
  colnames(heathP2G)=c('idxHealthyA','idxHealthyB','CorrelationHealthy','FDRHealthy','VarHealthyA','VarHealthyB')
  o=cbind(totalP2G,diseaseP2G,heathP2G)

  fixA <- "center"
  fixB <- "start"
  associationWindow <- 2 * 250*10^3 + 1 #+-250 Kb

  if(class(GTF)=="GRanges"){
    gtf <- GTF
  }
  else if(class(GTF)=="character"){
    gtf <- getGeneGTF(GTF)
  }else{
    stop("Invalid gtf!!!")
  }
  tssRNA <- resize(gtf, 1, "start")
  strand(tssRNA) <- "*"
  peakLinks <- rowRanges(seA)[o[,1]]
  geneLinks <- rowRanges(seB) %>% resize(1, "start") %>% {.[o[,2]]}
  mcolsLinks <- data.frame(geneLinks)[,c("seqnames","start","strand","gene_name","gene_id","exonLength")]
  colnames(mcolsLinks) <- c("gene_chr","gene_start","gene_strand","gene_name","gene_id","exonLength")
  mcolsLinks <- cbind(mcolsLinks, data.frame(o))
  mcolsLinks$nearestGene <- tssRNA$gene_name[subjectHits(distanceToNearest(peakLinks, tssRNA, ignore.strand=TRUE))]
  mcolsLinks$nearestGeneDistance <- mcols(distanceToNearest(peakLinks, tssRNA, ignore.strand=TRUE))$distance
  mcols(peakLinks) <- mcolsLinks
  peakLinks$peakName <- paste(seqnames(peakLinks), start(peakLinks), end(peakLinks), sep = "_")

  peakLinks$sigDisease <- peakLinks$CorrelationDisease >= corCutOff & peakLinks$FDRDisease <= fdrCutOff
  peakLinks$sigHealthy <- peakLinks$CorrelationHealthy >= corCutOff & peakLinks$FDRHealthy <= fdrCutOff


  linksSig <- peakLinks[which(peakLinks$sigDisease | peakLinks$sigHealthy)]
  linksDisease <- peakLinks[which(peakLinks$sigDisease & !peakLinks$sigHealthy)]
  linksHealthy <- peakLinks[which(!peakLinks$sigDisease & peakLinks$sigHealthy)]
  linksShared <- peakLinks[which(peakLinks$sigDisease & peakLinks$sigHealthy)]

  outMatch <- list(
    seA = seA[unique(mcols(linksSig)$peakName),],
    seB = seB[unique(mcols(linksSig)$gene_name),],
    linksDisease = linksDisease,
    linksHealthy = linksHealthy,
    linksShared = linksShared,
    linksSig = linksSig,
    linksAll = peakLinks

  )
  return(outMatch)
}
