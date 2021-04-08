#' apply ChipSeeker annote for summits result from MACS2
#' @param summitFiles summitFiles rds file from CallSummitsMACS
#' @param summitNames relativte name for summitFiles
#' @param filterChr whether file GRanges object
#' @param genome genome name,hg38 or mm10
#' @param plot whether plot
#' @param outdir output path to save result
#' @return peakAnnoList object
#' @export
summitsChipSeekerAnnote<-function(summitFiles=NULL,
                       summitNames=NULL,
                       filterChr=TRUE,
                       genome="hg38",
                       outdir="chipSeeker",
                       plot=FALSE){
  require(ggplot2)
  stopifnot(length(summitFiles)==length(summitNames))
  if(tolower(genome)=="hg38"){
    require("TxDb.Hsapiens.UCSC.hg38.knownGene")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  }else if(tolower(genome)=="mm10"){
    require("TxDb.Mmusculus.UCSC.mm10.knownGene")
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  }else{
    stop("Invalid genome!!!")
  }

  if(is.null(outdir)){
      outdir<-"result"
  }
  makedir(outdir)

  message("loading call Summits MACS2 result..")
  peakList<-list()
  for(i in seq_along(summitFiles)){
    Name<-summitNames[i]
    File<-summitFiles[i]

    cat(sprintf("Reading [ No.%d/%d ] -- %s -- callpeak from -- %s --\n",i,length(summitFiles),Name,File))
    gr<-readRDS(File)
    if(filterChr){
      gr<-keepGR(gr,genome=genome)
    }
    peakList[[Name]]<-gr
    if(plot){
      cat(sprintf("Covplot [ No.%d/%d ] -- %s\n",i,length(summitFiles),Name))
      p<-covplot(gr, weightCol="score")
      filename<-file.path(outdir,paste0(Name,"_covplot.pdf"))
      ggsave(filename=filename,plot=p,width = 12,height = 16)
      }
  }

  message("chipSeeker annotation..")
  peakAnnoList <- lapply(peakList, annotatePeak, TxDb=txdb,
                         tssRegion=c(-3000, 3000), verbose=TRUE)
  if(plot){
    message("---- plot list ----")
    message("AnnoBar plot..")
    p<-plotAnnoBar(peakAnnoList)
    filename<-file.path(outdir,"AnnoBar.pdf")
    ggsave(filename=filename,plot=p,width = 16,height = 12)

    message("DistToTSS plot..")
    p<-plotDistToTSS(peakAnnoList)
    filename<-file.path(outdir,"DistToTSS.pdf")
    ggsave(filename=filename,plot=p,width = 16,height = 12)

    message("vennplot : Overlap of peaks and annotated genes..")
    genes <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
    pdf(file = file.path(outdir,"vennplot.pdf"),width = 8,height = 8)
    print(vennplot(genes))
    dev.off()

    message("upsetplot ..")
    for(Name in names(peakAnnoList)){
      anno<-peakAnnoList[[Name]]
      p<-upsetplot(anno, vennpie=TRUE)
      filename<-file.path(outdir,paste0(Name,"_upsetplot.pdf"))
      ggsave(filename=filename,plot=p,width = 16,height = 12)
    }
  }
  return(peakAnnoList)
}

#' GO Pathway analysis for peakAnnoList from ChipSeeker Annotation
#' @param peakAnnoList peakAnnoList object from annotatePeak
#' @param fun fun: One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway"
#' @param genome genome name,hg38 or mm10
#' @param plot whether plot
#' @param outdir output path to save result
#' @return compareCluster object
#' @export
chipseekerGO<-function(peakAnnoList=NULL,
                       fun="enrichKEGG",
                       genome="hg38",
                       plot=TRUE,

                       outdir="chipSeeker"){
  require(ggplot2)
  require(clusterProfiler)
  if(tolower(genome)=="hg38"){
    require("org.Hs.eg.db")
    OrgDb <- "org.Hs.eg.db"
  }else if(tolower(genome)=="mm10"){
    require("org.Mm.eg.db")
    OrgDb <- "org.Mm.eg.db"
  }else{
    stop("Invalid genome!!!")
  }


  message("compareCluster GO analysis..")
  genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
  names(genes) = sub("_", "\n", names(genes))
  comp <- compareCluster(geneCluster = genes,
                         fun = fun,
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         OrgDb=OrgDb)
  if(plot){
    if(is.null(outdir)){
      outdir<-"result"
    }
    makedir(outdir)
    pdf(file = file.path(outdir,"vennplot.pdf"),width = 12,height = 16)
    print(dotplot(comp, showCategory = 15, title = "KEGG Pathway Enrichment Analysis"))
    dev.off()
  }
  return(comp)
}
