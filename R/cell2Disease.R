#' Granges to Features
#' @param gr Granges object
#' @export
grToFeature <- function(gr){
  require("GenomicRanges")
  #.validInput(input = gr, name = "gr", valid = c("GRanges"))
  peakinfo <- data.frame(
    row.names = paste(seqnames(gr),start(gr),end(gr),sep="_"),
    site_name = paste(seqnames(gr),start(gr),end(gr),sep="_"),
    chr = gsub("chr","",as.character(seqnames(gr))),
    bp1 = start(gr),
    bp2 = end(gr)
  )
  return(peakinfo)
}

#' feature to Granges
#' @param feature feature name
#' @export
featureToGR <- function(feature){
  require("GenomicRanges")
  featureSplit <- stringr::str_split(paste0(feature), pattern = "_", n = 3, simplify = TRUE)
  gr <- GRanges(featureSplit[,1],IRanges(as.integer(featureSplit[,2]),as.integer(featureSplit[,3])))
  return(gr)
}

#' RangedSummarizedExperiment objects to newCellDataSet
#' @param se RangedSummarizedExperiment Object
#' @param binarize whether binarize matrix
#' @export
makeCDS <- function(se, binarize = TRUE){
  require("SummarizedExperiment")
  require("monocle")

  #.validInput(input=se,name="se",valid = c("RangedSummarizedExperiment"))
  peakinfo <- grToFeature(se)
  mat <- assay(se)
  if(binarize){
    mat@x[which(mat@x > 0)] <- 1
  }
  cellinfo <- data.frame(colData(se))
  cellinfo$cells <- rownames(cellinfo)
  cds <-  suppressWarnings(newCellDataSet(mat,
                                          phenoData = methods::new("AnnotatedDataFrame", data = cellinfo),
                                          featureData = methods::new("AnnotatedDataFrame", data = peakinfo),
                                          expressionFamily=negbinomial.size(),
                                          lowerDetectionLimit=0))
  fData(cds)$chr <- as.character(fData(cds)$chr)
  fData(cds)$bp1 <- as.numeric(as.character(fData(cds)$bp1))
  fData(cds)$bp2 <- as.numeric(as.character(fData(cds)$bp2))
  cds <- cds[order(fData(cds)$chr, fData(cds)$bp1),]
  return(cds)
}


#' get TxDb genes
#' @param txdb txdb oject
#' @export
getTxDbGenes <- function(txdb = NULL, orgdb = NULL, gr = NULL, ignore.strand = TRUE){
  require("GenomicRanges")
  require("AnnotationDbi")
  if (is.null(genome)) {
    if (is.null(txdb) | is.null(orgdb)) {
      stop("If no provided genome then you need txdb and orgdb!")
    }
  }

  if (is.null(gr)) {
    genes <- GenomicFeatures::genes(txdb)
  }else {
    genes <- suppressWarnings(subsetByOverlaps(GenomicFeatures::genes(txdb), gr, ignore.strand = ignore.strand))
  }

  if (length(genes) > 1) {
    mcols(genes)$symbol <- suppressMessages(mapIds(orgdb,
                                                   keys = mcols(genes)$gene_id, column = "SYMBOL", keytype = "ENTREZID",
                                                   multiVals = "first"))
    genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)
    names(genes) <- NULL
    out <- genes
  }else {
    out <- GRanges(seqnames(gr), ranges = IRanges(0, 0), gene_id = 0, symbol = "none")[-1]
  }

  return(out)

}

######## peakset to cicero object(pass)
#' SummarizedExperiment to Cicero object
#' @param obj RangedSummarizedExperiment Object
#' @param umapMat umap Embedding matrix or data.frame
#' @param genome genome name(hg38 or mm10)
#' @param binarize whether binarize data
#' @export
seTocicero<-function(obj=NULL,
                     umapMat=NULL,
                     flank=250*10^3,
                     #genome="hg38",
                     binarize=TRUE,
                     save=FALSE,
                     outdir=NULL){
  require("cicero")
  require("monocle")
  require("Matrix")
  require("GenomicRanges")

  ###############
  stopifnot(inherits(obj,"RangedSummarizedExperiment"))
  stopifnot(class(umapMat)%in%c("matrix","data.frame"))

  peaks=rowRanges(obj)
  rownames(obj)=paste(seqnames(peaks),start(peaks),end(peaks),sep="_")

  if(inherits(umapMat,"matrix")){
    mat<-as.data.frame(umapMat)
  }else{
    mat<-umapMat
  }
  colnames(mat)=c("UMAP1","UMAP2")
  colData(obj)$UMAP1=mat$UMAP1
  colData(obj)$UMAP2=mat$UMAP2
  dimred <- data.frame(
    row.names = colnames(obj),
    colData(obj)$UMAP1,
    colData(obj)$UMAP2
  )

  obj2<-obj

  #Get CDS
  message("INFO : Get CDS")
  obj <- makeCDS(obj, binarize = binarize)
  obj <- detectGenes(obj)
  message("INFO : estimateSizeFactors")
  obj <- estimateSizeFactors(obj)
  message("INFO : make_cicero_cds ")
  ciceroObj <- make_cicero_cds(obj, k = 50, reduced_coordinates = dimred[colnames(obj),])

  #Compute Correlations
  message("INFO : Computing grouped correlations...")
  gr <- featureToGR(featureData(ciceroObj)[[1]])
  o <- suppressWarnings(as.matrix( findOverlaps(resize( resize(gr,1,"center"), 2*flank + 1, "center"), resize(gr,1,"center"), ignore.strand=TRUE) ))
  o <- data.table::as.data.table(data.frame(i = matrixStats::rowMins(o), j = matrixStats::rowMaxs(o)))
  o <- data.frame(o[!duplicated(o),])
  o <- o[o[,1]!=o[,2],]
  o$cor <- rowCorCpp(o[,1], o[,2], assayData(ciceroObj)$exprs, assayData(ciceroObj)$exprs)
  connections <- data.frame(
    Peak1 = featureData(ciceroObj)[[1]][o[,1]],
    Peak2 = featureData(ciceroObj)[[1]][o[,2]],
    coaccess = o[,3]
  )
  if(save){
    if(is.null(outdir)){
      outdir="connection"
    }
    if(!dir.exists(outdir)){
      dir.create(outdir,recursive = TRUE,showWarnings = TRUE)
    }
    #Save Output
    message("INFO : Save Output..")
    saveRDS(ciceroObj,file.path(outdir,"ciceroObj.rds"))
    saveRDS(connections, file.path(outdir,"Peaks-Co-Accessibility.rds"))
    saveRDS(obj2,file.path(outdir,"scATAC-Summarized-Experiment.rds"))
  }
  return(list("conn"=connections,"obj"=obj2))
}

#' connect SNP to scATAC Object
#' @param se RangedSummarizedExperiment object(eg,from function seTocicero)
#' @param snp SNP GRanges object
#' @param conn connections object from cicero computation
#' @param cutoff float value,filter connection
#' @param save whether save result
#' @export
gwasSNPconn<-function(se=NULL,
                      snp=NULL,
                      conn=NULL,
                      cutoff=0.35,
                      genome="hg38",
                      outdir=NULL,
                      save=TRUE){
  require("SummarizedExperiment")
  require("GenomicRanges")
  require("chromVAR")
  require("dplyr")

  stopifnot(inherits(se,"RangedSummarizedExperiment"))
  stopifnot(inherits(snp,"GRanges"))
  stopifnot("Disease"%in%colnames(mcols(snp)))
  stopifnot(tolower(genome)%in%c("hg38","mm10"))

  gr <- snp
  conn<-na.omit(conn)
  conn <- conn[conn[,3] >= cutoff,]
  peaknames <- paste(seqnames(se),start(se),end(se),sep="_")
  conn[,4] <- match(paste0(conn[,1]), peaknames)
  conn[,5] <- match(paste0(conn[,2]), peaknames)
  connMat <- Matrix::sparseMatrix(i=conn[,4],j=conn[,5],x=rep(1,nrow(conn)),dims=c(nrow(se),nrow(se)))

  if(tolower(genome)=="hg38"){
    bsgenome<-BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  }else{
    bsgenome<-BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  }
  #Add Bias
  message("INFO : Add Bias")
  se <- addGCBias(se, genome =bsgenome)

  #Overlap GWAS SNPs
  message("INFO Overlap GWAS SNPs")
  o <- lapply(split(gr, gr$Disease), function(x){
    #extend snp +- 500 bp
    overlapsAny(se, resize(x,1001,"center"), ignore.strand = TRUE)
  }) %>% Reduce("cbind",.)
  ow <- which(o > 0, arr.ind=TRUE)
  matches <- Matrix::sparseMatrix(i=ow[,1],j=ow[,2],x=o[cbind(ow[,1],ow[,2])], dims = c(nrow(se),length(unique(gr$Disease))))
  colnames(matches) <- names(split(gr, gr$Disease))

  #Use connections mat!
  message("INFO : Use connections match")
  matches2 <- matches
  idx <- which(Matrix::rowSums(matches2)>0) #which peaks have a snp
  for(i in seq_along(idx)){
    if(i %% 100 == 0){message(sprintf("%s of %s",i,length(idx)))}
    #peaks co-accessible to peak with snp
    coi <- unique(c(which(connMat[,idx[i]]>0),which(connMat[idx[i],]>0)))
    if(length(coi) > 0){
      #create sub mat
      mati <- as(t(replicate(length(coi), matches[idx[i],])),"dgCMatrix")
      #add it to sub mat of connected peaks
      matches2[coi,,drop=FALSE] <- matches2[coi,,drop=FALSE] + mati
    }
  }
  diff <- Matrix::colSums(matches2) - Matrix::colSums(matches)
  #print(diff)

  #Make Annotation SE
  message("INFO : Make Annotation SummarizedExperiment")
  anno_ix <- SummarizedExperiment(assays=SimpleList(motifMatches=matches2), rowRanges=rowRanges(se))

  #Compute Deviations
  message("INFO : Compute Deviations")
  assayNames(se)="counts"
  dev <- computeDeviations(se, anno_ix)

  #compute variability
  message("INFO : compute variability")
  metadata(dev)$Variability <- computeVariability(dev)

  #add matches
  metadata(dev)$gwasMatches <- anno_ix

  if(save){
    if(is.null(outdir)){
      outdir="connection"
    }
    if(!dir.exists(outdir)){
      dir.create(outdir,recursive = TRUE,showWarnings = TRUE)
    }
    #save output
    message("INFO : Save Result")
    saveRDS(dev, file.path(outdir,"GWAS-Co-Accessibility-chromVAR-Summarized-Experiment.rds"))
  }
  return(dev)
}

#' plot Deviations
#' @param dev Deviations object
#' @param groupBy column in metadata as group,and plot with Diease
#' @param plot whether plot
#' @export
DevHeatmap<-function(dev=NULL,
                     groupBy=NULL,
                     plot=TRUE){
  require(pheatmap)
  require(chromVAR)
  require(SummarizedExperiment)
  require(magrittr)

  dev_score=deviationScores(dev)
  metadata=as.data.frame(colData(dev))
  Clusters=metadata[[groupBy]]

  score_name=rownames(dev_score)
  score_list=list()
  for(name in score_name){
    cat(sprintf("INOF : get :%s median deviation score\n",name))
    score=dev_score[name,]
    X=data.frame("Score"=score,"Cluster"=Clusters)
    z=as.data.frame(X %>% group_by(Cluster) %>% summarise(v = median(Score,na.rm=TRUE)))
    rownames(z)=z$Cluster
    z=subset(z,select="v")
    colnames(z)=name
    score_list[[name]]=z
  }

  scores=do.call(cbind,score_list)
  scores[is.na(scores)]=0.0
  if(plot){
    pdf(paste0("Median-DeviationScores-Across-",groupBy,".pdf"),width=12,height=16)
    pheatmap(t(scores),
             scale="column",
             cluster_cols=TRUE,
             cluster_rows = TRUE,
             fontsize = 10,
             fontsize_row = 10,
             fontsize_col=10
    )
    dev.off()
  }
  return(scores)
}

#' main pipeline
#' @param obj RangedSummarizedExperiment object(eg,from function seTocicero)
#' @param snp SNP GRanges object
#' @param umapMat umap Embedding matrix or data.frame
#' @param genome genome name(hg38 or mm10)
#' @param binarize whether binarize data
#' @param cutoff float value,filter connection
#' @param save whether save result
#' @export
chromVARSNPpipeline<-function(obj=NULL,
                              umapMat=NULL,
                              SNP=NULL,
                              binarize=TRUE,
                              groupBy=NULL,
                              cutoff=0.35,
                              genome="hg38",
                              outdir=NULL,
                              save=FALSE){

  message("step 1..")
  res<-seTocicero(obj=obj,
                  umapMat=umapMat,
                  binarize=binarize,
                  save=save,
                  outdir=outdir)
  conn<-res[["conn"]]
  obj2<-res[["obj"]]

  message("step 2..")
  dev<-gwasSNPconn(se=obj2,
                   conn=conn,
                   snp=SNP,
                   genome = genome,
                   cutoff=cutoff,
                   save = save,
                   outdir=outdir
                   )
  message("step 3..")
  score<-DevHeatmap(dev=dev,
                    groupBy = groupBy,
                    plot=TRUE
                    )
  return(score)
}
