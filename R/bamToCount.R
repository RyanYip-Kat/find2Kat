#' bam file to count matrix
#' @param bamfiles bame files after aligment from fastq file
#' @param bamNames relative bam file names
#' @param peaks GRanges object
#' @param addCol colname add into colData
#' @param paired whether paried-seq
#' @return ‘RangedSummarizedExperiment-class’ object
#' @export
bam2Count<-function(bamfiles=NULL,
                    bamNames=NULL,
                    peaks=NULL,
                    addCol="Sample",
                    paired =  TRUE){
  #if(length(peaks)!=length(bamfiles)){
  #  if(length(peaks)!=1){
  #    stop("Length of Peak files must equal to 1 or length of bamfiles!!!")
  # }
  #}
  stopifnot(length(bamfiles)==length(bamNames))

  metaData<-DataFrame("Sample"=bamNames)
  colnames(metaData)<-addCol
  fragment_counts <- getCounts(bamfiles,
                               peaks=peaks,
                               paired =  paired,
                               by_rg = TRUE,
                               format = "bam",
                               colData = metaData)
  return(fragment_counts)
}
