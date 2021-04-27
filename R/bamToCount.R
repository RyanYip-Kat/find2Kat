#' bam file to count matrix
#' @param files bame files after aligment from fastq file
#' @param Names relative bam file names
#' @param peaks GRanges object
#' @param addCol colname add into colData
#' @param paired whether paried-seq
#' @param by_rg use RG tags in bam to separate groups
#' @param format bed or bam format
#' @return ‘RangedSummarizedExperiment-class’ object
#' @export
bam2Count<-function(files=NULL,
                    Names=NULL,
                    peaks=NULL,
                    format="bam",
                    addCol="Sample",
                    paired =  TRUE,
                    by_rg=TRUE){
  #if(length(peaks)!=length(bamfiles)){
  #  if(length(peaks)!=1){
  #    stop("Length of Peak files must equal to 1 or length of bamfiles!!!")
  # }
  #}
  stopifnot(length(files)==length(Names))

  metaData<-DataFrame("Sample"=Names)
  colnames(metaData)<-addCol
  fragment_counts <- getCounts(files,
                               peaks=peaks,
                               paired =  paired,
                               by_rg = by_rg,
                               format = format,
                               colData = metaData)
  return(fragment_counts)
}
