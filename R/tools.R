#######################

#' export GRanges peakset object into BW file
#'
#' @param gr GRanges Object
#' @export
#'
PeakSet2BW<-function(gr=NULL,outdir="bwFile"){
  require("GenomicRanges")
  require("rtracklayer")

  stopifnot(inherits(gr,"GRanges"))

  if(!dir.exists(outdir)){
    dir.create(outdir,recursive = TRUE,showWarnings = TRUE)
  }
  message("INFO : Export Into BW File")
  filename=file.path(outtir,"PeakSet.bw")
  rtracklayer::export(gr,con=filename,format="BigWig")
  message("INFO : Done!")
}


featureName <- function(gr){
  paste(GenomicRanges::seqnames(gr),GenomicRanges::start(gr),GenomicRanges::end(gr),sep="_")
}

#' softmax
#' @export
softmax=function(V){
  expV=exp(V)
  return(expV/sum(expV))
}

#' make dir
#' @export
makedir<-function(path){
  if(!dir.exists(path)){
    dir.create(path,recursive=TRUE)
  }
}

#' Max Min Scale
#' @export
MinMaxScale<-function(x){
  stopifnot(is.vector(x))
  v=(x-min(x))/(max(x)-min(x)+1e-4)
  return(v)
}

#' rowZscores
#' @export
rowZscores<-function (m = NULL, min = -2, max = 2, limit = FALSE)
{
  z <- sweep(m - Matrix::rowMeans(m), 1, matrixStats::rowSds(m), `/`)
  if (limit) {
    z[z > max] <- max
    z[z < min] <- min
  }
  return(z)
}


requirePackage<-function (x = NULL, load = TRUE, installInfo = NULL, source = NULL)
{
  if (x %in% rownames(installed.packages())) {
    if (load) {
      suppressPackageStartupMessages(require(x, character.only = TRUE))
    }
    else {
      return(0)
    }
  }
  else {
    if (!is.null(source) & is.null(installInfo)) {
      if (tolower(source) == "cran") {
        installInfo <- paste0("install.packages(\"",
                              x, "\")")
      }
      else if (tolower(source) == "bioc") {
        installInfo <- paste0("BiocManager::install(\"",
                              x, "\")")
      }
      else {
        stop("Unrecognized package source, available are cran/bioc!")
      }
    }
    if (!is.null(installInfo)) {
      stop(paste0("Required package : ", x, " is not installed/found!\n  Package Can Be Installed : ",
                  installInfo))
    }
    else {
      stop(paste0("Required package : ", x, " is not installed/found!"))
    }
  }
}

#' function for matrix to doheatmap in archr Style
#' @param seurat An `Seurat` object.
#' @param mat an matrix as sparse matrix
#' @param scale bool value,wether to do scale in heatmap.
#' @param limits limit the min and max value in mat.
#' @param colData provide cell bar in top heatmap if it is not null.
#' @param clusterRows bool value,whether do cluster in row.
#' @param clusterCols same as clusterRows.
#' @param customRowLabel custom label in row want to show.
#' ...
#' @export
ArchRHeatmap <- function(
  mat = NULL,
  scale = FALSE,
  limits = c(min(mat), max(mat)),
  colData = NULL,
  color = paletteContinuous(set = "solarExtra", n = 100),
  clusterCols = TRUE,
  clusterRows = FALSE,
  labelCols = FALSE,
  labelRows = FALSE,
  colorMap = NULL,
  useRaster = TRUE,
  rasterQuality = 5,
  split = NULL,
  fontSizeRows = 10,
  fontSizeCols = 10,
  fontSizeLabels = 8,
  colAnnoPerRow = 4,
  showRowDendrogram = FALSE,
  showColDendrogram = FALSE,
  customRowLabel = NULL,
  customRowLabelIDs = NULL,
  customColLabel = NULL,
  customColLabelIDs = NULL,
  customLabelWidth = 0.75,
  rasterDevice = "png",
  padding = 45,
  borderColor = NA,
  draw = TRUE,
  name = "Heatmap"
){

  #Packages
  requirePackage("ComplexHeatmap", source = "bioc")
  requirePackage("circlize", source = "cran")

  #Z-score
  if (scale) {
    message("Scaling Matrix..")
    mat <- rowZscores(mat, limit = FALSE)
    name <- paste0(name," Z-Scores")
  }

  #Get A Color map if null
  if (is.null(colorMap)) {
    colorMap <- colorMapAnno(colData)
  }

  #Prepare ColorMap format for Complex Heatmap
  if (!is.null(colData)){
    colData = data.frame(colData)
    colorMap <- colorMapForCH(colorMap, colData) #change
    showLegend <- checkShowLegend(colorMap[match(names(colorMap), colnames(colData))]) #change
  }else {
    colorMap <- NULL
    showLegend <- NULL
  }

  #Prepare Limits if needed
  breaks <- NULL
  if (!is.null(limits)) {
    mat[mat > max(limits)] <- max(limits)
    mat[mat < min(limits)] <- min(limits)
  }else{
    limits <- c(round(min(mat),2), round(max(mat),2))
  }

  #Scale Values 0 - 1
  mat <- (mat - min(limits)) / (max(limits) - min(limits))
  breaks <- seq(0, 1, length.out = length(color))
  color <- circlize::colorRamp2(breaks, color)

  if(exists('anno_mark', where='package:ComplexHeatmap', mode='function')){
    anno_check_version_rows <- ComplexHeatmap::anno_mark
    anno_check_version_cols <- ComplexHeatmap::anno_mark
  }else{
    anno_check_version_rows <- ComplexHeatmap::row_anno_link
    anno_check_version_cols <- ComplexHeatmap::column_anno_link
  }

  #Annotation Heatmap
  if(!is.null(colData) & !is.null(customColLabel)){
    message("Adding Annotations..")
    if(is.null(customColLabelIDs)){
      customColLabelIDs <- colnames(mat)[customColLabel]
    }
    ht1Anno <- HeatmapAnnotation(
      df = colData,
      col = colorMap,
      show_legend = showLegend,
      show_annotation_name = TRUE,
      gp = gpar(col = "NA"),
      annotation_legend_param =
        list(
          nrow = min(colAnnoPerRow, max(round(nrow(colData)/colAnnoPerRow), 1))
        ),
      foo = anno_check_version_cols(at = customColLabel, labels = customColLabelIDs, labels_gp = gpar(fontsize = fontSizeLabels))
    )

  }else if(!is.null(colData)){
    message("Adding Annotations..")
    ht1Anno <- HeatmapAnnotation(
      df = colData,
      col = colorMap,
      show_legend = showLegend,
      show_annotation_name = TRUE,
      gp = gpar(col = "NA"),
      annotation_legend_param =
        list(
          nrow = min(colAnnoPerRow, max(round(nrow(colData)/colAnnoPerRow), 1))
        )
    )
  }else if(is.null(colData) & !is.null(customColLabel)){
    if(is.null(customColLabelIDs)){
      customColLabelIDs <- colnames(mat)[customColLabel]
    }
    message("Adding Annotations..")
    #print(customColLabel)
    #print(customColLabelIDs)
    #ht1Anno <- columnAnnotation(foo = anno_check_version_cols(
    #   at = customColLabel, labels = customColLabelIDs),
    #   width = unit(customLabelWidth, "cm") + max_text_width(customColLabelIDs))
    #ht1Anno <- HeatmapAnnotation(foo = anno_mark(at = c(1:4, 20, 60, 1097:1100), labels = month.name[1:10]))
    ht1Anno <- HeatmapAnnotation(foo = anno_check_version_cols(at = customColLabel, labels = customColLabelIDs, labels_gp = gpar(fontsize = fontSizeLabels)))
  }else{
    ht1Anno <- NULL
  }

  message("Preparing Main Heatmap..")
  ht1 <- Heatmap(

    #Main Stuff
    matrix = as.matrix(mat),
    name = name,
    col = color,

    #Heatmap Legend
    heatmap_legend_param =
      list(
        at = c(0, 1),
        labels = c(round(min(limits),2), round(max(limits),2)),
        color_bar = "continuous",
        legend_direction = "horizontal",
        legend_width = unit(3, "cm")
      ),
    rect_gp = gpar(col = borderColor),

    #Column Options
    show_column_names = labelCols,
    cluster_columns = clusterCols,
    show_column_dend = showColDendrogram,
    clustering_method_columns = "ward.D2",
    column_names_gp = gpar(fontsize = fontSizeCols),
    column_names_max_height = unit(100, "mm"),

    #Row Options
    show_row_names = labelRows,
    row_names_gp = gpar(fontsize = fontSizeRows),
    cluster_rows = clusterRows,
    show_row_dend = showRowDendrogram,
    clustering_method_rows = "ward.D2",
    split = split,

    #Annotation
    top_annotation = ht1Anno,

    #Raster Info
    use_raster = useRaster,
    raster_device = rasterDevice,
    raster_quality = rasterQuality
  )

  if(!is.null(customRowLabel)){
    if(is.null(customRowLabelIDs)){
      customRowLabelIDs <- rownames(mat)[customRowLabel]
    }
    ht1 <- ht1 + rowAnnotation(link =
                                 anno_check_version_rows(at = customRowLabel, labels = customRowLabelIDs, labels_gp = gpar(fontsize = fontSizeLabels)),
                               width = unit(customLabelWidth, "cm") + max_text_width(customRowLabelIDs))
  }

  if(draw){
    draw(ht1,
         padding = unit(c(padding, padding, padding, padding), "mm"),
         heatmap_legend_side = "bot",
         annotation_legend_side = "bot")
  }else{
    ht1
  }

}


colorMapForCH <- function(colorMap = NULL, colData = NULL){
  colorMap <- colorMap[which(names(colorMap) %in% colnames(colData))]
  colorMapCH <- lapply(seq_along(colorMap), function(x){
    if(attr(colorMap[[x]],"discrete")){
      colorx <- colorMap[[x]]
    }else{
      vals <- colData[[names(colorMap)[x]]][!is.na(colData[[names(colorMap)[x]]])]
      s <-  seq(min(vals), max(vals), length.out = length(colorMap[[x]]))
      colorx <- circlize::colorRamp2(s, colorMap[[x]])
    }
    if(any(is.na(names(colorx)))){
      names(colorx)[is.na(names(colorx))] <- paste0("NA",seq_along(names(colorx)[is.na(names(colorx))]))
    }
    return(colorx)
  })
  names(colorMapCH) <- names(colorMap)
  return(colorMapCH)
}

checkShowLegend <- function(colorMap = NULL, max_discrete = 30){
  show <- lapply(seq_along(colorMap), function(x){
    if(attr(colorMap[[x]],"discrete") && length(unique(colorMap[[x]])) > max_discrete){
      sl <- FALSE
    }else{
      sl <- TRUE
    }
    return(sl)
  }) %>% unlist
  names(show) <- names(colorMap)
  return(show)
}

colorMapAnno <- function(colData = NULL, customAnno = NULL, discreteSet = "stallion", continuousSet = "solarExtra"){
  discreteCols <- sapply(colData,function(x) !is.numeric(x))
  if(!is.null(customAnno)){
    colorMap <- lapply(seq_along(discreteCols),function(x){
      if(discreteCols[x]){
        colors <- paletteDiscrete(values = colData[[names(discreteCols[x])]], set = discreteSet)
        names(colors) <- unique(colData[[names(discreteCols[x])]])
        attr(colors, "discrete") <- TRUE
      }else{
        colors <- paletteContinuous(set = continuousSet)
        attr(colors, "discrete") <- FALSE
      }
      if(length(which(customAnno[,1] %in% names(discreteCols[x]))) > 0){
        if(length(which(customAnno[,2] %in% names(colors))) > 0){
          customAnnox <- customAnno[which(customAnno[,2] %in% names(colors)),]
          colors[which(names(colors) %in% customAnnox[,2])] <- paste0(customAnnox[match(names(colors),customAnnox[,2]),3])
        }
      }
      return(colors)
    })
    names(colorMap) <- colnames(colData)
    return(colorMap)
  }else{
    colorMap <- lapply(seq_along(discreteCols), function(x){
      if(discreteCols[x]){
        colors <- paletteDiscrete(values = colData[[names(discreteCols[x])]], set = discreteSet)
        names(colors) <- unique(colData[[names(discreteCols[x])]])
        attr(colors, "discrete") <- TRUE
      }else{
        colors <- paletteContinuous(set = continuousSet)
        attr(colors, "discrete") <- FALSE
      }
      return(colors)
    })
    names(colorMap) <- colnames(colData)
    return(colorMap)
  }

}

binarySort <- function(m = NULL, scale = FALSE, cutOff = 1, lmat = NULL, clusterCols = TRUE){

  if(is.null(lmat)){
    #Compute Row-Zscores
    if(scale){
      lmat <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m), `/`)
    }else{
      lmat <- m
    }
    lmat <- lmat >= cutOff
  }

  #Transpose
  m <- t(m)
  lmat <- t(lmat)

  #Identify Column Ordering
  if(clusterCols){
    hc <- hclust(dist(m))
    colIdx <- hc$order
    m <- t(m[colIdx,])
    lmat <- t(lmat[colIdx,])
  }else{
    m <- t(m)
    lmat <- t(lmat)
    hc <- NULL
  }

  #Identify Row Ordering
  rowIdx <- do.call("order", c(as.data.frame(lmat)[seq_len(ncol(lmat))], list(decreasing = TRUE)))
  m <- t(m[rowIdx,])
  lmat <- t(lmat[rowIdx,])

  #Transpose
  m <- t(m)
  lmat <- t(lmat)

  return(list(mat = m, hclust = hc))

}

ArchRPalettes <- list(

  #DISCLOSURE: This is a collection of palettes that includes some original palettes and some palettes originally
  #implemented by others in other packages.
  #They are included here for convenience because they help improve plot aesthetics.

  #NOTE: all palettes included in the "Primarily Continuous Palettes" section should also work for discrete usage but not vice versa.
  #Each continuous palette has been ordered by color to generate a visually appealing discrete palette.

  #---------------------------------------------------------------
  # Primarily Discrete Palettes
  #---------------------------------------------------------------

  #20-colors
  stallion = c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
               "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D"),

  stallion2 = c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
                "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767"),

  calm = c("1"="#7DD06F", "2"="#844081", "3"="#688EC1", "4"="#C17E73", "5"="#484125", "6"="#6CD3A7", "7"="#597873","8"="#7B6FD0", "9"="#CF4A31", "10"="#D0CD47",
           "11"="#722A2D", "12"="#CBC594", "13"="#D19EC4", "14"="#5A7E36", "15"="#D4477D", "16"="#403552", "17"="#76D73C", "18"="#96CED5", "19"="#CE54D1", "20"="#C48736"),

  kelly = c("1"="#FFB300", "2"="#803E75", "3"="#FF6800", "4"="#A6BDD7", "5"="#C10020", "6"="#CEA262", "7"="#817066", "8"="#007D34", "9"="#F6768E", "10"="#00538A",
            "11"="#FF7A5C", "12"="#53377A", "13"="#FF8E00", "14"="#B32851", "15"="#F4C800", "16"="#7F180D", "17"="#93AA00", "18"="#593315", "19"="#F13A13", "20"="#232C16"),

  #16-colors
  bear = c("1"="#faa818", "2"="#41a30d","3"="#fbdf72", "4"="#367d7d",  "5"="#d33502", "6"="#6ebcbc", "7"="#37526d",
           "8"="#916848", "9"="#f5b390", "10"="#342739", "11"="#bed678","12"="#a6d9ee", "13"="#0d74b6",
           "14"="#60824f","15"="#725ca5", "16"="#e0598b"),

  #15-colors
  ironMan = c("9"='#371377',"3"='#7700FF',"2"='#9E0142',"10"='#FF0080', "14"='#DC494C',"12"="#F88D51","1"="#FAD510","8"="#FFFF5F","4"='#88CFA4',
              "13"='#238B45',"5"="#02401B", "7"="#0AD7D3","11"="#046C9A", "6"="#A2A475", "15"='grey35'),

  circus = c("1"="#D52126", "2"="#88CCEE", "3"="#FEE52C", "4"="#117733", "5"="#CC61B0", "6"="#99C945", "7"="#2F8AC4", "8"="#332288",
             "9"="#E68316", "10"="#661101", "11"="#F97B72", "12"="#DDCC77", "13"="#11A579", "14"="#89288F", "15"="#E73F74"),

  #12-colors
  paired = c("9"="#A6CDE2","1"="#1E78B4","3"="#74C476","12"="#34A047","11"="#F59899","2"="#E11E26",
             "10"="#FCBF6E","4"="#F47E1F","5"="#CAB2D6","8"="#6A3E98","6"="#FAF39B","7"="#B15928"),

  #11-colors
  grove = c("11"="#1a1334","9"="#01545a","1"="#017351","6"="#03c383","8"="#aad962","2"="#fbbf45","10"="#ef6a32","3"="#ed0345","7"="#a12a5e","5"="#710162","4"="#3B9AB2"),

  #7-colors
  summerNight = c("1"="#2a7185", "2"="#a64027", "3"="#fbdf72","4"="#60824f","5"="#9cdff0","6"="#022336","7"="#725ca5"),

  #5-colors
  zissou = c("1"="#3B9AB2", "4"="#78B7C5", "3"="#EBCC2A", "5"="#E1AF00", "2"="#F21A00"), #wesanderson
  darjeeling = c("1"="#FF0000", "2"="#00A08A", "3"="#F2AD00", "4"="#F98400", "5"="#5BBCD6"), #wesanderson
  rushmore = c("1"="#E1BD6D", "5"="#EABE94", "2"="#0B775E", "4"="#35274A" , "3"="#F2300F"), #wesanderson
  captain = c("1"="grey","2"="#A1CDE1","3"="#12477C","4"="#EC9274","5"="#67001E"),

  #---------------------------------------------------------------
  # Primarily Continuous Palettes
  #---------------------------------------------------------------

  #10-colors
  horizon = c("1"='#000075',"4"='#2E00FF', "6"='#9408F7', "10"='#C729D6', "8"='#FA4AB5', "3"='#FF6A95', "7"='#FF8B74', "5"='#FFAC53', "9"='#FFCD32', "2"='#FFFF60'),

  #9-colors
  horizonExtra =c("1"="#000436","4"="#021EA9","6"="#1632FB","8"="#6E34FC","3"="#C732D5","9"="#FD619D","7"="#FF9965","5"="#FFD32B","2"="#FFFC5A"),
  blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D"),
  sambaNight = c("6"='#1873CC',"2"='#1798E5',"8"='#00BFFF',"5"='#4AC596',"1"='#00CC00',"4"='#A2E700',"9"='#FFFF00',"7"='#FFD200',"3"='#FFA500'), #buencolors
  solarExtra = c("5"='#3361A5', "7"='#248AF3', "1"='#14B3FF', "8"='#88CEEF', "9"='#C1D5DC', "4"='#EAD397', "3"='#FDB31A',"2"= '#E42A2A', "6"='#A31D1D'),  #buencolors
  whitePurple = c("9"='#f7fcfd',"6"='#e0ecf4',"8"='#bfd3e6',"5"='#9ebcda',"2"='#8c96c6',"4"='#8c6bb1',"7"='#88419d',"3"='#810f7c',"1"='#4d004b'),
  whiteBlue = c("9"='#fff7fb',"6"='#ece7f2',"8"='#d0d1e6',"5"='#a6bddb',"2"='#74a9cf',"4"='#3690c0',"7"='#0570b0',"3"='#045a8d',"1"='#023858'),
  whiteRed = c("1"="white", "2"="red"),
  comet = c("1"="#E6E7E8","2"="#3A97FF","3"="#8816A7","4"="black"),

  #7-colors
  greenBlue = c("4"='#e0f3db',"7"='#ccebc5',"2"='#a8ddb5',"5"='#4eb3d3',"3"='#2b8cbe',"6"='#0868ac',"1"='#084081'),

  #6-colors
  beach = c("4"="#87D2DB","1"="#5BB1CB","6"="#4F66AF","3"="#F15F30","5"="#F7962E","2"="#FCEE2B"),

  #5-colors
  coolwarm = c("1"="#4858A7", "4"="#788FC8", "5"="#D6DAE1", "3"="#F49B7C", "2"="#B51F29"),
  fireworks = c("5"="white","2"="#2488F0","4"="#7F3F98","3"="#E22929","1"="#FCB31A"),
  greyMagma = c("2"="grey", "4"="#FB8861FF", "5"="#B63679FF", "3"="#51127CFF", "1"="#000004FF"),
  fireworks2 = c("5"="black", "2"="#2488F0","4"="#7F3F98","3"="#E22929","1"="#FCB31A"),
  purpleOrange = c("5"="#581845", "2"="#900C3F", "4"="#C70039", "3"="#FF5744", "1"="#FFC30F")
)

paletteDiscrete <- function(
  values = NULL,
  set = "stallion",
  reverse = FALSE
){

  values <- gtools::mixedsort(values)
  n <- length(unique(values))
  pal <- ArchRPalettes[[set]]
  palOrdered <- pal[gtools::mixedsort(names(pal))] #mixed sort gets 1,2,3,4..10,11,12

  if(n > length(palOrdered)){
    message("Length of unique values greater than palette, interpolating..")
    palOut <- colorRampPalette(pal)(n)
  }else{
    palOut <- palOrdered[seq_len(n)]
  }

  if(reverse){
    palOut <- rev(palOut)
  }

  names(palOut) <- unique(values)

  return(palOut)

}

paletteContinuous <- function(
  set = "solarExtra",
  n = 256,
  reverse = FALSE
){

  pal <- ArchRPalettes[[set]]
  palOut <- colorRampPalette(pal)(n)

  if(reverse){
    palOut <- rev(palOut)
  }

  return(palOut)

}

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


#' Convert scRNA SummarizedExperiment Object into RangedSummarizedExperiment Object
#' @param seRNA SummarizedExperiment Object
#' @param gtfFile gene gtf file
#' @param tssWindow tssWindow's number,defaut 2500
#' @return RangedSummarizedExperiment object
#' @export
seRNA2Rse<-function(seRNA=NULL,
                    gtf="gene.gtf",
                    tssWindow =2500){
  #require("GenomicRanges")
  #require("SummarizedExperiment")

  #genes <- getGeneGTF(gtfFile) %>% GenomicRanges::resize(1,"start") %>% GenomicRanges::resize(tssWindow * 2 + 1, "center")
  if(class(gtf)=="GRanges"){
    genes <- gtf %>% GenomicRanges::resize(1,"start") %>% GenomicRanges::resize(tssWindow * 2 + 1, "center")
  }
  else if(class(gtf)=="character"){
    genes <- getGeneGTF(gtf) %>% GenomicRanges::resize(1,"start") %>% GenomicRanges::resize(tssWindow * 2 + 1, "center")
  }else{
    stop("Invalid gtf!!")
  }
  ##seRNA TO RseRNA
  names(genes) <- genes$gene_name
  features<-rownames(seRNA)
  features<-features[features%in%mcols(genes)$gene_name]
  seRNA <-seRNA[features,]
  rowranges <- genes[rownames(seRNA),]
  rowRanges(seRNA)<-rowranges
  return(seRNA)
}

#' rename RangedSummarizedExperiment seATAC
#' @param seATAC RangedSummarizedExperiment Object(eg,from getMatrixFromProject )
gr2Feature<-function(seATAC=NULL){
  chr<-SummarizedExperiment::seqnames(seATAC)
  st<-SummarizedExperiment::start(seATAC)
  ed<-SummarizedExperiment::end(seATAC)
  peak<-paste(chr,st,ed,sep="_")
  rownames(seATAC)<-peak
  seATAC
}

#' Convert Seurat object into SummarizedExperiment Object
#' @param object seurat object
#' @return SummarizedExperiment Object
#' @export
seuratTose<-function(object=NULL){
  #require("Seurat")
  metadata<-object@meta.data
  logcount<-GetAssayData(object,"data")
  count<-GetAssayData(object,"counts")

  print("Converting...")
  se<-SummarizedExperiment::SummarizedExperiment(assays=list(counts=count,logcounts=logcount),
                           colData=metadata)
  return(se)

}

.safelapply<-function (..., threads = 1, preschedule = FALSE)
{
  require(parallel)
  if (tolower(.Platform$OS.type) == "windows") {
    threads <- 1
  }
  if (threads > 1) {
    o <- mclapply(..., mc.cores = threads, mc.preschedule = preschedule)
    errorMsg <- list()
    for (i in seq_along(o)) {
      if (inherits(o[[i]], "try-error")) {
        capOut <- utils::capture.output(o[[i]])
        capOut <- capOut[!grepl("attr\\(\\,|try-error",
                                capOut)]
        capOut <- head(capOut, 10)
        capOut <- unlist(lapply(capOut, function(x) substr(x,
                                                           1, 250)))
        capOut <- paste0("\t", capOut)
        errorMsg[[length(errorMsg) + 1]] <- paste0(c(paste0("Error Found Iteration ",
                                                            i, " : "), capOut), "\n")
      }
    }
    if (length(errorMsg) != 0) {
      errorMsg <- unlist(errorMsg)
      errorMsg <- head(errorMsg, 50)
      errorMsg[1] <- paste0("\n", errorMsg[1])
      stop(errorMsg)
    }
  }
  else {
    o <- lapply(...)
  }
  o
}
