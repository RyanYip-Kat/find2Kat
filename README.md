# find2Kat
some functions for bed or Bam files ,or Granges object
### Description
Find2Kat this package, can implement calleak and annote bam, bed file, but also can implement the ArchR object DORC score calculation (generally only on the share-seq data), as well as can achieve the correlation between cell type and disease and so on
...

### installtation
```r
install.packages("find2Kat_0.2.0.tar.gz")
```
you also need to install depend packages via *require-packages.txt*

### callpeak for BAM or BED file and annotation
you can use function **annoteFile** reliaze
```r
annoteFile(File="yourPathOfFiles",Names="relativeNames",blacklist=find2Kat::blacklist[["hg38"]],genome="hg38",...)
```
and the result like follow 

![annotation](Figures/annoteFile.png)
