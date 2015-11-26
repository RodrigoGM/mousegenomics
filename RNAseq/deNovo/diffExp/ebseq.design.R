## using ebseq for matrix.* from RSEM
library(EBSeq)

## general experimental design
## We create a design table based on the samples, the Rg*_series.txt files are a tab delimited file that contains the
## path/to/sample, and treatment or additional information needed
qDesign <- do.call(rbind, lapply(c("../Rg1_series.txt", "../Rg2_series.txt", "../Rg3_series.txt"), read.table, row.names = 1))
colnames(qDesign) <- c("sex", "trt", "group")
qDesign$region <- factor(c(rep("Rg1", 11), rep("Rg2", 12), rep("Rg3", 12)))


## gene count matrix
## Create a read counts matrix of the "trinity genes" from all the the samples.  RSEM creates the gene count matrix, here we just read it in.
GeneMat <- data.matrix(read.delim("genes/genes.counts.matrix", row.names = 1)[,-6])
( colnames(GeneMat) <- gsub("NGS15.D[0-9]{3}_([MR].*)", "\\1", colnames(GeneMat)) )  ## changing sample names to match those in qDesign
nrow(GeneMat)  ## total number of rows pre selection

GeneMat <- GeneMat[rowSums(GeneMat) >= 100, ]  ## keep those whose total read count is greater than 100
GeneMat <- GeneMat + 1 ## adding 1 to avoid zeros

nrow(GeneMat)  ## number of rows after selection

( Sizes <- MedianNorm(GeneMat) )

## Isoform level
IsoMat <- data.matrix(read.delim("isoforms/isoforms.counts.matrix", row.names = 1)[,-6])
( colnames(IsoMat) <- gsub("NGS15.D[0-9]{3}_([MR].*)", "\\1", colnames(IsoMat)) ) ## changing sample names to mactch those in qDesign

IsoMat <- IsoMat[rowSums(IsoMat) >= 100, ]
IsoMat <- IsoMat + 1

IsoNames <- gsub("(TR.*)\\|(c.*)", "\\2", rownames(IsoMat))
IsoGeneNames <- gsub("(TR.*)\\|(c.*)", "\\1", rownames(IsoMat))
IsosGeneNames <- IsoGeneNames

IsoSizes <- MedianNorm(IsoMat)
NgList <- GetNg(IsoNames, IsoGeneNames)
IsoNgTrun <- NgList$IsoformNgTrun

save.image("qEB.data.rda")


## region specific subsets
Rg1 <- GeneMat[ , grep("Rg1", colnames(GeneMat))]
Rg2 <- GeneMat[ , grep("Rg2", colnames(GeneMat))]
Rg3 <- GeneMat[ , grep("Rg3", colnames(GeneMat))]

Rg1Iso <- IsoMat[ , grep("Rg1", colnames(IsoMat))]
Rg2Iso <- IsoMat[ , grep("Rg2", colnames(IsoMat))]
Rg3Iso <- IsoMat[ , grep("Rg3", colnames(IsoMat))]

## qDesign subset
qRG1 <- qDesign[grep("Rg1", rownames(qDesign)) , ]
qRG2 <- qDesign[grep("Rg2", rownames(qDesign)) , ]
qRG3 <- qDesign[grep("Rg3", rownames(qDesign)) , ]

## estimating sizes
Rg1.Sizes <- MedianNorm(Rg1)
Rg2.Sizes <- MedianNorm(Rg2)
Rg3.Sizes <- MedianNorm(Rg3)

## estimating sizes
Rg1Iso.Sizes <- MedianNorm(Rg1Iso)
Rg2Iso.Sizes <- MedianNorm(Rg2Iso)
Rg3Iso.Sizes <- MedianNorm(Rg3Iso)

save.image("qEB.data.rda")

sessionInfo()
