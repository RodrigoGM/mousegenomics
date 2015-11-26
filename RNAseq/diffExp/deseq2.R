#############################################################
##  This code is inteded to analyze differential expression
## from RNAseq mouse tissue data. Data was mapped to mm10 with
## tophat2, then sorted with samtools sort -n, inexed with
## samtools index.  All files were analyzed with picard-tools
## CollectRnaseqMetrics.jar, then ran by cufflinks.  Gene
## counts were performed with HTSeq using -m insertion-strict
## -s 'reverse'
#############################################################

library(GenomicFeatures)
library(DESeq2)
library(RColorBrewer)
library(biomaRt)
cols <- brewer.pal(5, "Set1")[c(4:5,3,1,2)]

## Grcm38 Esembl 75 Genes 
# mme <- makeTranscriptDbFromGFF("~/mmu/Grcm38/Annotation/Genes/Grcm38.genes.gtf", format = "gtf", species="mus_musculus")

## creating count matrix for each tissue.  Uses the file.path to select
##  for each tissue/ the library/ and specifically the 2pass/ star.gene.htseq.txt
##  count file.  The count.files is sorted to have CTRL libraries at the
##  end of the batch
tissues <- c("pituitary", "liver", "heart", "embryo", "placenta")
count.files <- unlist(lapply(tissues, function(TSS) file.path(TSS, list.files(path = TSS, pattern = ".*L.*"), "2pass", "star.gene.htseq.txt")))[c(1:16,19:20,17:18,21:26,29:30,27:28,31:42,47:50,43:46,51:62,67:70,63:66)]

## confirming that all count.files exist
if(any(!sapply(count.files, file.exists))) warning("You have a missing count file in your expected data set")

## reading all count.files into a data.frame. We assume that all count.files
##  have the genes rows in the same order.
htcount <- data.frame(do.call(cbind, lapply(count.files, read.delim, header = FALSE, row.names = 1)))

## Adding appropriate names i.e. the library name as column names.
colnames(htcount) <- gsub("_L0.*/2pass/star.gene.htseq.txt", "", count.files)
colnames(htcount) <- gsub(".*/", "", colnames(htcount))
colnames(htcount) <- gsub("ctrl", "CTRL", colnames(htcount))

## removing the stats provided by HTSeq-count
htcount.stats <- htcount[grep("^__", rownames(htcount)), ]
htcount <- htcount[-grep("^__", rownames(htcount)), ]
htcount.stats <- rbind(htcount.stats, colSums(htcount))
rownames(htcount.stats)[6] <- "__total_aligned"
colSums(htcount) / (colSums(htcount.stats) + colSums(htcount))

## cssDesign : a data.frame containing the experiment metadata e.g. library names, strains, library-type, tissues, etc...
##      libraryName : name of each library as <StrainABV>_<Tissue>-<Pool>
##      strain : Long Strain Name as factor
##      libType : library-type as described in tophat
##      tissue :  tissue from library
strain <- factor(gsub("ctrl", "CTRL", gsub("_PL", "", gsub(".*/(.*[1Ll])[-_].*", "\\1", count.files))),
                 labels = c("B6.A-15BC1", "B6.A-17BC1", "B6.A-19BC1", "B6.CTRL", "B6.A-XF1"),
                 ordered = FALSE)

cssDesign <- data.frame(libraryName = colnames(htcount), strain = strain, libType = rep("fr-firststrand", 70), tissue = rep(tissues, c(10, 10, 10, 20, 20)), File = count.files, method=rep('2pass', 70))
rownames(cssDesign)  <- gsub(".*/(.*)/2pass/star.gene.htseq.txt", "\\1", count.files)

## Creating the class: DESeqDataSet objects for the complete raw count data from HTSeq
## ddsFCT : ddsFullCountTable : complete raw count data from HTSeq,
##  with metadata, and experimental design.  Formula/Design is `~ strain`, as it will
##  be separated by tissue at a later step
## dds : collapsed data set.  This data adds both technical replicates
##  from embryo and placenta thus from 70 to 50 libraries as desired
ddsFCT <- DESeqDataSetFromMatrix(countData = htcount,  colData = cssDesign, design = ~ strain)
ddsCR <- collapseReplicates(ddsFCT, groupby = ddsFCT$libraryName)

## RE-orders the levels of the strain factos, such that B6.CTRL is the first one.
## recomended in the beginner's guide to DESeq2
ddsCR$strain <- relevel(ddsCR$strain, "B6.CTRL")
as.data.frame(colData(ddsCR))

## Separating the collapsed data set into individual tissues,
## Removing non-used levels from the tissue column
## Confirming that B6.CTRL is the first level as the first level is used
##  for comparison 
dds <- lapply(tissues, function(TSS) ddsCR[ , ddsCR$tissue == TSS ])
sapply(1:5, function(i) dds[[i]]$tissue <<- droplevels(dds[[i]]$tissue))
sapply(1:5, function(i) dds[[i]]$strain <<- relevel(dds[[i]]$strain, "B6.CTRL"))
names(dds) <- tissues

## Running the DESeq pipeline in one step : this function :
##  1. estimates size factors w/ `estimateSizeFactors`
##  2. estimates dispersions w/ `estimateDispertions`
##  3. runs GLM w/ `nbinomWaldTest`
invisible(sapply(tissues, function(i) dds[[i]] <<- DESeq(dds[[i]], fitType = "parametric")))
htcount.norm <- lapply(dds, counts, normalized = TRUE)

## Getting results by tissue in a per strain vs CTRL basis
pituitary.de <- lapply(as.character(unique(strain))[1:4], function(STR) results(dds[["pituitary"]], contrast = c("strain", STR, "B6.CTRL")))
liver.de <- lapply(as.character(unique(strain))[1:4], function(STR) results(dds[["liver"]], contrast = c("strain", STR, "B6.CTRL")))
heart.de <- lapply(as.character(unique(strain))[1:4], function(STR) results(dds[["heart"]], contrast = c("strain", STR, "B6.CTRL")))
embryo.de <- lapply(as.character(unique(strain))[1:4], function(STR) results(dds[["embryo"]], contrast = c("strain", STR, "B6.CTRL")))
placenta.de <- lapply(as.character(unique(strain))[1:4], function(STR) results(dds[["placenta"]], contrast = c("strain", STR, "B6.CTRL")))

## adding MGI gene symbols to each results objects
##  uses biomaRt package
##
gtf.genes <- sapply(strsplit(rownames(pituitary.de[[1]]), split = "\\+"), "[", 1)
mmu <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
thesaurus <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
              filters = "ensembl_gene_id",
              values = gtf.genes,
              mart = mmu)
idx <- match(gtf.genes, thesaurus$ensembl_gene_id)

invisible(
sapply(1:4, function(s) {
  pituitary.de[[s]]$mgi.symbol <<- thesaurus$mgi_symbol[idx]
  liver.de[[s]]$mgi.symbol <<- thesaurus$mgi_symbol[idx]
  heart.de[[s]]$mgi.symbol <<- thesaurus$mgi_symbol[idx]
  embryo.de[[s]]$mgi.symbol <<- thesaurus$mgi_symbol[idx]
  placenta.de[[s]]$mgi.symbol <<- thesaurus$mgi_symbol[idx]
})
)
save.image("css.deseq2.rda")


for(i in 1:4) {
  pituitary.de[[i]]$comparison <- gsub(".*strain", "", attr(pituitary.de[[i]], "elementMetadata")[2,2])
  pituitary.de[[i]]$tissue <- "pituitary"
  pituitary.de[[i]]$id <- rownames(pituitary.de[[i]])
}
pituitary.df <- do.call(rbind, pituitary.de)
write.csv(pituitary.df, file = "PituitaryDE.csv", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
for(i in 1:4) {
  liver.de[[i]]$comparison <- gsub(".*strain", "", attr(liver.de[[i]], "elementMetadata")[2,2])
  liver.de[[i]]$tissue <- "liver"
  liver.de[[i]]$id <- rownames(liver.de[[i]])
  
}
liver.df <- do.call(rbind, liver.de)
write.csv(liver.df, file = "LiverDE.csv", quote = FALSE, row.names = TRUE, col.names = TRUE)

for(i in 1:4) {
  heart.de[[i]]$comparison <- gsub(".*strain", "", attr(heart.de[[i]], "elementMetadata")[2,2])
  heart.de[[i]]$tissue <- "heart"
  heart.de[[i]]$id <- rownames(heart.de[[i]])
}
heart.df <- do.call(rbind, heart.de)
write.csv(heart.df, file = "HeartDE.csv", quote = FALSE, row.names = TRUE, col.names = TRUE)

for(i in 1:4) {
  embryo.de[[i]]$comparison <- gsub(".*strain", "", attr(embryo.de[[i]], "elementMetadata")[2,2])
  embryo.de[[i]]$tissue <- "embryo"
  embryo.de[[i]]$id <- rownames(embryo.de[[i]])
}
embryo.df <- do.call(rbind, embryo.de)
write.csv(embryo.df, file = "EmbryoDE.csv", quote = FALSE, row.names = TRUE, col.names = TRUE)

for(i in 1:4) {
  placenta.de[[i]]$comparison <- gsub(".*strain", "", attr(placenta.de[[i]], "elementMetadata")[2,2])
  placenta.de[[i]]$tissue <- "placenta"
  placenta.de[[i]]$id <- rownames(placenta.de[[i]])
}
placenta.df <- do.call(rbind, placenta.de)
write.csv(placenta.df, file = "PlacentaDE.csv", quote = FALSE, row.names = TRUE, col.names = TRUE)

## merge all
cssDE <- rbind(pituitary.df, liver.df, heart.df, embryo.df, placenta.df)
write.csv(cssDE, file = "cssDE.csv", quote = FALSE, row.names = FALSE)

sigDE <- cssDE[cssDE$id <= 0.00001 , ]
## read in significant genes
#sigDE <- read.csv("Significant_cssDE.csv")

save.image("css.deseq2.rda")

## 
all.genes <- list(pituitary.df, liver.df, heart.df, embryo.df, placenta.df)
names(all.genes) <- tissues

##
sig.genes <- lapply(tissues, function(i) as.character(sigDE$id[sigDE$tissue == i]))
names(sig.genes) <- tissues

## 
sig.genes <- lapply(tissues, function(SG) all.genes[[SG]][as.character(all.genes[[SG]]$id) %in% sig.genes[[SG]], ])
sig.genes <- do.call(rbind, data.frame(sig.genes))

write.csv(sig.genes, "ForExp_Matrix.csv")

