###############################################################
##  This code is inteded to visualize differential expression
## from RNAseq mouse tissue data. Data was mapped to mm10 with
## tophat2, then sorted with samtools sort -n, inexed with
## samtools index.  All files were analyzed with picard-tools
## CollectRnaseqMetrics.jar, then ran by cufflinks.  Gene
## counts were performed with HTSeq using -m insertion-strict
## -s 'reverse'
###############################################################

library(DESeq2)
library(RColorBrewer)
#library(venneuler)
library(biomaRt)
library(gplots)
library(fields)

## loading data or source
load("css.deseq2.rda")
##source("deseq2.R")

## list of genes from first analysis
gg <- c("character", rep("numeric", 6), rep("character", 3))
mf <- read.delim("MichelsFavorites.txt", skip = 2, header = TRUE, colClasses = gg)
sig.genes <- read.delim("geneList.txt", skip = 2, header = TRUE, colClasses = gg)

## add genes as attributes to each tissue/strain combination
comparisons <- unique(pituitary.df$comparison)
sapply(1:4, function(i) {
  attr(pituitary.de[[i]], "sig.genes" ) <<- sig.genes$id[sig.genes$tissue == "pituitary" & sig.genes$comparison == comparisons[i]]
  attr(pituitary.de[[i]], "mf" ) <<- mf$id[mf$tissue == "pituitary" & mf$comparison == comparisons[i]]
})

sapply(1:4, function(i) {
  attr(liver.de[[i]], "sig.genes" ) <<- sig.genes$id[sig.genes$tissue == "liver" & sig.genes$comparison == comparisons[i]]
  attr(liver.de[[i]], "mf" ) <<- mf$id[mf$tissue == "liver" & mf$comparison == comparisons[i]]
})

sapply(1:4, function(i) {
  attr(heart.de[[i]], "sig.genes" ) <<- sig.genes$id[sig.genes$tissue == "heart" & sig.genes$comparison == comparisons[i]]
  attr(heart.de[[i]], "mf" ) <<- mf$id[mf$tissue == "heart" & mf$comparison == comparisons[i]]
})

sapply(1:4, function(i) {
  attr(embryo.de[[i]], "sig.genes" ) <<- sig.genes$id[sig.genes$tissue == "embryo" & sig.genes$comparison == comparisons[i]]
  attr(embryo.de[[i]], "mf" ) <<- mf$id[mf$tissue == "embryo" & mf$comparison == comparisons[i]]
})

sapply(1:4, function(i) {
  attr(placenta.de[[i]], "sig.genes" ) <<- sig.genes$id[sig.genes$tissue == "placenta" & sig.genes$comparison == comparisons[i]]
  attr(placenta.de[[i]], "mf" ) <<- mf$id[mf$tissue == "placenta" & mf$comparison == comparisons[i]]
})


## PercentMapped:  Counts are binned by htseq-count into six
##  different types of statistics : no_feature, ambiguous,
##  too_low_aQual, not_aligned, alignent_not_unique, and
##  total_aligned.  This plot summarizes the counted reads,
##  and the statistics
pdf("HTSeqPercentCounted.pdf", width = 200/25.4, height = 200/25.4, useDingbats = FALSE)
#png("HTSeqPercentCounted.png", width = 200, height = 200, units = "mm", res = 400)
par(mar = c(10, 5, 2, 2)+.1)
plot(colSums(htcount) / (colSums(htcount) + colSums(htcount.stats)),
  pch = rep(15:19, each = 2),  cex = 1, col = rep(cols, each = 2),
  xlab = "", ylab = "HTSeq Percent Counted",
  xaxt = "n", bty = "l")
abline(v=seq(from = 0.5, to = 70.5, by = 2), col = "grey")
axis(1, at = 1:70, labels = colnames(htcount), las = 2)
mtext("Pool", side = 1, line = 6.5)
dev.off()


## Quality control for each agregated replications
png("HTCount_Replications%03d.png", width = 320, height = 320/5, units = "mm", res = 500)
layout(mat = matrix(1:5, nrow = 1, byrow = TRUE))
par(mar = c(5, 4, 1, 1))
sapply(c("pituitary", "liver", "heart", "embryo", "placenta"), function(T) {
  sapply(seq(1, 9, 2), function(i) {
    plot(log2(htcount.norm[[T]][, i+1]) ~ log2(htcount.norm[[T]][, i]), col = "gray", pch = ".",
    xlab = paste("Log2 ", as.character(cssDesign[i, "strain"]), "A", sep = " "),
    ylab = paste("Log2 ", as.character(cssDesign[i, "strain"]), "B", sep = " "),
    xlim = c(0, 20), ylim = c(0, 20), xaxs = "i", yaxs = "i", xlog = TRUE, ylog = TRUE,
    )
    abline(a = 0, b = 1, col = 'red')
  })
})
dev.off()


## Plot dispersions for everything under the sun !
png("Dispersions_by_tissue.png", width = 290, height = 290/4.5, units = "mm", res = 450)
layout(mat = matrix(c(1:5), nrow = 1, byrow = TRUE))
sapply(1:5, function(t) {
  plotDispEsts(dds[[t]])
  title(main = paste("Dispersion for", tissues[t]))
})
dev.off()


## Generating MA-plot for each tissue per strain
png("MAplot_by_strain%03d.png", width = 320, height = 320/5, units = "mm", res = 450)
sapply(list(pituitary.de, liver.de, heart.de, embryo.de, placenta.de), function(TDE) {
layout(mat = matrix(c(1:4), nrow = 1, byrow = TRUE))
  par(mar = c(5, 4, 1, 1)+.1)
  sapply(TDE, function(T) {
    plotMA(T, ylim = c(-8, 8), xlab = "Mean Expression", ylab = expression(Log[2]~Fold~Change), xaxt = "n", yaxt = "n")
    points(x = T[attr(T, "sig.genes"), "baseMean"], T[attr(T, "sig.genes"), "log2FoldChange"], col = "blue", pch = 19, cex = .69)
    points(x = T[attr(T, "mf"), "baseMean"], T[attr(T, "mf"), "log2FoldChange"], col = "red", pch = 1, cex = 1)
    axis(side = 1, at = c(0.1, 10, 1000, 100000, 10000000),
      labels = c(expression(1%*%10^-1), expression(1%*%10^1), expression(1%*%10^3),
        expression(1%*%10^5), expression(1%*%10^7)))
    axis(side = 2, at = seq(from = -8, to = 8, by = 2), labels = TRUE)
  })
legend("topright", legend = c(expression(p<=1%*%10^-5), expression(FDR<=0.05)), col = c("blue", "red"), pch = c(19, 19), cex = c(1, 1))
})
dev.off()


## QQ plot generation :
##  Create a list of all *.de objects in the order we want plotted.
##    loop through each strain contrast on the *.de object and obtain
##    its observed p-value.  Create a vector with the same length of
##    a uniform distribution i.e. a sequence of numbers from 1:length(obs)
##    at uniform intervals.  Add detail i.g. uniform line, text, nice
##    legend font, etc...
png("QQplot_by_tissue%03d.png", width = 320, height = 320/5, units = "mm", res = 450)
sapply(list(pituitary.de, liver.de, heart.de, embryo.de, placenta.de), function(TDE) {
layout(mat = matrix(c(1:4), nrow = 1, byrow = TRUE))
  sapply(TDE, function(STR) {
    obs = STR$pvalue
    names(obs) = rownames(STR)
    obs = sort(obs)
    lexpp = log10(1:length(obs) / (length(obs)+1)) * -1
    par(mar = c(5, 4, 1, 1)+.1)
    plot(-log10(obs) ~ lexpp, xlab = expression(Expected~italic(-Log[10]~P-value)),
    	 ylab = expression(italic(-Log[10]~P-value)), pch = 20
#	 main = paste("QQ Plot :", gsub(".*strain ", "", attr(STR, "elementMetadata")[2,2] ))
    )
    points(x = STR[attr(STR, "sig.genes"), "baseMean"], STR[attr(STR, "sig.genes"), "log2FoldChange"], col = "blue", pch = 19, cex = .69)
    points(x = STR[attr(STR, "mf"), "baseMean"], STR[attr(STR, "mf"), "log2FoldChange"], col = "red", pch = 1, cex = 1)
    abline(a = 0, b = 1, col = "blue", lwd = 2)
    abline(h = -log10(0.00001), lty = 2)
    text(x = 0.1, y = -log10(0.00001)+.5, labels = expression(italic(p-value)==italic(-log[10])~5), pos = 4)
#    text(x = max(lexpp), y = max(-log10(obs)), labels = rownames(STR[min(STR$pvalue)]), cex =.3, pos = 2)
  })
})
dev.off()

png("QQfdr_by_tissue%03d.png", width = 320, height = 320/5, units = "mm", res = 450)
sapply(list(pituitary.de, liver.de, heart.de, embryo.de, placenta.de), function(TDE) {
layout(mat = matrix(c(1:4), nrow = 1, byrow = TRUE))
  sapply(TDE, function(STR) {
    obs = STR$padj
    names(obs) = rownames(STR)
    obs = sort(obs)
    lexpp = log10(1:length(obs) / (length(obs)+1)) * -1
    par(mar = c(5, 4, 1, 1)+.1)
    plot(-log10(obs) ~ lexpp, xlab = expression(Expected~italic(-Log[10]~q-value)),
    	 ylab = expression(italic(-Log[10]~q-value)), pch = 20
#	 main = paste("QQ Plot :", gsub(".*strain ", "", attr(STR, "elementMetadata")[2,2] ))
    )
    points(x = STR[attr(STR, "sig.genes"), "baseMean"], STR[attr(STR, "sig.genes"), "log2FoldChange"], col = "blue", pch = 19, cex = .69)
    points(x = STR[attr(STR, "mf"), "baseMean"], STR[attr(STR, "mf"), "log2FoldChange"], col = "red", pch = 1, cex = 1)
    abline(a = 0, b = 1, col = "red", lwd = 2)
    abline(h = -log10(0.05), lty = 2)
    text(x = 0.1, y = -log10(0.05)+.1, labels = expression(italic(p-value)==italic(FDR%<=%5)), pos = 4)
#    text(x = max(lexpp), y = max(-log10(obs)), labels = rownames(STR[min(STR$pvalue)]), cex =.3, pos = 2)
  })
})
dev.off()



## VennEuler diagrams for comparisons & tissues
##  First create the venneuler object for both comparisons and
##  tissues.  Also create a summary object with the
##  counts for each tissue by strain comparison
by.str <- venneuler(as.matrix(sigDE[, c("gene.id", "comparison")]))
by.tissue <- venneuler(as.matrix(sigDE[, c("gene.id", "tissue")]))
numbers <- table(sigDE$comparison, sigDE$tissue)

pdf("HTCount_StrainVennEuler.pdf")
plot(by.str)
text(x = by.str$centers[,1], y = by.str$centers[,2], labels = as.numeric(rowSums(numbers)[by.str$labels]), pos = 1)
dev.off()

pdf("HTCount_TissueVennEuler.pdf")
plot(by.tissue)
text(x = by.tissue$centers[,1], y = by.tissue$centers[,2], labels = as.numeric(colSums(numbers)[by.tissue$labels]), pos = 1)
dev.off()

## load expression matrix for heatmap
exp.mat <- read.csv("ExpressionMatrix_SignificantGenes.csv")

pdf("HTCount_heatmap.pdf", width = 297/25.4, height = 148/25.4)
heatmap.2(t(as.matrix(exp.mat[,2:5])), dendrogram = "none", Rowv = NA,
          Colv = NA, col = tim.colors(30), main = "",
          margins = c(8, 2), labRow = c("B6.A-15 N2 B/B", "B6.A-17 N2 B/B", "B6.A-19 N2 B/B", "B6.A-X F1 B/B"),
          labCol = as.character(exp.mat$id), xlab = "", cexRow = 2.8 , cexCol = .69,
          trace = "none", lmat = matrix(c(0, 2, 3, 1, 0, 4), nrow = 2, ncol = 3),
          lhei = c(.4, 3), lwid = c(.5, 4.5, 2), key = TRUE)
image.plot(zlim = c(-4, 4), legend.only = TRUE, col = tim.colors(30), legend.mar = 5, legend.width = 1.5)
dev.off()

means.htcount.norm <- lapply(file.path("RNAseq_ReadCounts", list.files("RNAseq_ReadCounts/", pattern = "means\\+range")), read.csv, row.names = 1)
names(means.htcount.norm) <- c("embryo", "heart", "liver", "pituitary", "placenta")

cols <- brewer.pal(5, "Set1")[c(4,5,3,1,2)]  ## colors for each strain
#colt <-                                    ## colors for each tissue
apply(means.htcount.norm$heart[sig.genes$id, ], 1, function(x) {
  barplot2(x[11:15], plot.ci = TRUE, ci.l = x[16:20], ci.u = x[21:25],
    main = sig.genes[rownames(x), "mgi.symbol"], col = c("#984EA3", "#FF7F00", "#4DAF4A", "#E41A1C", "#377EB8"))
})


ymax = sapply(rownames(means.htcount.norm[[1]]), function(g) max(c(as.numeric(means.htcount.norm[[1]][g, 21:25]), as.numeric(means.htcount.norm[[2]][g, 21:25]), as.numeric(means.htcount.norm[[3]][g, 21:25]), as.numeric(means.htcount.norm[[4]][g, 21:25]), as.numeric(means.htcount.norm[[5]][g, 21:25]))))


means.as.list.mat <- lapply(unique(sig.genes$id), function(g) as.matrix(do.call(rbind, lapply(means.htcount.norm, function(i) i[g, 11:15])))[tissues,])
lci.as.list.mat <- lapply(unique(sig.genes$id), function(g) as.matrix(do.call(rbind, lapply(means.htcount.norm, function(i) i[g, 16:20])))[tissues,])
uci.as.list.mat <- lapply(unique(sig.genes$id), function(g) as.matrix(do.call(rbind, lapply(means.htcount.norm, function(i) i[g, 21:25])))[tissues,])
names(means.as.list.mat) <- unique(sig.genes$id)
names(lci.as.list.mat) <- unique(sig.genes$id)
names(uci.as.list.mat) <- unique(sig.genes$id)


pdf("DE_barplots(topright).pdf", width = 297/25.4, height = (185)/25.4, useDingbats = FALSE)

barplot(table(mf$comparison, mf$tissue)[,tissues], beside = TRUE, col = cols[1:4], bty = "l", ylab = "Number of Differentially Expressed Genes", xlab = "Tissue / Stain Derived Cohort")
legend("topright", legend = unique(mf$comparison), col = cols[1:4], fill = cols[1:4])

barplot(t(table(mf$comparison, mf$tissue)[,tissues]), beside = TRUE, col = rep(cols[1:4], each = 5), bty = "l", ylab = "Number of Differentially Expressed Genes", xlab = "Tissue / Stain Derived Cohort")

sapply(unique(as.character(sig.genes$id)), function(g) {
  barplot2(t(means.as.list.mat[[g]]), beside = TRUE,
   plot.ci = TRUE, ci.l = t(lci.as.list.mat[[g]]), ci.u = t(uci.as.list.mat[[g]]),
   ylab = paste("Normalized Read Count of", unique(sig.genes$mgi.symbol[sig.genes$id == g])),
   col = cols, xaxt = "n")
  abline(h = 0)
  axis(1, at = c(3.5, 9.5, 15.5, 21.5, 27.5), labels = c("Pituitary", "Liver", "Heart", "Embryo", "Placenta"), tick = FALSE)
  legend("topright", legend = c("B6.A-15 N2 B/B", "B6.A-17 N2 B/B", "B6.A-19 N2 B/B", "B6.A-X F1 B/B", "B6.C"), fill = cols)
})

dev.off()

all.means.as.list.mat <- lapply(unique(all.genes), function(g) as.matrix(do.call(rbind, lapply(means.htcount.norm, function(i) i[g, 11:15])))[tissues,])
all.lci.as.list.mat <- lapply(unique(all.genes), function(g) as.matrix(do.call(rbind, lapply(means.htcount.norm, function(i) i[g, 16:20])))[tissues,])
all.uci.as.list.mat <- lapply(unique(all.genes), function(g) as.matrix(do.call(rbind, lapply(means.htcount.norm, function(i) i[g, 21:25])))[tissues,])
names(all.means.as.list.mat) <- unique(all.genes)
names(all.lci.as.list.mat) <- unique(all.genes)
names(all.uci.as.list.mat) <- unique(all.genes)

pdf("AllGenes_DE_barplots.pdf", width = 297/25.4, height = 210/25.4)
sapply(all.genes, function(g) {
  barplot2(t(all.means.as.list.mat[[g]]), beside = TRUE,
   plot.ci = TRUE, ci.l = t(all.lci.as.list.mat[[g]]), ci.u = t(all.uci.as.list.mat[[g]]),
   ylab = paste("Normalized Read Count of", g),
   col = cols, xaxt = "n")
  abline(h = 0)
  axis(1, at = c(3.5, 9.5, 15.5, 21.5, 27.5), labels = c("Pituitary", "Liver", "Heart", "Embryo", "Placenta"), tick = FALSE)
  legend("topright", legend = c("B6.A-15 N2 B/B", "B6.A-17 N2 B/B", "B6.A-19 N2 B/B", "B6.A-X F1 B/B", "B6.C"), fill = cols)
})
dev.off()



## EXPERIMENTAL ##
## attempt to implement key.xtickfun = function(...) {}
# heatmap.2(t(as.matrix(exp.mat[,2:5])), dendrogram = "none", Rowv = NA,
#          Colv = NA, col = tim.colors(30), main = "",
#          margins = c(8, 2), labRow = c("B6.A-15 B/B", "B6.A-17 B/B", "B6.A-19 B/B", "B6.A-X B/B"),
#          labCol = as.character(exp.mat$id), xlab = "", cexRow = 2.8 , cexCol = .69,
#          trace = "none", lmat = matrix(c(0, 2, 3, 1, 0, 4), nrow = 2, ncol = 3),
#          lhei = c(.4, 3), lwid = c(.5, 4.5, 2),
#          key = TRUE, keysize = .5, density.info = "histogram", denscol = "black",
#          key.par = list(lwd = 2, cex.lab = 2, cex.axis = 2, mar = c(15, 0, 20, 0)), key.title = "", key.ylab = "")

