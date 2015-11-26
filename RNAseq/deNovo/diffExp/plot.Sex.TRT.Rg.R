## Plotting results from DE analysis

## libraries
library(EBSeq)
library(ggplot2)
library(hexbin)

## data
load("qEB.Sex.TRT.Rg.rda")
plot(PPDE ~ Log2PostFC, aComp, pch = 19, cex = 0.3)
plot(Log2PostFC ~ log2(TxMean), aComp, pwd = ".")

pdf("DE_Summary.pdf", width = 420/25.4, height = 297/25.4)
## Volcano Plot with all the data
gg <- ggplot(aComp, aes(x = Log2PostFC, y = PPDE, color = Comparison, fill = Comparison))
## hex-binning, shows mutiple observations per hexagon
gg + stat_binhex(bins = 400)  + facet_grid(Region ~ Comparison) + geom_hline(yintercept = c(0.95, 0.99), color = "black", linetype = 2, size = 0.2) + xlim(-10, 10) + geom_vline(xintercept = c(-5, 5), linetype = 2, size = 0.2) + ggtitle("HexBin Volcano Plot - PPDE vs. Log2PostFC ( all transcripts)") + xlab(expression(italic(Log[2]~PostFC)))

## M-A plot with all the data
gg <- ggplot(aComp, aes(x = 1/2 * (log2(C1Mean) + log2(C2Mean)), y = Log2PostFC, color = Comparison, fill = Comparison))
## Use hex-bin for overplotting
gg + stat_binhex(bins = 400)  + facet_grid(Region ~ Comparison) + ggtitle("HexBin MA Plot - Log2PostFC vs. Log2 Transcript Mean (all transcripts)") + ylab(expression(italic(Log[2])~PostFC)) + xlab(expression(over(1,2)~(~italic(Log[2]~(C1Mean)) + italic(Log[2]~(C2Mean)))))

## Volcano plot with significant genes
gg <- ggplot(pwde, aes(x = Log2PostFC, y = PPDE, color = Comparison, fill = Comparison))
## hex-binning, shows mutiple observations per hexagon
gg + stat_binhex(bins = 400) + facet_grid(Region ~ Comparison)  + geom_hline(yintercept = c(0.95, 0.99), color = "black", linetype = 2, size = 0.2) + xlim(-10, 10) + geom_vline(xintercept = c(-5, 5), linetype = 2, size = 0.2, color = "black") + ggtitle("HexBin Volcano Plot - PPDE vs. Log2PostFC (FDR < 0.05)") + xlab(expression(italic(Log[2]~PostFC)))

## M-A plot with significant genes
gg <- ggplot(pwde, aes(x = 1/2 * (log2(C1Mean) + log2(C2Mean)), y = Log2PostFC, color = Comparison, fill = Comparison))
## Use hex-bin for overplotting
gg + stat_binhex(bins = 400)  + facet_grid(Region ~ Comparison) + ggtitle("HexBin MA Plot - Log2PostFC vs. Log2 Transcript Mean (FDR < 0.05)") + geom_hline(yintercept = 0) + ylab(expression(italic(Log[2])~PostFC)) + xlab(expression(over(1,2)~(~italic(Log[2]~(C1Mean)) + italic(Log[2]~(C2Mean)))))

dev.off()



## single point per observation
## gg + geom_point()  + facet_grid(Region ~ Comparison) + geom_hline(yintercept = c(0.95, 0.99), color = "black", linetype = 2, size = 0.2) + xlim(-10, 10) + geom_vline(xintercept = c(-5, 5), linetype = 2, size = 0.2, color = "black") + ggtitle("Volcano Plot - PPDE vs. Log2PostFC (FDR < 0.05)")


pdf("QQPlots_PairwiseComp.pdf", width = 297/25.4, height = 148/25.4)
layout(mat = matrix(1:2, ncol = 2))
sapply(list(Rg1.SMvsSF, Rg1.TMvsTF, Rg1.SMvsTM, Rg1.SFvsTF, Rg2.SMvsSF, Rg2.TMvsTF, Rg2.SMvsTM, Rg2.SFvsTF, Rg3.SMvsSF, Rg3.TMvsTF, Rg3.SMvsTM, Rg3.SFvsTF), function(PWDE) {QQP(PWDE)})
dev.off()


pdf("DensityHist_PairwiseComp.pdf", width = 297/25.4, height = 148/25.4)
layout(mat = matrix(1:2, ncol = 2))
sapply(list(Rg1.SMvsSF, Rg1.TMvsTF, Rg1.SMvsTM, Rg1.SFvsTF, Rg2.SMvsSF, Rg2.TMvsTF, Rg2.SMvsTM, Rg2.SFvsTF, Rg3.SMvsSF, Rg3.TMvsTF, Rg3.SMvsTM, Rg3.SFvsTF), function(PWDE) {DenNHist(PWDE)})
dev.off()
