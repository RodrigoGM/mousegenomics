## using ebseq for matrix.* from RSEM
library(EBSeq)

load("qEB.data.rda")

## conditions
Rg1.SexTRT <- paste(qRG1$sex, qRG1$trt, sep = "_")
Rg2.SexTRT <- paste(qRG2$sex, qRG2$trt, sep = "_")
Rg3.SexTRT <- paste(qRG3$sex, qRG3$trt, sep = "_")

## contrasts
SexTRT.PosPar <- GetPatterns( Rg1.SexTRT )

pdf("SexTRTpattern.pdf")
PlotPattern(SexTRT.PosPar)
dev.off()

## analysis for differential expression modeling sex and testosterone treatment
Rg1.SexTRT.out <- EBMultiTest(Data = as.matrix(Rg1), Conditions = Rg1.SexTRT, sizeFactors = Rg1.Sizes, maxround = 8)
Rg2.SexTRT.out <- EBMultiTest(Data = as.matrix(Rg2), Conditions = Rg2.SexTRT, sizeFactors = Rg2.Sizes, maxround = 8)
Rg3.SexTRT.out <- EBMultiTest(Data = as.matrix(Rg3), Conditions = Rg3.SexTRT, sizeFactors = Rg3.Sizes, maxround = 8)

Rg1.SexTRT.MultiPP <- GetMultiPP(Rg1.SexTRT.out)
Rg2.SexTRT.MultiPP <- GetMultiPP(Rg2.SexTRT.out)
Rg3.SexTRT.MultiPP <- GetMultiPP(Rg3.SexTRT.out)

Rg1.SexTRT.MultiFC <- GetMultiFC(Rg1.SexTRT.out)
Rg2.SexTRT.MultiFC <- GetMultiFC(Rg2.SexTRT.out)
Rg3.SexTRT.MultiFC <- GetMultiFC(Rg3.SexTRT.out)

save.image("qEB.SexTRT.Rg.rda")

sessionInfo()

## SUMMARIES TO SEND BRET
table(Rg1.SexTRT.MultiPP$MAP)
table(Rg2.SexTRT.MultiPP$MAP)
table(Rg3.SexTRT.MultiPP$MAP)

## table of genes with statistics : fold change, p-value, fdr, etc...
patternsOI <- paste("Pattern", c(2, 4, 5, 6, 9 , 10, 13, 14), sep = "")

## RG1
Rg1TOI <- names(Rg1.SexTRT.MultiPP$MAP[Rg1.SexTRT.MultiPP$MAP %in% patternsOI])

Rg1FC <- lapply(Rg1.SexTRT.MultiFC[1:5], function(FC) FC[Rg1TOI, ])
Rg1FC$ConditionOrder <- Rg1.SexTRT.MultiFC$ConditionOrder
Rg1FC$PP <- Rg1.SexTRT.MultiPP$PP[Rg1TOI, ]
Rg1FC$MAP <- Rg1.SexTRT.MultiPP$MAP[Rg1TOI]
Rg1FC$Patterns <- Rg1.SexTRT.MultiPP$Patterns

sapply(Rg1FC, write.table, file = "Summary_Rg1FC.csv", append = TRUE, sep = ",", quote = FALSE)

## RG2
Rg2TOI <- names(Rg2.SexTRT.MultiPP$MAP[Rg2.SexTRT.MultiPP$MAP %in% patternsOI])

Rg2FC <- lapply(Rg2.SexTRT.MultiFC[1:5], function(FC) FC[Rg2TOI, ])
Rg2FC$ConditionOrder <- Rg2.SexTRT.MultiFC$ConditionOrder
Rg2FC$PP <- Rg2.SexTRT.MultiPP$PP[Rg2TOI, ]
Rg2FC$MAP <- Rg2.SexTRT.MultiPP$MAP[Rg2TOI]
Rg2FC$Patterns <- Rg2.SexTRT.MultiPP$Patterns

sapply(Rg2FC, write.table, file = "Summary_Rg2FC.csv", append = TRUE, sep = ",", quote = FALSE)

## RG3
Rg3TOI <- names(Rg3.SexTRT.MultiPP$MAP[Rg3.SexTRT.MultiPP$MAP %in% patternsOI])

Rg3FC <- lapply(Rg3.SexTRT.MultiFC[1:5], function(FC) FC[Rg3TOI, ])
Rg3FC$ConditionOrder <- Rg3.SexTRT.MultiFC$ConditionOrder
Rg3FC$PP <- Rg3.SexTRT.MultiPP$PP[Rg3TOI, ]
Rg3FC$MAP <- Rg3.SexTRT.MultiPP$MAP[Rg3TOI]
Rg3FC$Patterns <- Rg3.SexTRT.MultiPP$Patterns

sapply(Rg3FC, write.table, file = "Summary_Rg3FC.csv", append = TRUE, sep = ",", quote = FALSE)

