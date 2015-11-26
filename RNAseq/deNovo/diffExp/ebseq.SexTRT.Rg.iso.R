## using ebseq for matrix.* from RSEM
library(EBSeq)

load("qEB.data.rda")

## conditions
Rg1.SexTRT <- paste(qRG1$sex, qRG1$trt, sep = "_")
Rg2.SexTRT <- paste(qRG2$sex, qRG2$trt, sep = "_")
Rg3.SexTRT <- paste(qRG3$sex, qRG3$trt, sep = "_")

## contrasts
SexTRT.PosPar <- GetPatterns( Rg1.SexTRT )

## analysis for isoform usage modeling sex and brain region
Rg1.SexTRT.Iso <- EBMultiTest(Data = Rg1Iso, NgVector = IsoNgTrun,
    Conditions = Rg1.SexTRT, sizeFactors = Rg1Iso.Sizes, maxround = 8) 

Rg2.SexTRT.Iso <- EBMultiTest(Data = Rg2Iso, NgVector = IsoNgTrun,
    Conditions = Rg2.SexTRT, sizeFactors = Rg2Iso.Sizes, maxround = 8) 


Rg3.SexTRT.Iso <- EBMultiTest(Data = Rg3Iso, NgVector = IsoNgTrun,
    Conditions = Rg3.SexTRT, sizeFactors = Rg3Iso.Sizes, maxround = 8) 

save.image("qEB.SexTRT.Rg.iso.rda")


sessionInfo()
