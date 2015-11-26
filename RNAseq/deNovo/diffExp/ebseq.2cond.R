## using ebseq for matrix.* from RSEM
library(EBSeq)

## need to run ebseq.design.R before this one
load("qEB.data.rda")

## analysis of DE for two conditions : general Male vs Female independent of TRT
Sex.out <- EBTest(Data = as.matrix(GeneMat), Conditions = qDesign$sex, sizeFactors = Sizes, maxround = 8)
Sex.out$Alpha
Sex.out$Beta
Sex.out$P

Sex.PPMat <- GetPPMat(Sex.out)
## Sex.DE <- Sex.PPMat[Sex.PPMat[, "PPDE"] >= 0.95, ]  ## FDR = 0.05  post prob of diff express^ >= .95
Sex.FC <- PostFC(Sex.out, SmallNum=0.01)

## PlotPostVsRawFC(Sex.out, Sex.FC)
## layout(mat = matrix(1:2, nrow = 1))
## QQP(Sex.out)
## DenNHist(Sex.out)

## analysis of DE for two conditions : general +T vs CTRL
TRT.out <- EBTest(Data = as.matrix(GeneMat), Conditions = qDesign$trt, sizeFactors = Sizes, maxround = 8)
TRT.out$Alpha
TRT.out$Beta
TRT.out$P

save.image("qEB.2cond.rda")

TRT.PPMat <- GetPPMat(TRT.out)
## TRT.DE <- TRT.PPMat[TRT.PPMat[, "PPDE"] >= 0.95, ]  ## FDR = 0.05  post prob of diff express^ >= .95
TRT.FC <- PostFC(TRT.out, SmallNum=0.01)

## analysis of DE for region main effect : pom , vmn , tna
Rg.out <- EBMultiTest(Data = as.matrix(GeneMat), Conditions = qDesign$region, sizeFactors = Sizes, maxround = 8)
Rg.out$Alpha
Rg.out$Beta
Rg.out$P

save.image("qEB.2cond.rda")

#Rg.PPMat <- GetPPMat(Rg.out)
#Rg.DE <- Rg.PPMat[Rg.PPMat[, "PPDE"] >= 0.95, ]  ## FDR = 0.05  post prob of diff express^ >= .95
#Rg.FC <- PostFC(Rg.out, SmallNum=0.01)


## ISOFORMS

## Treatment
TRT.Iso <- EBTest(Data = IsoMat, NgVector = IsoNgTrun,
    Conditions = qDesign$trt, sizeFactors = IsoSizes, maxround = 8) 
TRT.Iso$Alpha
TRT.Iso$Beta
TRT.Iso$P

save.image("qEB.2cond.rda")

TRT.IsoPP <- GetPPMat(TRT.Iso)
## TRT.IsoDE <- TRT.IsoPP[TRT.IsoPP[, "PPDE"] >= 0.95, ]
TRT.IsoFC <- PostFC(TRT.Iso)

# Sex
Sex.Iso <- EBTest(Data = IsoMat, NgVector = IsoNgTrun,
    Conditions = qDesign$sex, sizeFactors = IsoSizes, maxround = 8) 
Sex.Iso$Alpha
Sex.Iso$Beta
Sex.Iso$P

save.image("qEB.2cond.rda")

Sex.IsoPP <- GetPPMat(Sex.Iso)
## Sex.IsoDE <- Sex.IsoPP[Sex.IsoPP[, "PPDE"] >= 0.95, ]
Sex.IsoFC <- PostFC(Sex.Iso)

# Region
Rg.Iso <- EBMultiTest(Data = IsoMat, NgVector = IsoNgTrun,
    Conditions = qDesign$region, sizeFactors = IsoSizes, maxround = 8) 
Rg.Iso$Alpha
Rg.Iso$Beta
Rg.Iso$P

save.image("qEB.2cond.rda")

Rg.IsoPP <- GetPPMat(Rg.Iso)
##Rg.IsoDE <- Rg.IsoPP[Rg.IsoPP[, "PPDE"] >= 0.95, ]
Rg.IsoFC <- PostFC(Rg.Iso)


## layout(mat = matrix(1:3, nrow =1 ))
## TRTPolyFitValue = Vector("list", 3)
## for(i in 1:3)
##    TRTPolyFitValue[[i]] <- PolyFitPlot(TRT.Iso$C1Mean[[i]], TRT.Iso$C1EstVar[[i]], 5)

save.image("qEB.2cond.rda")

sessionInfo()
