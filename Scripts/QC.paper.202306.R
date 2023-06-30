### This script quality controlled the protein profiles measured on SOMAscan platform (Relative Fluorescence Unit, abbreviation of "RFU").
### After required normalization, batch adjustment, sample and protein filtering, series of quality assessment were performed. 
# Date: 26 JUN 2023
# Version: 3.0
# # Author(s): name = Dr.Yun Deng
#              email = yun.deng@kennedy.ox.ac.uk
#              name = Dr.Luke jostins 
#              email = luke.jostins@kennedy.ox.ac.uk
#              name = Dr.Thomas Perry
#              email = thomas.perry@ndorms.ox.ac.uk


###-----------------------------------------------------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------------------------------------------------
### user set the file path where the folders of Resource file and Scripts are stored downloaded from Github https://github.com/dengyun-git/STEpUp_QC_Paper
userPath <- "/Users/ydeng/Documents/QCpaper.Code/" 

###-----------------------------------------------------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------------------------------------------------
library(readxl) 
library(SomaDataIO)
library(nlme)
library(sva)
library(factoextra)
library(umap)
library(ggplot2)
library(GGally)
library(mixtools)
library(data.table)
library(openxlsx)
library(writexl)
library(ggpubr)
library(gdata)
library(cowplot)

dir.create(paste0(userPath,"Robj.Paper")) ### automatically create a folder for output from this script, used to restore objects for further plots and tables 

myCodeIn <- paste0(userPath,"Scripts/")  
pathIn <- paste0(userPath,"Resource File/")
pathOut <- paste0(userPath,"Robj.Paper/")

### load in two R scripts, which contain functions to be called
source(paste0(myCodeIn,"QCnorm.Paper.202306.R"))
source(paste0(myCodeIn,"QCassess.Paper.202306.R"))

###-----------------------------------------------------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------------------------------------------------
### DATA TRANSFORMATION1: Optimized Standardization -- select the optimal normalization steps from SOMAscan routine standardization
###-----------------------------------------------------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------------------------------------------------
### Read in RFU and protein meta information from adat file. RawM -- original RFU frame. ProMeta -- protein metadata. 
RFUlist <- prepareRFU(pathIn,"SS-228545_v4.1.20220224.adat")
RawM <- RFUlist[[1]]
ProMeta <- RFUlist[[2]]

### Clear up protein meta data table by removing redundant spaces, and use consistent delimiter "|"
for(Counter in 1:nrow(ProMeta)){
  ProMeta$UniProt[Counter] <- gsub("\\s+","\\|",ProMeta$UniProt[Counter]) 
  ProMeta$EntrezGeneID[Counter] <- gsub("\\s+","\\|",ProMeta$EntrezGeneID[Counter]) 
  ProMeta$EntrezGeneSymbol[Counter] <- gsub("\\s+","\\|",ProMeta$EntrezGeneSymbol[Counter])
}

### Prepare external immunoassay measures for the comparisons with SOMAscan measures
Test1 <- c("mcp1bl","il6bl","il8bl","mmp3bl","activinabl","tsg6bl","timp1bl","tgfb1bl","fgf2bl")
Test2 <- paste0("seq.",c("2578.67","4673.13","3447.64","2788.55","13738.8","5036.50","2211.9","2333.72","3025.50"))

immunoFile <- paste0(pathIn,"Masterlist.xlsx")
sandwich_master_xls_1 <- read_excel(immunoFile,sheet=1)
sandwich_master_xls_2 <- read_excel(immunoFile,sheet=2)
sandwich_master_xls_3 <- read_excel(immunoFile,sheet=3)

### process OA data (the file name "Ben"), averaging across replicates
temp1 <- data.frame(sandwich_master_xls_2)
temp2 <- (temp1[temp1$replicate == 1,-c(1:2)] + temp1[temp1$replicate == 2,-c(1:2)])/2
sandwich_master_Ben <- data.frame(PIN=temp1[temp1$replicate == 1,1],temp2)
rownames(sandwich_master_Ben) <- sandwich_master_Ben$PIN

#process injury data (the file name "Historic")
sandwich_master_Historic <- data.frame(sandwich_master_xls_3)
rownames(sandwich_master_Historic) <- sandwich_master_Historic$PIN

###-----------------------------------------------------------------------------------------------------------------------------------------
### Investigate different combinations of SOMAscan routine normalisation steps, and for each combination we store the corresponding evaluation criteria objects for further convenient reference
inputList <- list(c("RawM"),c("HYBNORM"),c("HYBNORM","MIDNORMcali","PLATESCALE"),c("HYBNORM","MIDNORMcali","PLATESCALE","MIDNORMsamp"),c("HYBNORM","MIDNORMcali","PLATESCALE","MIDNORMsamp","CALIBRATION"),c("HYBNORM","MIDNORMcali","PLATESCALE","CALIBRATION"))

for(UserInput in inputList){
  ### call getMySoma funnction to perform user selected standardization steps
  MySoma2022 = getMySoma(UserInput,RawM)
  exprDat <- MySoma2022[grep("STEP",MySoma2022$SampleId),which(colnames(MySoma2022)=="seq.10000.28"):ncol(MySoma2022)]
  rownames(exprDat) <- MySoma2022[grep("STEP",MySoma2022$SampleId),"SampleId"]
  
  ### Checks against calibrators.
  calib_norm <- as.matrix(MySoma2022[grep("POOL",MySoma2022$SampleId),-c(1:(which(colnames(MySoma2022)=="seq.10000.28")-1))])
  calibIDs <-  MySoma2022$SampleId[grep("POOL",MySoma2022$SampleId)]  ###calibID corresponding to calib_norm
  
  CVoa <- CalibratorCheck2(MySoma2022,"OA",exprDat,"OA")
  CVinj <- CalibratorCheck2(MySoma2022,"INJ",exprDat,"Injury")
  
  ### pool biological variance explained
  R2_norm1 = VarExp2(MySoma2022,"OA",exprDat,"")
  R2_norm2 = VarExp2(MySoma2022,"INJ",exprDat,"")
  
  CorData_norm1 <- ExtVal(sandwich_master_Ben,exprDat,Test1,Test2)
  CorData_norm2 <- ExtVal(sandwich_master_Historic,exprDat,Test1,Test2)
  
  ### for further comparisons across different standardization, we restore "CVoa","CVinj","R2_norm1","R2_norm2","CorData_norm1","CorData_norm2" by renaming them with distinguishable suffix in separate files
  l1 <- list(CVoa,CVinj,R2_norm1,R2_norm2,CorData_norm1,CorData_norm2)
  l2 <- c("CVoa","CVinj","R2_norm1","R2_norm2","CorData_norm1","CorData_norm2")
  for(VRct in 1:6){
    VR <- l1[[VRct]]
    gdata::mv(from="VR",to=paste(l2[VRct],paste(UserInput,collapse="."),sep="."))
  }
  
  DataTransformation1 <- sapply(paste(c("CVoa","CVinj","R2_norm1","R2_norm2","CorData_norm1","CorData_norm2"),paste(UserInput,collapse="."),sep="."),function(x){get(x)})
  save(DataTransformation1,file=paste0(pathOut,paste(UserInput,collapse="."),".DataTransformation1.Rdat"))
}


### After storing all the "DataTransformation1.Rdat" for different standardizations as above, it is convenient to directly read in all the CV/R2/CorData related R objects
for(x in list.files(pathOut,".DataTransformation1.Rdat")){
  load(paste0(pathOut,x))
  for(VRct in 1:6){
    VR <- DataTransformation1[[VRct]]
    gdata::mv(from="VR",to=names(DataTransformation1)[VRct])
  }
}

###-----------------------------------------------------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------------------------------------------------
### DATA TRANSFORMATION2: Batch correction and IPS adjustment
###-----------------------------------------------------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------------------------------------------------
### We defined our optimal standardization as "inputList[[6]]"
UserInput <- inputList[[6]]
MySoma2022 = getMySoma(UserInput,RawM)
exprDat <- MySoma2022[grep("STEP",MySoma2022$SampleId),which(colnames(MySoma2022)=="seq.10000.28"):ncol(MySoma2022)]
rownames(exprDat) <- MySoma2022[grep("STEP",MySoma2022$SampleId),"SampleId"]

### Combine clinical meta data with standardised MySoma, merge key as STEpUpId.
### MySoma and MySoma2022 are the same, only difference is: MySoma with  STEpUpId, for the convenience to match clinical meta data. MySoma2022 maintains "SampleId" as column name.
FrarmeList1 <- getFrarmeList(pathIn,MySoma2022,ProMeta)
MySoma <- FrarmeList1[[1]]
CombinedFrame <- FrarmeList1[[2]]
exprDat_norm <- FrarmeList1[[3]]
exprDat_norm.HM <- FrarmeList1[[4]]

### Dimensionality reduction by PCA, 80% variations are kept by top PCs based on standardized RFU
pcDat <- pcDat.standardised <- getTopPCdt(exprDat_norm.HM,80)[[1]]

### Investigate associations between predefined technical confounders and top PCs 
ConfouderCheck2 <- ConfouderCheck(CombinedFrame,pcDat.standardised)
ConfounderTable2 <- ConfouderCheck2[[1]]

### detect bimodal signal
dev2PC <- rownames(ConfounderTable2)[order(ConfounderTable2[,"sf_iknee_proc_batch"])[1:2]]
gmm_fit <- mixtools::normalmixEM(pcDat.standardised[,dev2PC[1]])
plot(gmm_fit,which=2,breaks=100,xlab2=dev2PC[1],main2="Bimodal Signal Distribution along PC2",cex.main=0.9)
bimodalLabel <- gmm_fit$posterior[,1] > 0.5
bimodalLabel <- ifelse(bimodalLabel==TRUE,"bimodal1","bimodal2")
names(bimodalLabel)=CombinedFrame$STEpUpId

### Perform batch correction -- obtained exprDat.batchDone
BimodalBatch <- as.factor(bimodalLabel)
PlateId <- levels(as.factor(CombinedFrame$PlateId))
PlateBimodal <- unlist(sapply(1:length(BimodalBatch),function(x){ifelse(BimodalBatch[x]=="bimodal2",which(PlateId==CombinedFrame$PlateId[x])+22,which(PlateId==CombinedFrame$PlateId[x]))}))
exprDat.batchDone <- exp(t(sva::ComBat(t(log(exprDat_norm)), PlateBimodal, mod=NULL, par.prior = TRUE, prior.plots = FALSE)))

## Extract out spun/unspun paired samples to investigate PC1 driver
spun_sams <- grep("-SP",MySoma$STEpUpId,value=T)
unspun_sams <- grep("-UN",MySoma$STEpUpId,value=T)
spun_unspun_f <- CombinedFrame[c(spun_sams,unspun_sams),]
exp_spun_unspun <- spun_unspun_f[,which(colnames(spun_unspun_f)=="seq.10000.28"):ncol(spun_unspun_f)]
names(exp_spun_unspun) <- colnames(exprDat_norm)
spun_unspun_label <- c(rep(1,length(spun_sams)),rep(0,length(unspun_sams))) ### 1: spun; 0:unspun

rowmap <- 1:length(MySoma$PlateId)
names(rowmap) <- MySoma$STEpUpId
MySoma_spun <- MySoma[rowmap[spun_sams],]
MySoma_unspun <- MySoma[rowmap[unspun_sams],]

MySoma_spun_unspun <- MySoma[rowmap[c(spun_sams,unspun_sams)],]
spun_unspun_d <- sapply(which(colnames(MySoma_spun)=="seq.10000.28"):ncol(MySoma_spun),function(i) mean(log(MySoma_unspun[,i])-log(MySoma_spun[,i]))/sd(log(c(MySoma_spun[,i],MySoma_unspun[,i]))))
names(spun_unspun_d) <- colnames(exprDat_norm)

### Intracellular protein score calculation using 18 paired spun/unspun samples
IPScore <- apply(t(t(log(exprDat_norm))*spun_unspun_d),1,sum) # we used transpose matrix to avoid the errors of element arrangement orders brought by row vector calculation  

### Perform intracellular protein score adjustment and batch correction, obtained exprDat.batchDone2
exprDat.RegDone <- exp(t(limma::removeBatchEffect(t(log(exprDat_norm)), covariates=IPScore)))
exprDat.batchDone2 <- exp(t(sva::ComBat(t(log(exprDat.RegDone)), PlateBimodal, mod=NULL, par.prior = TRUE, prior.plots = FALSE)))

### Dimensionality reduction by PCA, 80% variations are kept by top PCs based on batch corrected/non-IPS adjusted data, and batch corrected and IPS adjusted data 
tempList <- getTopPCdt(exprDat.batchDone,80)
pcDat.batchDone <- tempList[[1]]
eig.val.batchDone <- tempList[[2]]
tempList <- getTopPCdt(exprDat.batchDone2,80)
pcDat.batchDone2 <- tempList[[1]]
eig.val.batchDone2 <- tempList[[2]]

###-----------------------------------------------------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------------------------------------------------
### DATA FILTERING (Sample & Protein Filtering )
###-----------------------------------------------------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------------------------------------------------
### Sample filtering: PCA outliers
centroid <- colMeans(pcDat)
distanceToCenter = vector(mode="numeric",length=nrow(pcDat))
names(distanceToCenter) = rownames(pcDat)
for (mK in 1: nrow(pcDat)){
  distanceToCenter[mK] <- sqrt(sum((pcDat[mK,]-centroid)^2))
}
PCAout <- which(distanceToCenter>(mean(distanceToCenter) + 5*sd(distanceToCenter)))
RemoveS1 = names(PCAout)

### Sample filtering: total protein outliers
totalProtein_norm = TotalProCheck(CombinedFrame,exprDat_norm.HM,"Standardised")
TotalProteinOut <- names(totalProtein_norm)[which(totalProtein_norm>(mean(totalProtein_norm) +5*sd(totalProtein_norm)))]
RemoveS2 = TotalProteinOut

### Sample filtering: limit of detection -- below LoD per protein or above uLoD
CompM1 = LoDdetection(MySoma) ###CompM1 < 0, RFU below limit of Detection
CompM2 = uLoDdetection(MySoma) ###CompM2 > 0, RFU above upper limit of Detection
SampRatio1 = apply(CompM1,1,function(x){length(which(x<0))/ncol(CompM1)})
SampRatio2 = apply(CompM2,1,function(x){length(which(x>0))/ncol(CompM2)})
RemoveS3 = c(MySoma[names(which(SampRatio1>0.25)),"STEpUpId"],MySoma[names(which(SampRatio2>0.25)),"STEpUpId"])

### Somascan inhouse filter Hybridization Scale Factor > 2.5
RemoveS4F = prepareRFU(pathIn,"SS-228545_v4.1.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.20220224.adat")[[1]]
# all(rownames(MySoma)==rownames(RemoveS4F)) ### always check the orders of data frames which should be matching each other 
RemoveS4 = MySoma$STEpUpId[which(RemoveS4F$HybControlNormScale>2.5)]

### Protein filterinng: non human proteins
NonHuman <- rownames(ProMeta)[which(ProMeta$Organism !="Human" | ProMeta$Type !="Protein")]

### Protein filtering: proteins majorly with technical variations based on OA/Injury Pooled sample replicates
R2_norm1 = VarExp2(MySoma2022,"OA",exprDat_norm.HM,"")
R2_norm2 = VarExp2(MySoma2022,"INJ",exprDat_norm.HM,"")
RemovePro1 = names(which(R2_norm1<0.5))
RemovePro2 = names(which(R2_norm2<0.5))

### Protein filtering: significantly associated with technical confounders. Batch corrected/non-IPS adjusted, and batch corrected/IPS adjusted are investigated seperately 
################################################################################################################
### Batch corrected data
bimodalP1<- bimodalTest(exprDat.batchDone,bimodalLabel) # association between proteins and bimodal signal 
###all(rownames(CombinedFrame)==rownames(exprDat.batchDone))
ConfouderCheck1.batchDone <- ConfouderCheck(CombinedFrame,exprDat.batchDone) # associationn between proteins and technical confounders
ConfounderTable1.batchDone <- ConfouderCheck1.batchDone[[1]]

### Table summarizing the filters for batch corrected data
removeTableList.BatchCorrected <- getRemoveTable(NonHuman,RemoveS1,RemoveS2,RemoveS3,RemoveS4,RemovePro1,RemovePro2,bimodalP1,ConfounderTable1.batchDone[,c(7,10)],exprDat.batchDone)
removeTable.BatchCorrected <- removeTableList.BatchCorrected[[1]]
exprDatFiltered.BatchCorrected <- removeTableList.BatchCorrected[[2]]

tempList <- getTopPCdt(exprDatFiltered.BatchCorrected,80)
pcDat.batchCorrected.Filtered <- tempList[[1]]
eig.val.batchCorrected.Filtered <- tempList[[2]]

################################################################################################################
### Intracellular protein score adjusted, batch corrected data
bimodalP2<- bimodalTest(exprDat.batchDone2,bimodalLabel) # association between proteins and bimodal signal 
###all(rownames(CombinedFrame)==rownames(exprDat.batchDone2))
ConfouderCheck1.batchDone2 <- ConfouderCheck(CombinedFrame,exprDat.batchDone2) # associationn between proteins and technical confounders
ConfounderTable1.batchDone2 <- ConfouderCheck1.batchDone2[[1]]

### Table summarizing the filters for IPS adjusted, batch corrected data
removeTableList.IPSregressed <- getRemoveTable(NonHuman,RemoveS1,RemoveS2,RemoveS3,RemoveS4,RemovePro1,RemovePro2,bimodalP2,ConfounderTable1.batchDone2[,c(7,10)],exprDat.batchDone2)
removeTable.IPSregressed <- removeTableList.IPSregressed[[1]]
exprDatFiltered.IPSregressed <- removeTableList.IPSregressed[[2]]

tempList <- getTopPCdt(exprDatFiltered.IPSregressed,80)
pcDat.IPSreg.Filtered <- tempList[[1]]
eig.val.IPSreg.Filtered <- tempList[[2]]

tempList <- getTopPCdt(exprDat_norm,80)
pcDat.standardised.all <- tempList[[1]]
eig.val.standardised.all <- tempList[[2]]

### TableS4 in QC paper shows removeTable.BatchCorrected and removeTable.IPSregressed 

### store all the RFUs, top PCs(80% variation explained), confounder table at different stages for the convenience of future references
save(ProMeta, MySoma2022, MySoma, CombinedFrame, exprDat_norm,exprDat.batchDone,exprDat.batchDone2,exprDatFiltered.BatchCorrected,exprDatFiltered.IPSregressed,file=paste0(pathOut,"All.RUFs.Rdat"))
save(pcDat.standardised,pcDat.batchDone,pcDat.batchDone2,pcDat.batchCorrected.Filtered,pcDat.IPSreg.Filtered,pcDat.standardised.all,eig.val.batchDone,eig.val.batchDone2,eig.val.batchCorrected.Filtered,eig.val.IPSreg.Filtered,eig.val.standardised.all,file=paste0(pathOut,"pcDat.Rdat"))
save(NonHuman,RemoveS1,RemoveS2,RemoveS3,RemoveS4,RemovePro1,RemovePro2,bimodalP1,ConfounderTable1.batchDone,bimodalP2,ConfounderTable1.batchDone2,file=paste0(pathOut,"Filters.Rdat"))

mn###-----------------------------------------------------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------------------------------------------------
### Quality assessment and storing R objects for visualization by plots and tables in paper, which run in QC.plot.202306.R
###-----------------------------------------------------------------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------------------------------------------------------------
################################################################################################################
### 1. Comparisons for different standardizations 
################################################################################################################
### construct variables for convenient labelling in plots
AllStandardisation <- c("RawM","HYBNORM","HYBNORM.MIDNORMcali.PLATESCALE","HYBNORM.MIDNORMcali.PLATESCALE.MIDNORMsamp","HYBNORM.MIDNORMcali.PLATESCALE.MIDNORMsamp.CALIBRATION","HYBNORM.MIDNORMcali.PLATESCALE.CALIBRATION")
StandardisationLabel = c("Original RFU","Hybridisation Normalisation","Hybridisation Normalisation+\nPlate Scaling","Hybridisation Normalisation+\nPlate Scaling+\nMedian Normalsation on SF","Hybridisation Normalisation+\nPlate Scaling+\nMedian Normalsation on SF+\nPlate Calibration","Hybridisation Normalisation+\nPlate Scaling+\nPlate Calibration")

tempV<- vector("numeric",length=length(AllStandardisation))
names(tempV) <- StandardisationLabel
meanCorOA <- meanCorINJ <- meanCV.OA <- meanCV.INJ <- meanR2.OA <- meanR2.INJ <- tempV
for(x in 1:length(AllStandardisation)){
  meanCV.OA[x] <- signif(mean(get(paste0("CVoa.",AllStandardisation[x]))),4)
  meanCV.INJ[x] <- signif(mean(get(paste0("CVinj.",AllStandardisation[x]))),4)
  meanR2.OA[x] <- signif(mean(get(paste0("R2_norm1.",AllStandardisation[x]))),4)
  meanR2.INJ[x] <- signif(mean(get(paste0("R2_norm2.",AllStandardisation[x]))),4)
  meanCorOA[x] <- signif(mean(as.numeric(na.omit(get(paste0("CorData_norm1.",AllStandardisation[x]))$cor))),4)
  meanCorINJ[x] <- signif(mean(as.numeric(na.omit(get(paste0("CorData_norm2.",AllStandardisation[x]))$cor))),4)
}

### Mean CV
meanCV <- data.frame(cbind(rep(c("OA","injury"),each=length(StandardisationLabel)),rep(1:length(AllStandardisation),2),c(meanCV.OA,meanCV.INJ)))
colnames(meanCV) <- c("DiseaseGroup","NormalisationSteps","meanCV")

### MeanR2
meanR2 <- data.frame(cbind(rep(c("OA","injury"),each=length(StandardisationLabel)),rep(1:length(StandardisationLabel),2),c(meanR2.OA,meanR2.INJ)))
colnames(meanR2) <- c("DiseaseGroup","NormalisationSteps","meanR2")

### Correlation coefficients between SOMAscan and external immunoassay
externalCheckOA <- externalCheckINJ <- vector()
for(x in 1:length(AllStandardisation)){
  externalCheckOA <- cbind(externalCheckOA,signif(as.numeric(get(paste0("CorData_norm1.",AllStandardisation[x]))$cor,4)))
  externalCheckINJ <- cbind(externalCheckINJ,signif(as.numeric(get(paste0("CorData_norm2.",AllStandardisation[x]))$cor,4)))
}
externalCheckOA <- cbind(CorData_norm1.RawM[,1:2],externalCheckOA)
externalCheckINJ <- cbind(CorData_norm2.RawM[,1:2],externalCheckINJ)
externalCheckOA[,"SandwichName"] <- externalCheckINJ[,"SandwichName"] <-  c("MCP1","IL6","IL8","MMP3","Activin A","TSG6","TIMP1","TGFb1","FGF2")

CorDatY = as.matrix(externalCheckOA[,3],ncol=1)
CorDatX = rep(as.matrix(externalCheckOA[,1],ncol=1),(ncol(externalCheckOA)-2))
CorC = rep(seq(1:(ncol(externalCheckOA)-2)),rep(nrow(externalCheckOA),ncol(externalCheckOA)-2))

for(k in 4:ncol(externalCheckOA)){
  CorDatY =  rbind(CorDatY,as.matrix(externalCheckOA[,k],ncol=1))
}
colnames(CorDatY) = "CorDatY"
CorDatP.OA = as.data.frame(cbind(CorDatX,CorC,CorDatY))

CorDatY = as.matrix(externalCheckINJ[,3],ncol=1)
CorDatX = rep(as.matrix(externalCheckINJ[,1],ncol=1),(ncol(externalCheckINJ)-2))
CorC = rep(seq(1:(ncol(externalCheckINJ)-2)),rep(nrow(externalCheckINJ),ncol(externalCheckINJ)-2))

for(k in 4:ncol(externalCheckINJ)){
  CorDatY =  rbind(CorDatY,as.matrix(externalCheckINJ[,k],ncol=1))
}
colnames(CorDatY) = "CorDatY"
CorDatP.INJ = as.data.frame(cbind(CorDatX,CorC,CorDatY))

save(meanCV,meanR2,CorDatP.OA,CorDatP.INJ,file=paste0(pathOut,"Figure1.Rdat"))

################################################################################################################
### 2. PC1 driver
################################################################################################################
### PCA space for all the different RUFs (80% variation explained)
### Pearson correlation coefficients between each protein and PC1
corPerPro.BatchCorrected <- getCorPro(exprDat.batchDone,ProMeta,pcDat.batchDone)
corPerPro.reg <- getCorPro(exprDat.batchDone2,ProMeta,pcDat.batchDone2)

### PC1 of 18 paired spun/unspun samples 
unspun_sams <- sort(grep("-UN",MySoma$STEpUpId,value=T))
spun_sams <- sort(grep("-SP",MySoma$STEpUpId,value=T))

spin.col <- sapply(1:18,function(x){ifelse(pcDat.batchDone[spun_sams[x],1] > pcDat.batchDone[unspun_sams[x],1],"green","blue")})
spinFrame <-data.frame(cbind(c(rep("Unspun",length(unspun_sams)),rep("Spun",length(spun_sams))),c(pcDat.batchDone[unspun_sams,1],pcDat.batchDone[spun_sams,1])))
colnames(spinFrame) <- c("Spin","PC1")

### PCA object for standardized,  batch corrected and IPS adjusted data
pcs_standardised <- prcomp(log10(exprDat_norm),scale=TRUE)
pcs_batchDone2 <- prcomp(log10(exprDat.batchDone2),scale=TRUE)

save(IPScore,corPerPro.BatchCorrected,corPerPro.reg,spinFrame,spin.col,eig.val.standardised.all,eig.val.batchDone2,file=paste0(pathOut,"Figure2.Rdat"))

################################################################################################################
### 3. PC2 driver
################################################################################################################
### The strong bimdal signal marker "TSG101"
TSGseq <- rownames(ProMeta)[grep("TSG101",ProMeta$EntrezGeneSymbol)]
keepID <- which(CombinedFrame$sf_iknee_proc_order>0)
CombinedFrameKeep <- CombinedFrame[keepID,]
exprDat_normKeep <- exprDat_norm[keepID,]

STEP1800IDs <- c(grep("STEP1800-F-V1-HT1",MySoma2022$SampleId),grep("STEP1800-F-V1-A-HT1",MySoma2022$SampleId))
STEP1801IDs <- c(grep("STEP1801-F-V1-HT1",MySoma2022$SampleId),grep("STEP1801-F-V1-A-HT1",MySoma2022$SampleId))
STEP1807IDs <- c(grep("STEP1807-F-V1-HT1",MySoma2022$SampleId),grep("STEP1807-F-V1-A-HT1",MySoma2022$SampleId))

reprocessFrame <- data.frame(cbind(c(rep("Original",3),rep("Reprocessed",3)),rep(c("STEP1800","STEP1801","STEP1807"),2),c(MySoma[c(STEP1800IDs[1],STEP1801IDs[1],STEP1807IDs[1]),TSGseq],MySoma[c(STEP1800IDs[2],STEP1801IDs[2],STEP1807IDs[2]),TSGseq])))
colnames(reprocessFrame) <- c("Processing","Sample","TSG101") 

### UMAP visualization on the data before and after bimodal batch correction
myUmap1 = umap::umap(exprDat_norm,random_state=123) #set seeds for UMAP
myUmap2 = umap::umap(exprDat.batchDone,random_state=123)

save(gmm_fit,bimodalLabel,myUmap1,myUmap2,TSGseq,CombinedFrameKeep,exprDat_normKeep,reprocessFrame,file=paste0(pathOut,"Figure3.Rdat"))

################################################################################################################
### 4. External check with immunoassay for IPS adjustment and batch correction
################################################################################################################### 
### Correlation with immunoassay for batch corrected non-IPS adjusted and IPS adjusted data
CorData_norm1.noLysis.done <- ExtVal(sandwich_master_Ben,exprDat.batchDone,Test1,Test2)
CorData_norm2.noLysis.done <- ExtVal(sandwich_master_Historic,exprDat.batchDone,Test1,Test2)

CorData_norm1.Lysis.done <- ExtVal(sandwich_master_Ben,exprDat.batchDone2,Test1,Test2)
CorData_norm2.Lysis.done <- ExtVal(sandwich_master_Historic,exprDat.batchDone2,Test1,Test2)

AllStandardisation <- c("RawM","HYBNORM.MIDNORMcali.PLATESCALE.CALIBRATION","noLysis.done","Lysis.done")
StandardisationLabel = c("Raw Data","Optimized Standardisation","Processed without IPS Adjustment","Processed with IPS Adjustment")

externalCheckOA.IPS <- externalCheckINJ.IPS <- vector()
for(x in 1:length(AllStandardisation)){
  externalCheckOA.IPS <- cbind(externalCheckOA.IPS,signif(as.numeric(get(paste0("CorData_norm1.",AllStandardisation[x]))$cor,4)))
  externalCheckINJ.IPS <- cbind(externalCheckINJ.IPS,signif(as.numeric(get(paste0("CorData_norm2.",AllStandardisation[x]))$cor,4)))
}
externalCheckOA.IPS <- cbind(CorData_norm1.RawM[,1:2],externalCheckOA.IPS)
externalCheckINJ.IPS <- cbind(CorData_norm2.RawM[,1:2],externalCheckINJ.IPS)
externalCheckOA.IPS[,"SandwichName"] <- externalCheckINJ.IPS[,"SandwichName"] <-  c("MCP1","IL6","IL8","MMP3","Activin A","TSG6","TIMP1","TGFb1","FGF2")
colnames(externalCheckOA.IPS)[3:ncol(externalCheckOA.IPS)] <- colnames(externalCheckINJ.IPS)[3:ncol(externalCheckINJ.IPS)] <- StandardisationLabel

CorDatP.OA.IPS <- plotExternal9(externalCheckOA.IPS)
CorDatP.INJ.IPS <- plotExternal9(externalCheckINJ.IPS)

save(CorDatP.OA.IPS,CorDatP.INJ.IPS,sandwich_master_Ben,sandwich_master_Historic,Test1,file=paste0(pathOut,"Figure4.Rdat"))

### Correlation between the 9 immunoassay proteins and IPS score. As shown in TableS6
KeepOA <- unlist(sapply(sandwich_master_Ben$PIN,function(x){names(IPScore)[grep(x,names(IPScore))]}))
sapply(Test1,function(x) {cor.test(as.numeric(sandwich_master_Ben[,x]),IPScore[KeepOA],use="pairwise.complete.obs")$estimate})
sapply(Test1,function(x) {cor.test(as.numeric(sandwich_master_Ben[,x]),IPScore[KeepOA],use="pairwise.complete.obs")$p.value})

KeepINJ <- unlist(sapply(sandwich_master_Historic$PIN,function(x){names(IPScore)[grep(x,names(IPScore))]}))
sapply(Test1[-6],function(x) {cor.test(as.numeric(sandwich_master_Historic[,x]),IPScore[KeepINJ],use="pairwise.complete.obs")$estimate})
sapply(Test1[-6],function(x) {cor.test(as.numeric(sandwich_master_Historic[,x]),IPScore[KeepINJ],use="pairwise.complete.obs")$p.value})

################################################################################################################
### 5. Technical confounders association with top 10 PCs 
################################################################################################################
### for standardised data
ConfouderCheck2.Standerdised <- ConfouderCheck(CombinedFrame,pcDat.standardised.all)
ConfounderTable2.Standerdised <- ConfouderCheck2.Standerdised[[1]]
outEXL(pathOut,"ConfounderTable2.Standerdised",signif(ConfounderTable2.Standerdised,4))

### for batch corrected, non-IPS adjusted, non-filtered data
ConfouderCheck2.batchCorrected <- ConfouderCheck(CombinedFrame[rownames(pcDat.batchDone),],pcDat.batchDone)
ConfounderTable2.batchCorrected <- ConfouderCheck2.batchCorrected[[1]]
outEXL(pathOut,"ConfounderTable2.batchCorrected",signif(ConfounderTable2.batchCorrected,4))

### for batch corrected, IPS adjusted, non-filtered data
ConfouderCheck2.IPSreg <- ConfouderCheck(CombinedFrame[rownames(pcDat.batchDone2),],pcDat.batchDone2)
ConfounderTable2.IPSreg <- ConfouderCheck2.IPSreg[[1]]
outEXL(pathOut,"ConfounderTable2.IPSreg",signif(ConfounderTable2.IPSreg,4))

### for batch corrected, non-IPS adjusted, filtered data
ConfouderCheck2.batchCorrected.Filtered <- ConfouderCheck(CombinedFrame[rownames(pcDat.batchCorrected.Filtered),],pcDat.batchCorrected.Filtered)
ConfounderTable2.batchCorrected.Filtered <- ConfouderCheck2.batchCorrected.Filtered[[1]]
outEXL(pathOut,"ConfounderTable2.BC.Filtered",signif(ConfounderTable2.batchCorrected.Filtered,4))

### for batch corrected, IPS adjusted, filtered data
ConfouderCheck2.IPSreg.Filtered <- ConfouderCheck(CombinedFrame[rownames(pcDat.IPSreg.Filtered),],pcDat.IPSreg.Filtered)
ConfounderTable2.IPSreg.Filtered <- ConfouderCheck2.IPSreg.Filtered[[1]]
outEXL(pathOut,"ConfounderTable2.IPS.Filtered",signif(ConfounderTable2.IPSreg.Filtered,4))

save(ConfounderTable2.Standerdised,ConfounderTable2.batchCorrected,ConfounderTable2.IPSreg,ConfounderTable2.batchCorrected.Filtered,ConfounderTable2.IPSreg.Filtered,file=paste0(pathOut,"TechVS10pc.Rdat"))

################################################################################################################
### 6. UMAP for data after filtering, on IPS adjusted/batch corrected data, and batch corrected/non-IPS adjusted data
################################################################################################################ 
myUmap.BC.final = umap::umap(pcDat.batchCorrected.Filtered,random_state=123) # set seeds for UMAP
myUmap.IPS.final = umap::umap(pcDat.IPSreg.Filtered,random_state=123)
CombinedFrame1 <- CombinedFrame[rownames(pcDat.batchCorrected.Filtered),]
CombinedFrame2 <- CombinedFrame[rownames(pcDat.IPSreg.Filtered),]
keepID1 <- which(!is.na(CombinedFrame1$sf_iknee_qc_group))
keepID2 <- which(!is.na(CombinedFrame2$sf_iknee_qc_group))
myUmap.BC.final.F <- data.frame(cbind(myUmap.BC.final$layout[keepID1,],CombinedFrame1$sf_iknee_qc_group[keepID1]))
myUmap.IPS.final.F <- data.frame(cbind(myUmap.IPS.final$layout[keepID2,],CombinedFrame2$sf_iknee_qc_group[keepID2]))
colnames(myUmap.BC.final.F) <- colnames(myUmap.IPS.final.F) <- c("D1","D2","DiseaseGroup")
save(myUmap.BC.final.F,myUmap.IPS.final.F,file=paste0(pathOut,"UmapDis.filtered.Rdat"))

################################################################################################################
### 7. %CV distribution based on pooled sample replicates stratified by disease group
################################################################################################################
OA.CV <- CVbreak(calib_norm,calibIDs,"OA")
INJ.CV <- CVbreak(calib_norm,calibIDs,"INJ")
  
save(OA.CV,INJ.CV,file=paste0(pathOut,"CV.OA.Repeats.FreezeThaw.Reprocess.Rdat"))

