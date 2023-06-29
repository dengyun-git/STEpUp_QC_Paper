# This script contains functions for quality assessment on the investigated RFU. To be called by QC.paper.202306.R.
# Version: 3.0
# Author(s): name = Dr.Yun Deng
#            email = yun.deng@kennedy.ox.ac.uk
#            name = Dr.Luke jostins 
#            email = luke.jostins@kennedy.ox.ac.uk
#            name = Dr.Thomas Perry
#            email = thomas.perry@ndorms.ox.ac.uk

#__________________________________________________________________________________________________________
# function to check %CV based on synovial fluid pooled samples.
# input: calib_norm -- concentration frame for pooled samples and QC plasma calibrators; calibIDs: sample ID for calibrators; clinicType: disease group.
#        MySoma2022 -- RFU frame; exprDat -- protein concentration frame.
# output: temp3 -- a vector of CV for interplate, intraplate, spun/unspun, freeze thaw, different hyaluronidase treatments.
#         plots of accumulative distribution of CVs marked at 80 percentile.
CalibratorCheck2 <- function(MySoma2022,clinicType,exprDat,titleMessage){
  if(clinicType=="OA"){calibIDs <- paste0("OA POOL-HT-",c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23),"/29")
  }else if(clinicType=="INJ"){calibIDs <- paste0("INJ POOL-HT-",c(1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23),"/25")
  }else if(clinicType=="FreezeThawOA"){calibIDs <- paste0("OA POOL-HT-",c(12,25),"/29")
  }else if(clinicType=="FreezeThawINJ"){calibIDs <- paste0("INJ POOL-HT-",c(13,25),"/25")
  }else{calibIDs <- paste0("UNSP POOL-HT-",c(5,1,3,2,4,6),"/8")}
  
  calib_normClinic3 <- MySoma2022[which(MySoma2022$SampleId %in% calibIDs),colnames(exprDat)]
  temp3 <- apply(calib_normClinic3,2,sd)/apply(calib_normClinic3,2,mean)
  
  ### plot the cumulative distribution of %CV 
  plot(100*quantile(temp3,seq(0,1,length.out=100)),seq(0,1,length.out=100),type="l",xlim=c(0,30),xlab=paste("%CV"),ylab="Cumulative Distribution",main=titleMessage,cex.main=1.5,cex.lab=1.3,font.lab=2)
  lines(c(0,quantile(temp3,0.8))*100,c(0.8,0.8),lty=2,col="red")
  lines(c(quantile(temp3,0.8),quantile(temp3,0.8))*100,c(0.8,0),lty=2,col="red")
  axis(1, font=2)
  axis(2, font=2)
  
  text(100*quantile(temp3,seq(0,1,length.out=100))[80]-1,0.85,formatC(100*quantile(temp3,seq(0,1,length.out=100))[80],digits=4),pos=4,col="red",cex=1.5,font=2)
  return(temp3)
}

#__________________________________________________________________________________________________________
# function to perform biological signal strength evaluation in the RFU (namely biological variance explained in the investigated RFU)
# input: calib_norm -- concentration frame for pooled samples; calibIDs: sample ID for calibrators; clinicType: disease group; 
#        exprDat_norm: protein concentration only frame; MySoma -- RFU frame.
# output: R2_norm -- a vector of R2 based on synovial fluid pooled samples or QC plasma calibrators.
#         plots of accumulative distribution of CVs marked at 80 percentile.
VarExp2 <- function(MySoma,clinicType,exprDat_norm,titleMessage){
  if(clinicType=="OA"){
    calibIDs <- paste0("OA POOL-HT-",c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23),"/29")
  }else if(clinicType=="INJ"){
    calibIDs <- paste0("INJ POOL-HT-",c(1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23),"/25")
  }else if(clinicType=="UNSP"){
    calibIDs <- paste0("UNSP POOL-HT-",c(5,1,3,2,4,6),"/8")
  }else if(clinicType=="FreezeThawOA"){
    calibIDs <- paste0("OA POOL-HT-",c(12,25),"/29")
  }else if(clinicType=="FreezeThawINJ"){
    calibIDs <- paste0("INJ POOL-HT-",c(13,25),"/25")
  }else{
    calibIDs <- clinicType}
  
  calib_normClinic <- MySoma[which(MySoma$SampleId %in% calibIDs),colnames(exprDat_norm)]
  temp <- apply(calib_normClinic,2,var)/apply(exprDat_norm,2,var)
  temp[temp > 1] <- 1
  R2_norm <- (1 - temp)^2
  plot(quantile(R2_norm,1 - seq(0,1,length.out=100)),seq(0,1,length.out=100),xlim=c(1,0),type="l",xlab=paste("R2"),ylab="Cumulative Distribution",main=titleMessage,cex.main=1.5,cex.lab=1.3,font.lab=2)
  axis(1, font=2)
  axis(2, font=2)
  abline(h=0.8,lty=2)
  abline(v=quantile(R2_norm,1 - 0.8),lty=2)
  text(quantile(R2_norm,1 - seq(0,1,length.out=100))[80],0.85,formatC(quantile(R2_norm,1 - seq(0,1,length.out=100))[80],digits=4),pos=4,col="red",cex=1.5,font=2)
  
  return(R2_norm)
}

#__________________________________________________________________________________________________________
# function for external comparison with immunoassay measures
# input: sandwich_master -- measurements for the same set of synovial fluid samples from immunoassay; exprDat_norm -- protein concentration only frame; 
#       Test1, Test2 -- protein names for the immunoassay and somologic platform.  
# output: CorData_norm -- correlation matrix with correlation coefficients, p values and confidence intervals.
ExtVal <- function(sandwich_master,exprDat_norm,Test1,Test2){
  SomaName <- Test2
  colKeep <- unlist(sapply(SomaName,function(x){grep(x,colnames(exprDat_norm))}))
  rowKeep <- unlist(sapply(sandwich_master$PIN,function(x){which(grepl(x,rownames(exprDat_norm)))}))
  
  SomaValue <- exprDat_norm[rowKeep,colKeep]
  
  sandwichKeepRow <- sub("-.*$","",rownames(SomaValue))
  
  Test1 <- Test1[unlist(sapply(colnames(SomaValue),function(x){which(Test2==x)}))]
  
  #check Ben (this is all OA) / Historic data (this is all knee injury)
  CorData_norm = data.frame(matrix(NA,nrow=length(Test1),ncol=7))
  names(CorData_norm) <- c("SandwichName","SomaName","cor","Pvalue","conf.int.lower","conf.int.upper","N")
  
  for(testCounter in 1:length(Test1)){
    if (sum(!(is.na(sandwich_master[sandwichKeepRow,Test1[testCounter]] + SomaValue[,testCounter]))) == 0){CorData_norm[testCounter,]=(c(Test1[testCounter],Test2[testCounter],NA,NA,NA,NA,0))
    }else{md <- cor.test(sandwich_master[sandwichKeepRow,Test1[testCounter]],SomaValue[,testCounter])
    CorData_norm[testCounter,] = c(Test1[testCounter],colnames(SomaValue)[testCounter],unlist(md[c("estimate","p.value","conf.int")]),2+md$parameter)}
  }
  
  rownames(CorData_norm) <- CorData_norm$SandwichName
  return(CorData_norm)
}

#__________________________________________________________________________________________________________
# function to obtain the combined data frame between RFU and clinical meta information.
# input: pathIn -- file path of clinicl meta files; MySoma -- RFU frames; ProMeta -- protein meta table derived from somalogic adat file; 
# output: list contains combined data frame (CombinedFrame), RFU frame with observations only have matching clinical meta information (MySoma); protein concentration only frame (exprDat_norm).
# note: clinical meta files "discovery_QApheno_1_120522.csv", "Additional data for Yun_TP_27-05-22.xlsx", contain variables such as freeze thaw, processing batch, age ect. 
#       "adat_to_redcap_SIN_map.csv" is a self generated file which map between adat STEpUp ID and redcap STEpUp ID
getFrarmeList <- function(pathIn,MySoma,ProMeta){
  ### read in clinic metadata
  ClinicFrame <- read.csv(paste0(pathIn,'discovery_QApheno_1_120522.csv'))
  ClinicFrameTempSup <- read_excel(paste0(pathIn,'Additional data for Yun_TP_27-05-22.xlsx'))
  ClinicFrame <- merge(ClinicFrame,ClinicFrameTempSup[,c("sf_iknee_sample_id_number","sf_iknee_proc_treat_date","age")],by="sf_iknee_sample_id_number")
  
  ###SIN mapper
  SINmap <- read.csv(paste0(pathIn,'adat_to_redcap_SIN_map.csv'),sep=" ")
  
  ### change  SIN in adat frame MySoma to match ClinicFrame
  MySoma$SampleId <- unlist(sapply(MySoma$SampleId,function(x){
    if(x%in%SINmap$adat_SINs){SINmap$RedCap_SINs[which(SINmap$adat_SINs==x)]
    }else{x}
  }))
  
  ### combine MySoma with ClinicFrame, based on STEpUPId
  colnames(MySoma)[which(colnames(MySoma)=="SampleId")] <- "STEpUpId" ### After Combat, returned RFU frame only includes human proteins.
  colnames(ClinicFrame)[which(colnames(ClinicFrame)=="sf_iknee_sample_id_number")] <- "STEpUpId"
  CombinedFrame <- merge(ClinicFrame,MySoma,by="STEpUpId")
  rownames(CombinedFrame) <- CombinedFrame$STEpUpId
  
  ### set class for clinical features: for regression.
  CombinedFrame$sex <- as.factor(CombinedFrame$sex)
  CombinedFrame$sf_iknee_bloodstaining <- sub(5,NA,CombinedFrame$sf_iknee_bloodstaining) ### 5 in blood staining means missing data
  CombinedFrame$cohort_name <- as.factor(CombinedFrame$cohort_name)
  CombinedFrame$sf_iknee_proc_batch <- as.factor(sub(-999,NA,CombinedFrame$sf_iknee_proc_batch)) ### -999 in processing batch means missing data
  CombinedFrame$PlateId <- as.factor(CombinedFrame$PlateId)
  CombinedFrame$sf_iknee_qc_group <- as.factor(CombinedFrame$sf_iknee_qc_group)
  CombinedFrame$PlatePosition <- as.factor(CombinedFrame$PlatePosition)
  CombinedFrame$PlateRunDate <- as.factor(CombinedFrame$PlateRunDate)
  CombinedFrame$sl_tranche_number <- as.factor(CombinedFrame$sl_tranche_number)
  CombinedFrame$sf_iknee_prev_freeze_thaw <- as.factor(sub(2,NA,CombinedFrame$sf_iknee_prev_freeze_thaw)) ### 2 in freeze thaw cycle means missing data
  CombinedFrame$sf_iknee_proc_treat_date <- as.factor(CombinedFrame$sf_iknee_proc_treat_date)
  CombinedFrame$sample_age <- log(CombinedFrame$sample_age)
  CombinedFrame$sf_iknee_volume <- log(CombinedFrame$sf_iknee_volume)
  
  ### extract the data only zone: exprDat_norm. only human proteins are involved in our downstream analysis.
  exprDat_norm.HM <- CombinedFrame[,which(colnames(CombinedFrame)=="seq.10000.28")-1+which(ProMeta$Type=="Protein" & ProMeta$Organism=="Human")]
  exprDat_norm <- CombinedFrame[,which(colnames(CombinedFrame)=="seq.10000.28"):ncol(CombinedFrame)]
  rownames(exprDat_norm) <- CombinedFrame$STEpUpId
  
  return(list(MySoma,CombinedFrame,exprDat_norm,exprDat_norm.HM))
}

#__________________________________________________________________________________________________________
### function to obtain top PCs which explain the "thresh" percentage variations
### Input: exprDat -- RFU; thresh -- threshold to define top PCs
### Output: pcDat -- matrix of the top PCs; varianceExp -- variation explained per PC
getTopPCdt <- function(exprDat,thresh){
  pc_norm <- prcomp(log10(as.matrix(exprDat)),scale = TRUE)
  eig.val <- get_eigenvalue(pc_norm) 
  topPC <- which(eig.val$cumulative.variance.percent>thresh)[1]
  pcDat <- pc_norm$x[,1:topPC]
  colnames(pcDat) <- paste0("PC",1:topPC)
  varianceExp <- eig.val$variance.percent
  return(list(pcDat,varianceExp))
}

#__________________________________________________________________________________________________________
# function to generate pvalue matrix/effect size matrix, showing association between technical confounders and PCs or proteins.
# Input: CombinedFrame -- combined frame between RFU and clinical meta information; whichDat -- data frame for confounder check, per protein level, PC level or total protein level.
# output: list of a matrix of p values and a matrix of effect size, representing associations with confounders.
ConfouderCheck <- function(CombinedFrame,whichDat){
  
  confList1 = c("PlateId","sf_iknee_qc_group","PlatePosition","PlateRunDate")
  confList2 = c("sl_tranche_number","sf_iknee_prev_freeze_thaw","sf_iknee_freezethaw_cycles","sf_iknee_proc_batch","sf_iknee_proc_treat_date","sample_age","sf_iknee_volume","sf_iknee_bloodstaining")
  confList = c(confList1,confList2)
  
  ### p container based on each protein, each PC, total protein
  ### varContainer correspondings to pContainer, used to accormodate effect size for each test.
  pContainer = matrix(NA,nrow=ncol(whichDat),ncol=length(confList))
  rownames(pContainer) = colnames(whichDat)
  colnames(pContainer) = confList
  varContainer <- pContainer
  
  for (confounder in confList){
    
    ### consider each protein 
    for(ProCt in 1:ncol(whichDat)){
      print(ProCt)
      proName = colnames(whichDat)[ProCt]
      
      if(confounder %in% confList1){ # some predictors are conditioning on cohort
        f <- as.formula(paste0(proName,"~",confounder))
      }else{f <- as.formula(paste0(proName,"~","cohort_name+",confounder))}
      
      tempdata=data.frame(cbind(whichDat[,proName],CombinedFrame[,c(confounder,"cohort_name")]))
      tempdata_noNA <- tempdata[!apply(tempdata,1,function(x) any(is.na(x))),]
      colnames(tempdata_noNA)[1] <- proName
      model <- lm(formula=f,data=tempdata_noNA)
      tempObj <- anova(model)
      
      pContainer[proName,confounder] <- tempObj[confounder,"Pr(>F)"]
      varContainer[proName,confounder] <- tempObj[confounder,"Sum Sq"]/sum(tempObj[confounder,"Sum Sq"]+tempObj["Residuals","Sum Sq"])
    }
    print(paste0(Sys.time(),", confounder ",confounder, " rergression done.")) # monitor the progress
  }
  
  return(list(pContainer,varContainer))  
}

#__________________________________________________________________________________________________________
# function to perform total protein check. 
# Input: CombinedFrame -- combined data frame between RFU and clinical meta information; exprDat_norm -- the protein concentration frame; titleMessage -- for the convenience to generate plot titles.
# output: total protein distribution as well as broken down by disease groups; correlation plot between blood staining and total protein concentration.
#         totalProtein_norm -- total protein concentration for the investigated RFU frame.
TotalProCheck <- function(CombinedFrame,exprDat_norm,titleMessage){
  bloodStain <- CombinedFrame$sf_iknee_bloodstaining
  diseaseGroup <- unlist(sapply(as.numeric(CombinedFrame$sf_iknee_qc_group),function(x){base::switch(x,"OA","injury","control","RA")}))
  keepID <- which(!is.na(CombinedFrame$sf_iknee_qc_group))
  
  totalProtein_norm <- apply(exprDat_norm,1,sum)
  hist(totalProtein_norm,breaks=100,xlab="total protein amount for each sample",ylab="frequency",main="Histogram of Total Protein",cex.main=1.6,cex.lab=1.5)
  plot(as.factor(bloodStain),totalProtein_norm,xlab="blood staining grade",ylab="total protein",main=paste0("Correlation between Blood Stain and Total Protein\n",titleMessage),cex.main=1.6,cex.lab=1.5)
  abline(lm(totalProtein_norm ~ as.numeric(bloodStain)),col="red",lwd=2)
  text(1,max(totalProtein_norm)*0.5,paste("p =",formatC(cor.test(totalProtein_norm, as.numeric(bloodStain))$p.val,format="e",digits=3)),pos=4,col="red",cex=1.5)
  text(1,max(totalProtein_norm)*0.55,paste("cor =",formatC(cor.test(totalProtein_norm, as.numeric(bloodStain))$estimate,format="e",digits=3)),pos=4,col="red",cex=1.5)
  boxplot(totalProtein_norm[keepID] ~ diseaseGroup[keepID],main=paste("p =",signif(kruskal.test(totalProtein_norm[keepID] ~ diseaseGroup[keepID])$p.val,digits=4)),xlab="disease group",ylab="total protein within group",cex.main=2,cex.lab=1.5)
  return(totalProtein_norm)
}

#__________________________________________________________________________________________________________
# function to perform limit of detection check. Lower limit of detection (LoD) and upper limit of detection (uLoD).
# input: RawM -- investigated RFU frame.
# output: CompM -- matrix with rows as samples and columns as proteins, each entry is the difference between the protein measures and limit of detection threshholds.
#         histogram of proportions of unqualified proteins and samples. 
LoDdetection <- function(RawM){
  ###for calculation convenience, extract data zone only
  
  ### limit of detection calculation LoD = median(Buffer) + 4.9*MAD(Buffer) (median absolute deviation)
  MBuffer = RawM[which(RawM$SampleType == "Buffer"),]
  datBuffer = MBuffer[,which(colnames(RawM)=="seq.10000.28"):ncol(MBuffer)]
  BuffMedian = apply(datBuffer,2,median)
  LoD = BuffMedian + 4.9* apply(t(apply(datBuffer,1,function(x){abs(x - BuffMedian)})),2,median)
  names(LoD) = "LoD"
  
  MSample = RawM[grep("Sample",RawM$SampleType),]
  DatSamp = MSample[,which(colnames(MSample)=="seq.10000.28"):ncol(MSample)]
  
  ### which RFU is below LoD? CompM
  CompM = t(apply(DatSamp,1,function(x){x-LoD}))
  
  ### histgram of proportion below LoD per protein
  ProteinRatio = apply(CompM,2,function(x){length(which(x<0))/nrow(CompM)})
  ### histgram of proportion below LoD per sammple
  SampRatio = apply(CompM,1,function(x){length(which(x<0))/ncol(CompM)})
  
  hist(ProteinRatio,main="Protein Expression Level \nBelow the Lower Limit of Detection Across Proteins",xlab="ratio among all the proteins",breaks=50,cex.main=1)
  hist(SampRatio,main="Protein Expression Level \nBelow the Lower Limit of Detection Across Samples",xlab="ratio among all the samples",breaks=50,cex.main=1)
  
  return(CompM)
}

uLoDdetection <- function(RawM){
  
  MSample = RawM[grep("Sample",RawM$SampleType),]
  DatSamp = MSample[,which(colnames(MSample)=="seq.10000.28"):ncol(MSample)]
  
  ### which RFU is above 80,000
  CompM = as.matrix(DatSamp - 80000)
  
  ### histgram of proportion above 80,000 per protein
  ProteinRatio = apply(CompM,2,function(x){length(which(x>0))/nrow(CompM)})
  ### histgram of proportion below LoD per sammple
  SampRatio = apply(CompM,1,function(x){length(which(x>0))/ncol(CompM)})
  
  hist(ProteinRatio,main="Protein Expression Level \nAbove the Upper Limit of Detection Accross Proteins",xlab="ratio among all the proteins",breaks=50,cex.main=1)
  hist(SampRatio,main="Protein Expression Level \nAbove the Upper Limit of Detection Accross Samples",xlab="ratio among all the samples",breaks=50,cex.main=1)
  
  return(CompM)
}

#__________________________________________________________________________________________________________
# function to calculate associations between proteins and bimodal signal status
# input: exprDat -- protein expression profile. bimodalLabel -- bimodal signal status.
# output: bimodalP -- pvalue showing associations between protein and bimodal signal status.
bimodalTest <- function(exprDat,bimodalLabel){
  bimodalP=vector(mode="numeric",length=ncol(exprDat))
  names(bimodalP) = colnames(exprDat)
  for(ProCt in 1:ncol(exprDat)){
    print(ProCt) # monitor the progress
    proName <- colnames(exprDat)[ProCt]
    model0 <- lm(exprDat[,ProCt]~1)
    model1 <- lm(exprDat[,ProCt]~as.factor(bimodalLabel))
    tempObj <- anova(model0,model1)
    bimodalP[proName] <- tempObj$`Pr(>F)`[2]
  }
  return(bimodalP)
}

#__________________________________________________________________________________________________________
# function to obtain summary of filters for proteins and samples
# input: previously calculated removed items; ConfounderTable1: p values showing associations between proteins and confounders; exprDat: RFU
# output: removeTable -- summary table for filtering; exprDatFiltered --  RFU contains the remaining proteins and samples.
getRemoveTable <- function(NonHuman,RemoveS1,RemoveS2,RemoveS3,RemoveS4,RemovePro1,RemovePro2,bimodalP,ConfounderTable1,exprDat){
  RemovePro3 = names(which(bimodalP<0.05/7289)) ### bimodal signal test
  RemoveConf <- sapply(colnames(ConfounderTable1),function(x){rownames(ConfounderTable1)[which(ConfounderTable1[,x] < 0.05/7289)]})
  removeprotein <- unique(c(unlist(RemoveConf),NonHuman,RemovePro1,RemovePro2,RemovePro3)) 
  removesample = unique(c(RemoveS1,RemoveS2,RemoveS3,RemoveS4))
  exprDatFiltered <- exprDat[which(!(rownames(exprDat) %in% removesample)),which(!(colnames(exprDat) %in% removeprotein))]
  
  removeTable1 <- rbind(t(data.frame("total remove protein"=length(removeprotein),"total remove sample"=length(removesample),
                                     "remain protein"=ncol(exprDatFiltered),"remain sample"=nrow(exprDatFiltered),
                                     "NonHuman" = length(NonHuman),
                                     "OA repeats"=length(RemovePro1),"Injury repeats"=length(RemovePro2),"bimodal signal"=length(RemovePro3),
                                     "PCAout"=length(RemoveS1),"TotalProteinOut"=length(RemoveS2),"LoD|uLoD(S)"=length(RemoveS3),"SomaLogicInHouse"=length(RemoveS4))))
  colnames(removeTable1) <- "moved amount"
  
  removeTable2 <- data.frame(sapply(names(RemoveConf),function(x){x=length(RemoveConf[[x]])}))
  colnames(removeTable2) <- "moved amount"
  
  removeTable <- rbind(removeTable1,removeTable2)
  
  return(list(removeTable,exprDatFiltered))
}

#__________________________________________________________________________________________________________
# function to find proteins for a particular subcellular location
# input: enrichSource -- online resource data for proteins and their subcellular locations; enrichMessage -- annotation of which sublocations are inverstigated
# output: a list;  enrichSourceDat -- dataframe of gene name and matching main location; enrichSourceType-- character vector showing all the locations; GeneSets -- list, location and corresponding gene names)
diffProCombineOnline <- function(enrichSource,enrichMessage){
  
  enrichSourceDat <- enrichSource[,c("Gene.name",enrichMessage)]
  enrichSourceTypeR <- levels(as.factor(enrichSource[,enrichMessage]))
  
  ### extract all the types for enrichment test. delimiter: "sepSym" should be selected by looking at the raw data frame
  if(length(grep("location",enrichMessage))!=0){sepSym = ";"
  }else if(length(grep("Tissue",enrichMessage))!=0){sepSym=", "
  }else if(length(grep("Cell.type",enrichMessage))!=0){sepSym="/"} 
  else{sepSym = ""} ### leave the chance to add more enrichMessage types  
  
  ### process element with multiple entries
  if(sepSym != ""){
    k=1
    enrichSourceType=vector()
    for (locationCouter in 1:length(enrichSourceTypeR)){
      templocal = enrichSourceTypeR[locationCouter]
      if(!grepl(sepSym,templocal)){enrichSourceType[k]=templocal
      k=k+1}
      else{templocal2 <- strsplit(templocal,sepSym)[[1]]
      for(splitC in 1:length(templocal2)){enrichSourceType[k+splitC-1]=templocal2[splitC]}
      k=k+length(templocal2)
      }
    }
  }
  
  enrichSourceTypeA = names(table(enrichSourceType))
  
  ###construct reference genesets, with elements of type name and matching EntrezGene name. 
  GeneSets=vector(mode="list",length=length(enrichSourceTypeA))
  for (subCounter in 1:length(enrichSourceTypeA)){
    GeneSets[[subCounter]] = enrichSourceDat[grep(enrichSourceTypeA[subCounter],enrichSourceDat[,enrichMessage]),"Gene.name"]
    names(GeneSets)[subCounter] = enrichSourceTypeA[subCounter]
  } 
  
  return(list(enrichSourceDat,enrichSourceType,GeneSets))
}

#__________________________________________________________________________________________________________
# function to calculate correlation coefficients between each protein and PC1
# input: exprHere -- RFU; ProMeta -- protein meta information; pcDat -- top PCs (80% variations explained)
# output: a list;  corPerPro -- Pearson correlation coefficient between protein and PC1, mean.abundacne -- mean abundance of proteins, adjusted for dilution bins)
getCorPro <- function(exprHere,ProMeta,pcDat){
  abundance <- sapply(1:ncol(exprHere),function(x){exprHere[,x]/ProMeta[colnames(exprHere)[x],"Dilution2"]})
  colnames(abundance) <- colnames(exprHere)
  mean.abundacne <- apply(abundance,2,mean)
  keepseq <- rownames(ProMeta)[which(ProMeta$Organism=="Human" & ProMeta$Type=="Protein")]
  corPerPro <- sapply(1:ncol(exprHere), function(x){cor(pcDat[,1],log(exprHere[,x]))})
  names(corPerPro) <- colnames(exprHere)
  return(list(corPerPro,mean.abundacne))
}

#__________________________________________________________________________________________________________
# function to obtain correlation  coefficients between RFU and external immunoassay measures among different standardization steps
# input: externalCheck -- correlation test results
# output: CorDatP -- data frame showing correlation between RFU and immunoassay measures
plotExternal9 <- function(externalCheck){
  CorDatY = as.matrix(externalCheck[,3],ncol=1)
  CorDatX = rep(as.matrix(externalCheck[,1],ncol=1),(ncol(externalCheck)-2))
  CorC = rep(seq(1:(ncol(externalCheck)-2)),rep(nrow(externalCheck),ncol(externalCheck)-2))
  
  for(k in 4:ncol(externalCheck)){
    CorDatY =  rbind(CorDatY,as.matrix(externalCheck[,k],ncol=1))
  }
  colnames(CorDatY) = "CorDatY"
  CorDatP = as.data.frame(cbind(CorDatX,CorC,CorDatY))
  return(CorDatP)
}

#__________________________________________________________________________________________________________
# function to easily write R table object into excel.
# input: sheetname -- excel file name; whichX -- which R object to be put out.
# output: excel file contain aimed R object.
outEXL <- function(pathOut,sheetname,whichX){
  wb <- createWorkbook()
  addWorksheet(wb, sheet = sheetname, gridLines = TRUE)
  writeData(wb, sheet = sheetname, x=whichX, rowNames = TRUE, withFilter = FALSE)
  saveWorkbook(wb, paste0(pathOut,sheetname,".xlsx"), overwrite = TRUE)
  return()
}

#__________________________________________________________________________________________________________
# function to visualize technical confounders on their most associated two PCs
# input: CombinedFrame -- combined frame between RFU and clinical data meta information; pcDatHere -- PCA space (80% variation explained); 
# output: cowplot
plotTechPC <- function(CombinedFrame,pcDatHere){
  ConfouderCheck2.QCed <- ConfouderCheck(CombinedFrame[rownames(pcDatHere),],pcDatHere)
  ConfounderTable2.QCed <- ConfouderCheck2.QCed[[1]]
  
  ### extract and record the top2 PCs
  which2PC.QCed = matrix(NA,ncol=2,nrow=ncol(ConfounderTable2.QCed))
  colnames(which2PC.QCed) <- c("topPC1","topPC2")
  keepID.QCed = vector(mode="list",length=ncol(ConfounderTable2.QCed))
  rownames(which2PC.QCed) <- names(keepID.QCed) <- colnames(ConfounderTable2.QCed)
  
  CombinedFrameKeep <- CombinedFrame[rownames(pcDatHere),]
  confounders <- colnames(ConfounderTable2.QCed)
  for(confounder in confounders){
    keepID.QCed[[confounder]] <- which(!is.na(CombinedFrameKeep[,confounder]))
    which2PC.QCed[confounder,] <- rownames(ConfounderTable2.QCed)[order(ConfounderTable2.QCed[1:10,confounder])[c(1,2)]]
  }
  
  plate.position.means <- aggregate(pcDatHere[keepID.QCed[["PlatePosition"]],which2PC.QCed["PlatePosition","topPC1"]],list(CombinedFrameKeep[keepID.QCed[["PlatePosition"]],"PlatePosition"]),FUN=mean)
  
  #calculate the mean value of the most significant PC by well
  plate.position.means$row <- match(substr(plate.position.means$Group.1,1,1),LETTERS)
  plate.position.means$col <- as.numeric(substr(plate.position.means$Group.1,2,10))
  
  #make the plate plot (it is on its side)
  p.platePosition.plan.QC <- ggplot(data = plate.position.means) + 
    geom_circle(aes(x0 = col, y0 = row, r = 0.45, fill = x)) + 
    coord_equal() + 
    scale_x_continuous(breaks = 1:12, expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(breaks = 1:8, labels = LETTERS[1:8], expand = expansion(mult = c(0.01, 0.01)), trans = reverse_trans()) + 
    scale_fill_distiller(type="div",palette ="RdYlGn") +
    theme_bw() + xlab("Column") + ylab("Row") + labs(fill = paste0("Mean of \n",which2PC.QCed["PlatePosition",1]," for well")) + 
    theme(axis.title.x = element_text(size=11,face="bold"),legend.title =element_text(size = 8,face="bold"), legend.text = element_text(size = 8,face="bold"), 
          axis.title.y =element_text(size=11,face="bold"), legend.position = "bottom")
  
  
  p.plate.QCed <- ggplot() + geom_point(aes(x=pcDatHere[keepID.QCed[["PlateId"]],which2PC.QCed["PlateId","topPC1"]],y=pcDatHere[keepID.QCed[["PlateId"]],which2PC.QCed["PlateId","topPC2"]],color=as.factor(CombinedFrameKeep[keepID.QCed[["PlateId"]],"PlateId"])),size=0.8) +
    labs(color="Plate Number") + xlab(which2PC.QCed["PlateId",1]) + ylab(which2PC.QCed["PlateId",2]) + scale_color_discrete(labels=seq(1:22)) + theme_bw() + 
    theme(axis.title.x = element_text(size=11,face="bold"),legend.title =element_text(size = 8,face="bold"), legend.text = element_text(size = 8,face="bold"), 
          axis.title.y =element_text(size=11,face="bold"))
  
  p.sampleAge.QCed <- ggplot() + geom_point(aes(x=pcDatHere[keepID.QCed[["sample_age"]],which2PC.QCed["sample_age","topPC1"]],y=pcDatHere[keepID.QCed[["sample_age"]],which2PC.QCed["sample_age","topPC2"]],color=as.numeric(exp(CombinedFrameKeep[keepID.QCed[["sample_age"]],"sample_age"]))),size=0.8) +
    scale_colour_distiller(type="div",palette ="RdYlGn") + labs(color="Sample\nAge") + xlab(which2PC.QCed["sample_age",1]) + ylab(which2PC.QCed["sample_age",2]) + theme_bw() + 
    theme(axis.title.x = element_text(size=11,face="bold"),legend.title =element_text(size = 8,face="bold"), legend.text = element_text(size = 8,face="bold"), 
          axis.title.y =element_text(size=11,face="bold"),legend.position = "bottom")
  
  p.FreezeThawCycle.QCed <- ggplot() + geom_point(aes(x=pcDatHere[keepID.QCed[["sf_iknee_freezethaw_cycles"]],which2PC.QCed["sf_iknee_freezethaw_cycles","topPC1"]],y=pcDatHere[keepID.QCed[["sf_iknee_freezethaw_cycles"]],which2PC.QCed["sf_iknee_freezethaw_cycles","topPC2"]],color=as.numeric(CombinedFrameKeep[keepID.QCed[["sf_iknee_freezethaw_cycles"]],"sf_iknee_freezethaw_cycles"])),size=0.8) +
    scale_colour_continuous(low = "lightblue",high="red")  + labs(color="Freeze-Thaw\nCycles") + xlab(which2PC.QCed["sf_iknee_freezethaw_cycles",1]) + ylab(which2PC.QCed["sf_iknee_freezethaw_cycles",2]) + theme_bw() + 
    theme(axis.title.x = element_text(size=11,face="bold"),legend.title =element_text(size = 8,face="bold"), legend.text = element_text(size = 8,face="bold"), 
          axis.title.y =element_text(size=11,face="bold"),legend.position = "bottom")
  
  DiseaseSize=ifelse(CombinedFrameKeep[keepID.QCed[["sf_iknee_qc_group"]],"sf_iknee_qc_group"]==3,"Big",ifelse(CombinedFrameKeep[keepID.QCed[["sf_iknee_qc_group"]],"sf_iknee_qc_group"]==2,"Mid","Small"))
  p.disease.QCed <- ggplot() + geom_point(aes(x=pcDatHere[keepID.QCed[["sf_iknee_qc_group"]],which2PC.QCed["sf_iknee_qc_group","topPC1"]],y=pcDatHere[keepID.QCed[["sf_iknee_qc_group"]],which2PC.QCed["sf_iknee_qc_group","topPC2"]],color=as.factor(CombinedFrameKeep[keepID.QCed[["sf_iknee_qc_group"]],"sf_iknee_qc_group"]),size=DiseaseSize)) +
    labs(color="Disease\nGroup") + xlab(which2PC.QCed["sf_iknee_qc_group",1]) + ylab(which2PC.QCed["sf_iknee_qc_group",2]) + scale_colour_manual(values = c("orange","green","black","purple"),labels=c('OA', 'Injury', 'Healthy\ncontrol',"Inflammatory\ncontrol")) + 
    scale_size_manual (values= c(1.5,1,0.03)) + guides(size=FALSE,color = guide_legend(override.aes = list(size = 10))) + theme_bw() + 
    theme(axis.title.x = element_text(size=11,face="bold"),legend.title =element_text(size = 8,face="bold"), legend.text = element_text(size = 8,face="bold"), 
          axis.title.y =element_text(size=11,face="bold"),legend.position = "bottom")  + guides(color=guide_legend(nrow=2, byrow=TRUE)) 
  
  p.sampleVol.QCed <- ggplot() + geom_point(aes(x=pcDatHere[keepID.QCed[["sf_iknee_volume"]],which2PC.QCed["sf_iknee_volume","topPC1"]],y=pcDatHere[keepID.QCed[["sf_iknee_volume"]],which2PC.QCed["sf_iknee_volume","topPC2"]],color=as.numeric(exp(CombinedFrameKeep[keepID.QCed[["sf_iknee_volume"]],"sf_iknee_volume"]))),size=0.8) +
    scale_colour_distiller(type="div",palette ="RdYlGn",limits=c(0,30),oob = squish)  + labs(color="Sample\nVolume") + xlab(which2PC.QCed["sf_iknee_volume",1]) + ylab(which2PC.QCed["sf_iknee_volume",2]) + theme_bw() + 
    theme(axis.title.x = element_text(size=11,face="bold"),legend.title =element_text(size = 8,face="bold"), legend.text = element_text(size = 8,face="bold"), 
          axis.title.y =element_text(size=11,face="bold"),legend.position = "bottom") 
  
  p.blood.QCed <- ggplot() + geom_point(aes(x=pcDatHere[keepID.QCed[["sf_iknee_bloodstaining"]],which2PC.QCed["sf_iknee_bloodstaining","topPC1"]],y=pcDatHere[keepID.QCed[["sf_iknee_bloodstaining"]],which2PC.QCed["sf_iknee_bloodstaining","topPC2"]],color=as.numeric(CombinedFrameKeep[keepID.QCed[["sf_iknee_bloodstaining"]],"sf_iknee_bloodstaining"])),size=0.8) +
    scale_colour_distiller(type="div",palette ="RdYlGn")  + labs(color="Blood\nStaining") + xlab(which2PC.QCed["sf_iknee_bloodstaining",1]) + ylab(which2PC.QCed["sf_iknee_bloodstaining",2]) + theme_bw() + 
    theme(axis.title.x = element_text(size=11,face="bold"),legend.title =element_text(size = 8,face="bold"), legend.text = element_text(size = 8,face="bold"), 
          axis.title.y =element_text(size=11,face="bold"),legend.position = "bottom")
  
  p.platePosition.QCed <- ggplot() + geom_point(aes(x=pcDatHere[keepID.QCed[["PlatePosition"]],which2PC.QCed["PlatePosition","topPC1"]],y=pcDatHere[keepID.QCed[["PlatePosition"]],which2PC.QCed["PlatePosition","topPC2"]],color=as.factor(CombinedFrameKeep[keepID.QCed[["PlatePosition"]],"PlatePosition"])),size=0.8) +
    labs(color="Plate Position") + xlab(which2PC.QCed["PlatePosition",1]) + ylab(which2PC.QCed["PlatePosition",2]) + theme_bw() + 
    theme(axis.title.x = element_text(size=11,face="bold"),legend.title =element_text(size = 8,face="bold"), legend.text = element_text(size = 8,face="bold"), 
          axis.title.y =element_text(size=11,face="bold"))
  
  p.plateRunDate.QCed <- ggplot() + geom_point(aes(x=pcDatHere[keepID.QCed[["PlateRunDate"]],which2PC.QCed["PlateRunDate","topPC1"]],y=pcDatHere[keepID.QCed[["PlateRunDate"]],which2PC.QCed["PlateRunDate","topPC2"]],color=as.factor(CombinedFrameKeep[keepID.QCed[["PlateRunDate"]],"PlateRunDate"])),size=0.8) +
    labs(color="Plate Run Date") + xlab(which2PC.QCed["PlateRunDate",1]) + ylab(which2PC.QCed["PlateRunDate",2]) + theme_bw() + 
    theme(axis.title.x = element_text(size=11,face="bold"),legend.title =element_text(size = 8,face="bold"), legend.text = element_text(size = 8,face="bold"), 
          axis.title.y =element_text(size=11,face="bold"))
  
  p.tranche.QCed <- ggplot() + geom_point(aes(x=pcDatHere[keepID.QCed[["sl_tranche_number"]],which2PC.QCed["sl_tranche_number","topPC1"]],y=pcDatHere[keepID.QCed[["sl_tranche_number"]],which2PC.QCed["sl_tranche_number","topPC2"]],color=as.factor(CombinedFrameKeep[keepID.QCed[["sl_tranche_number"]],"sl_tranche_number"])),size=0.8) +
    labs(color="Tranche Number") + xlab(which2PC.QCed["sl_tranche_number",1]) + ylab(which2PC.QCed["sl_tranche_number",2]) + theme_bw() + 
    theme(axis.title.x = element_text(size=11,face="bold"),legend.title =element_text(size = 8,face="bold"), legend.text = element_text(size = 8,face="bold"), 
          axis.title.y =element_text(size=11,face="bold"),legend.position = "bottom")
  
  p.previous.FreezeThaw.QCed <- ggplot() + geom_point(aes(x=pcDatHere[keepID.QCed[["sf_iknee_prev_freeze_thaw"]],which2PC.QCed["sf_iknee_prev_freeze_thaw","topPC1"]],y=pcDatHere[keepID.QCed[["sf_iknee_prev_freeze_thaw"]],which2PC.QCed["sf_iknee_prev_freeze_thaw","topPC2"]],color=as.factor(CombinedFrameKeep[keepID.QCed[["sf_iknee_prev_freeze_thaw"]],"sf_iknee_prev_freeze_thaw"])),size=0.8) +
    labs(color="Previous Freeze Thaw") + xlab(which2PC.QCed["sf_iknee_prev_freeze_thaw",1]) + ylab(which2PC.QCed["sf_iknee_prev_freeze_thaw",2]) + theme_bw() + 
    theme(axis.title.x = element_text(size=11,face="bold"),legend.title =element_text(size = 8,face="bold"), legend.text = element_text(size = 8,face="bold"), 
          axis.title.y =element_text(size=11,face="bold"))
  
  p.procBatch.QCed <- ggplot() + geom_point(aes(x=pcDatHere[keepID.QCed[["sf_iknee_proc_batch"]],which2PC.QCed["sf_iknee_proc_batch","topPC1"]],y=pcDatHere[keepID.QCed[["sf_iknee_proc_batch"]],which2PC.QCed["sf_iknee_proc_batch","topPC2"]],color=as.factor(CombinedFrameKeep[keepID.QCed[["sf_iknee_proc_batch"]],"sf_iknee_proc_batch"])),size=0.8) +
    labs(color="Processing Batch") + xlab(which2PC.QCed["sf_iknee_proc_batch",1]) + ylab(which2PC.QCed["sf_iknee_proc_batch",2]) + theme_bw() + 
    theme(axis.title.x = element_text(size=11,face="bold"),legend.title =element_text(size = 8,face="bold"), legend.text = element_text(size = 8,face="bold"), 
          axis.title.y =element_text(size=11,face="bold"))
  
  p.procTreatDate.QCed <- ggplot() + geom_point(aes(x=pcDatHere[keepID.QCed[["sf_iknee_proc_treat_date"]],which2PC.QCed["sf_iknee_proc_treat_date","topPC1"]],y=pcDatHere[keepID.QCed[["sf_iknee_proc_treat_date"]],which2PC.QCed["sf_iknee_proc_treat_date","topPC2"]],color=as.factor(CombinedFrameKeep[keepID.QCed[["sf_iknee_proc_treat_date"]],"sf_iknee_proc_treat_date"])),size=0.8) +
    labs(color="Process Treat Date") + xlab(which2PC.QCed["sf_iknee_proc_treat_date",1]) + ylab(which2PC.QCed["sf_iknee_proc_treat_date",2]) + theme_bw() + 
    theme(axis.title.x = element_text(size=11,face="bold"),legend.title =element_text(size = 8,face="bold"), legend.text = element_text(size = 8,face="bold"), 
          axis.title.y =element_text(size=11,face="bold"))
  
  cowplot::plot_grid(p.platePosition.plan.QC, p.blood.QCed, p.sampleVol.QCed, p.sampleAge.QCed, p.FreezeThawCycle.QCed,p.disease.QCed,
                     labels = c("A", "B","C","D", "E","F"),nrow=2,ncol=3,align = "v")
  
}

#__________________________________________________________________________________________________________
# functions to breakdown CVs for disease group, spun/unspun, using synovial fluid pooled samples.
# input: calib_norm -- concentration frame for pooled samples and QC plasma calibrators; calibIDs: sample ID for calibrators; clinicType: disease group
# output: list of %CV distributions 
CVbreak <- function(calib_norm,calibIDs,clinicType){
  
  if(clinicType=="OA"){
    suffix = "/29"
    acrossNUM = c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
    acrossPlates <- calibIDs %in% paste0(clinicType," POOL-HT-",acrossNUM,suffix)
    temp1 <- apply(calib_norm[acrossPlates,],2,sd)/apply(calib_norm[acrossPlates,],2,mean)
    
    freezeThawNUM = c(12,25) ###12-25;
    freezeThaw <-  calibIDs %in% paste0(clinicType," POOL-HT-",freezeThawNUM,suffix)
    temp3 <- apply(calib_norm[freezeThaw,],2,sd)/apply(calib_norm[freezeThaw,],2,mean)
    
    freshNUM = c(1,3,5,26,27,28) ###1-26; 3-27;5-28
    freshHT <-  calibIDs %in% paste0(clinicType," POOL-HT-",freshNUM,suffix)
    temp5 <- apply(calib_norm[freshHT,],2,sd)/apply(calib_norm[freshHT,],2,mean)
    
    tempList <- list(temp1,temp3,temp5)
    
  }else if(clinicType=="INJ"){
    suffix = "/25"
    acrossNUM = c(1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
    acrossPlates <- calibIDs %in% paste0(clinicType," POOL-HT-",acrossNUM,suffix)
    temp1 <- apply(calib_norm[acrossPlates,],2,sd)/apply(calib_norm[acrossPlates,],2,mean)
    
    freezeThawNUM = c(13,25)
    freezeThaw <-  calibIDs %in% paste0(clinicType," POOL-HT-",freezeThawNUM,suffix)
    temp3 <- apply(calib_norm[freezeThaw,],2,sd)/apply(calib_norm[freezeThaw,],2,mean)
    
    tempList <- list(temp1,temp3)
    
  }else{
    suffix = "/8"
    acrossNUM = c(5,1,3,2,4,6)
    acrossPlates <- calibIDs %in% paste0(clinicType," POOL-HT-",acrossNUM,suffix)
    temp1 <- apply(calib_norm[acrossPlates,],2,sd)/apply(calib_norm[acrossPlates,],2,mean)
    
    freezeThawNUM = c(7) ### only one freeze thawed unspun sample and no unspun control on the same plate.
    freezeThaw <-  calibIDs %in% paste0(clinicType," POOL-HT-",freezeThawNUM,suffix)
    temp3 <- as.matrix(abs(calib_norm[freezeThaw,] - apply(calib_norm[which(acrossPlates),],2,mean))/apply(calib_norm[acrossPlates,],2,mean))
    tempList <- list(temp3)
  }
  
  return(tempList)
}