# This script includes the functions for processing the raw RFU data, before quality assessment, to be called by QC.paper.202306.R. 
# Date: 26 JUN 2023
# Version: 3.0
# Author(s): name = Dr.Yun Deng
#            email = yun.deng@kennedy.ox.ac.uk
#            name = Dr.Luke jostins 
#            email = luke.jostins@kennedy.ox.ac.uk
#            name = Dr.Thomas Perry
#            email = thomas.perry@ndorms.ox.ac.uk

#__________________________________________________________________________________________________________
### function to read in RFU frame and extract protein meta data
### Input: adata file name
### Output: RawM -- RFU frame; ProMeta -- protein meta data matrix
### Note: file "SS-228545_v4.1.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.20220224.adat" is provided by SOMAlogic, which contains the information of plasma plate scaling vector
prepareRFU <- function(pathIn,adatFile){
  RawM <- SomaDataIO::read_adat(paste0(pathIn,adatFile))
  ProMeta <- data.frame(attributes(RawM)$Col.Meta)
  rownames(ProMeta) <- colnames(RawM)[which(colnames(RawM)=="seq.10000.28"):ncol(RawM)]
  Type <<- ProMeta$Type
  Dilution <<- ProMeta$Dilution
  PlateScale_Reference <<- attributes(SomaDataIO::read_adat(paste0(pathIn,"somascan/adat/SS-228545_v4.1.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.20220224.adat")))$Col.Meta$PlateScale_Reference
  return(list(RawM,ProMeta))
}


#__________________________________________________________________________________________________________
# function to perform user defined standardization.
# input: UserInput -- user defined standardization steps, a character string with "," separated. 
# RawM -- original protein expression level measures (RFU).
# output: MySoma -- RFU frame derived from selected standardization steps.
getMySoma <- function(UserInput,RawM){
  userSt = strsplit(UserInput,",")
  if(userSt[[1]] =="RawM"){MySoma = RawM
  }else{
    numStep = length(userSt) 
    Funlist = vector(mode="list",length=numStep)
    for (counterStep in 1:numStep){
      Funlist[[counterStep]] = get(userSt[[counterStep]])
    }
    MySoma = UserNorm(Funlist,RawM)}
  
  return(MySoma)
}

UserNorm <- function(Funlist,RawM){
  for (FunCounter in 1:length(Funlist)){
    f <- Funlist[[FunCounter]]
    MySoma = f(RawM)
    RawM = MySoma
  }
  return(MySoma)
}

#__________________________________________________________________________________________________________
# function to match the order after each standardisation step. We normalize based on each plate, so the row order of RFU frame may change after each normalization step.
# input: RawM -- RFU frame before the standardisation step. 
# MySomaTemp -- RFU frame after the standardisation step.
# output: MySoma -- order adjusted RFU frame with the consistent order with input RawM.
matchOrder <- function(RawM,MySomaTemp){
  
  rowOrderName = rownames(RawM)
  rowOrder = unlist(sapply(rowOrderName,function(x){which(rownames(MySomaTemp)==x)}))
  
  colOrderName = colnames(RawM)
  colOrder = unlist(sapply(colOrderName,function(x){which(colnames(MySomaTemp)==x)}))
  
  MySoma = MySomaTemp[rowOrder,colOrder]
  
  return(MySoma)
}

#__________________________________________________________________________________________________________
# function to perform hybridization normalization performed on one plate. 
# Input: RawM -- RFU frame before hybridization normalisation.
# output MySoma -- RFU frame after hybridization normalisation.
HYBNORM <- function(RawM){
  
  PlateIdUni = levels(factor(RawM$PlateId))
  
  Platelist = list()
  
  DatStartId <- which(colnames(RawM)=="seq.10000.28")  ###for calculation convenience, extract data zone only
  
  HybId = which(Type == "Hybridization Control Elution")
  
  for (plateCounter in 1:length(PlateIdUni)){
    
    PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])
    
    RawMS = RawM[PlateIdSg,] ### single plate
    
    datZone = RawMS[,DatStartId:ncol(RawMS)] ###RFUs zone only
    
    rCmedian = matrix(apply(datZone[,HybId],2,median),nrow=1)
    
    HybNorm1 = t(apply(datZone[,HybId],1,function(x){rCmedian/x}))
    
    rRmedian = apply(HybNorm1,1,median)
    
    HybNorm = apply(datZone,2,function(x){x*rRmedian})
    
    RawMSDone = cbind(RawMS[,1:(DatStartId-1)],HybNorm)
    
    Platelist[[plateCounter]] = RawMSDone 
  }
  
  MySoma = do.call(rbind.data.frame,Platelist)  
  
  MySoma <- matchOrder(RawM,MySoma) # always pay attention whether the transformation changed the row order
  
  return(MySoma)
}

#__________________________________________________________________________________________________________
# function to perform median normalization on QC plasma calibrators.
# input: RawM -- RFU frame contain the original QC plasma calibartors.
# output: MySoma -- RFU frame contain the QC plasma calibrators after median normalisation.
MIDNORMcali = function(RawM){ ###caliCase control which sample type to be applied MidNorm
  
  PlateIdUni = levels(factor(RawM$PlateId))
  
  Platelist = list()
  
  for (plateCounter in 1:length(PlateIdUni)){
    
    PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])
    
    RawMS = RawM[PlateIdSg,] ### single plate
    DatStartIdP = which(colnames(RawMS) == "seq.10000.28")
    
    idHyb = DatStartIdP + which(Type == "Hybridization Control Elution") -1
    idNonHyb = c(1:(DatStartIdP-1),DatStartIdP + which(Type != "Hybridization Control Elution") -1)
    RawMS1 = RawMS[,idNonHyb] ###RawM1 to track "Ratio of Normalization Median to Sample Value"
    DatStartId = which(colnames(RawMS1) =="seq.10000.28")
    
    sampType = c("Calibrator")
    
    for (sampTypeCounter in 1:length(sampType)){
      
      idSamp = which(RawMS1$SampleType == sampType[sampTypeCounter])
      
      datZoneP = RawMS1[,DatStartId:ncol(RawMS1)] ### datazone excludes "HybControlElution"
      
      datZone = RawMS1[idSamp,DatStartId:ncol(RawMS1)]  ###dataZone only include interested sample type
      
      if(length(idSamp)==1) {SampTypeRFU = matrix(datZone,nrow=1)
      SampTypeMedian = SampTypeRFU
      SampMedianNorm = as.matrix(SampTypeMedian/SampTypeRFU)
      }else {SampTypeRFU = datZone
      SampTypeMedian = apply(SampTypeRFU,2,median)
      SampMedianNorm = t(apply(SampTypeRFU,1,function(x){SampTypeMedian/x}))}
      
      Dilute = Dilution[which(Dilution!="0")] ###SampMedianNorm & Dilute: ncol(SampMedianNorm)=length(Dilute)
      uniqDilute = levels(factor(Dilute))
      
      for (idDilute in (1:length(uniqDilute))){
        
        DataDiluteID = which(Dilute==uniqDilute[idDilute])
        
        MedianDiluteSingle = matrix(apply(SampMedianNorm[,DataDiluteID],1,median),ncol=1)
        
        DataDiluteNormT = apply(datZone[,DataDiluteID],2,function(x){MedianDiluteSingle*x})
        
        if (idDilute==1){DataDiluteNorm = DataDiluteNormT}
        else{DataDiluteNorm = cbind(DataDiluteNorm,DataDiluteNormT)}
      }
      
      if (sampTypeCounter==1) {RawMStemp = cbind(RawMS[idSamp,1:(DatStartIdP-1)],RawMS[idSamp,idHyb],DataDiluteNorm)
      }else{
        RawMStempT = cbind(RawMS[idSamp,1:(DatStartIdP-1)],RawMS[idSamp,idHyb],DataDiluteNorm)
        RawMStemp = rbind(RawMStemp,RawMStempT)
      }
      RawMStemp2 = rbind(RawMStemp,RawMS[which(RawMS$SampleType!="Calibrator"),])
    }
    
    Platelist[[plateCounter]] = RawMStemp2
  }
  
  MySomaTemp = do.call(rbind.data.frame,Platelist)
  
  MySoma <- matchOrder(RawM,MySomaTemp) # always pay attention whether the transformation changed the row order
  
  return(MySoma)
}

#__________________________________________________________________________________________________________
# function to perform plate scaling.
# input: RawM -- RFU frame before plate scaling.
# output: MySoma -- RFU frame after plate scaling.
PLATESCALE <- function(RawM){
  
  PlateIdUni = levels(factor(RawM$PlateId))
  
  Platelist = list()
  
  DatStartId <- which(colnames(RawM)=="seq.10000.28")  ###for calculation convenience, extract data zone only
  
  for (plateCounter in 1:length(PlateIdUni)){
    
    PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])
    
    RawMS = RawM[PlateIdSg,] ### single plate
    
    idCaliborator = which(RawMS$SampleType=="Calibrator")
    
    datZone = RawMS[,DatStartId:ncol(RawMS)]
    
    CaliboratorM = datZone[idCaliborator,]
    
    calibratorMedian = apply(CaliboratorM,2,median)
    
    PlateScaleRatio = PlateScale_Reference/calibratorMedian
    
    PlateScaleScalar = median(PlateScaleRatio)
    
    datZone2 = datZone*PlateScaleScalar
    
    RawMSDone = cbind(RawMS[,1:(DatStartId-1)],datZone2)
    
    Platelist[[plateCounter]] = RawMSDone 
  }
  
  MySoma = do.call(rbind.data.frame,Platelist)  
  
  MySoma <- matchOrder(RawM,MySoma) # always pay attention whether the transformation changed the row order
  
  return(MySoma)
}

#__________________________________________________________________________________________________________
# function to perform plate calibration.
# input: RawM -- RFU frame before plate calibration.
# output: MySoma -- RFU frame after plate calibation.
CALIBRATION <- function(RawM){
  
  PlateIdUni = levels(factor(RawM$PlateId))
  
  Platelist = list()
  
  DatStartId <- which(colnames(RawM)=="seq.10000.28")  ###for calculation convenience, extract data zone only
  
  for (plateCounter in 1:length(PlateIdUni)){
    
    PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])
    
    RawMS = RawM[PlateIdSg,] ### single plate
    
    idCaliborator = which(RawMS$SampleType=="Calibrator")
    
    datZone = RawMS[,DatStartId:ncol(RawMS)]
    
    CaliboratorM = datZone[idCaliborator,]
    
    CalibratorMedian = apply(CaliboratorM,2,median)
    
    CalSet=PlateScale_Reference/CalibratorMedian
    
    datZone2 = t(apply(datZone,1,function(x) x*CalSet))
    
    RawMSDone = cbind(RawMS[,1:DatStartId-1],datZone2)
    
    Platelist[[plateCounter]] = RawMSDone 
  }
  
  MySoma = do.call(rbind.data.frame,Platelist)  
  
  MySoma <- matchOrder(RawM,MySoma) # always pay attention whether the transformation changed the row order
  
  return(MySoma)
}

#__________________________________________________________________________________________________________
# function to perform median normalisation on samples
# input: RawM -- RFU frame before plate calibration.
# output: MySoma -- RFU frame after plate median normalisation on samples
MIDNORMsamp = function(RawM){ 
  
  PlateIdUni = levels(factor(RawM$PlateId))
  
  Platelist = list()
  
  for (plateCounter in 1:length(PlateIdUni)){
    
    PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])
    
    RawMS = RawM[PlateIdSg,] ### single plate
    DatStartIdP = which(colnames(RawMS) == "seq.10000.28")
    
    idHyb = DatStartIdP + which(Type == "Hybridization Control Elution") -1
    idNonHyb = c(1:(DatStartIdP-1),DatStartIdP + which(Type != "Hybridization Control Elution") -1)
    RawMS1 = RawMS[,idNonHyb] ###RawM1 to track "Ratio of Normalization Median to Sample Value"
    DatStartId = which(colnames(RawMS1) =="seq.10000.28")
    
    sampType = c("Sample")
    
    for (sampTypeCounter in 1:length(sampType)){
      
      idSamp =grep(sampType[sampTypeCounter],RawMS1$SampleType) 
      
      datZoneP = RawMS1[,DatStartId:ncol(RawMS1)] ### datazone excludes "HybControlElution"
      
      datZone = RawMS1[idSamp,DatStartId:ncol(RawMS1)]  ###dataZone only include interested sample type
      
      if(length(idSamp)==1) {SampTypeRFU = matrix(as.matrix(datZone),nrow=1)
      SampTypeMedian = SampTypeRFU
      SampMedianNorm = as.matrix(SampTypeMedian/SampTypeRFU)
      }else {SampTypeRFU = datZone
      SampTypeMedian = apply(SampTypeRFU,2,median)
      SampMedianNorm = t(apply(SampTypeRFU,1,function(x){SampTypeMedian/x}))}
      
      Dilute = Dilution[which(Dilution!="0")] ###SampMedianNorm & Dilute: ncol(SampMedianNorm)=length(Dilute)
      uniqDilute = levels(factor(Dilute))
      
      for (idDilute in (1:length(uniqDilute))){
        
        DataDiluteID = which(Dilute==uniqDilute[idDilute])
        
        MedianDiluteSingle = matrix(apply(SampMedianNorm[,DataDiluteID],1,median),ncol=1)
        
        DataDiluteNormT = apply(datZone[,DataDiluteID],2,function(x){MedianDiluteSingle*x})
        
        if (idDilute==1){DataDiluteNorm = DataDiluteNormT
        }else{DataDiluteNorm = cbind(DataDiluteNorm,DataDiluteNormT)}
      }
      
      if (sampTypeCounter==1) {RawMStemp = cbind(RawMS[idSamp,1:(DatStartIdP-1)],RawMS[idSamp,idHyb],DataDiluteNorm)
      }else{
        RawMStempT = cbind(RawMS[idSamp,1:(DatStartIdP-1)],RawMS[idSamp,idHyb],DataDiluteNorm)
        RawMStemp = rbind(RawMStemp,RawMStempT)
      }
      RawMStemp2 = rbind(RawMStemp,RawMS[which(RawMS$SampleType =="Calibrator"),])
    }
    
    Platelist[[plateCounter]] = RawMStemp2
  }
  
  MySomaTemp = do.call(rbind.data.frame,Platelist)
  
  MySoma <- matchOrder(RawM,MySomaTemp) # always pay attention whether the transformation changed the row order
 
  return(MySoma)
}
