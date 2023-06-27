### This script is specialized for making plots and tables in QC paper
### Before running this script, corresponding R objects need to be loaded, which were stored by running main script QC.paper.2023.06.R.
# Date: 26 JUN 2023
# Version: 3.0
# Author(s): name = Dr.Yun Deng
#            email = yun.deng@kennedy.ox.ac.uk
#            name = Dr.Luke jostins 
#            email = luke.jostins@kennedy.ox.ac.uk
#            name = Dr.Thomas Perry
#            email = thomas.perry@ndorms.ox.ac.uk

library(ggplot2)
library(cowplot)
library(GGally)
library(ggpubr)
library(factoextra)
library(ggforce)
library(scales)

myCodeIn <- "/Users/ydeng/Documents/QCpaper.Code/"  
pathIn <- "/Users/ydeng/Documents/QCpaper.Code/STEPUP_DAG_rel002_discovery1/"
Robj.Path <- "/Users/ydeng/Documents/QCpaper.Code/Robj.Paper/"

source(paste0(myCodeIn,"QCassess.Paper.202306.R"))
load(file=paste0(Robj.Path,"All.RUFs.Rdat"))
load(file=paste0(Robj.Path,"pcDat.Rdat"))
load(file=paste0(Robj.Path,"Figure1.Rdat")) 
load(file=paste0(Robj.Path,"Figure2.Rdat")) 
load(file=paste0(Robj.Path,"Figure3.Rdat")) 
load(file=paste0(Robj.Path,"Figure4.Rdat")) 
load(file=paste0(Robj.Path,"UmapDis.filtered.Rdat"))
load(file=paste0(Robj.Path,"CV.OA.Repeats.FreezeThaw.Reprocess.Rdat")) 
load(file=paste0(Robj.Path,"Filters.Rdat"))
load(file=paste0(Robj.Path,"TechVS10pc.Rdat"))

###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
### Assessment of mean % CV and  mean R2 of each protein, correlation coefficient with immunoassay, comparing  across different standardisations 
###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
StandardisationLabel = c("Raw Data","HN","HN + PS","HN + PS + MN","HN + PS + MN + PC","HN + PS + PC") ### use the abbreviations for each normasliation step 

p.meanCV <- ggplot(data = meanCV) + geom_point(aes(x=NormalisationSteps,y=100*as.numeric(meanCV),group=DiseaseGroup,col=DiseaseGroup)) + geom_line(aes(x=NormalisationSteps,y=100*as.numeric(meanCV),group=DiseaseGroup,col=DiseaseGroup)) + theme_bw() +
  xlab("") + ylab("\n\nMean %CV") + labs(color = "Disease Group") + scale_x_discrete(breaks=seq(1:length(StandardisationLabel)),labels=StandardisationLabel) + 
  theme(axis.text.x = element_text(size=10, angle=30,face="bold",hjust=1),axis.text.y = element_text(size=10,face="bold"),
        legend.title =element_text(size = 9,face="bold"), legend.text = element_text(size = 9,face="bold"),
        axis.title.y =element_text(size=10,face="bold"),axis.title.x =element_text(size=10)) 

p.meanR2 <- ggplot(data = meanR2) + geom_point(aes(x=NormalisationSteps,y=100*as.numeric(meanR2),group=DiseaseGroup,col=DiseaseGroup)) + geom_line(aes(x=NormalisationSteps,y=100*as.numeric(meanR2),group=DiseaseGroup,col=DiseaseGroup)) +
  xlab("") + ylab(bquote(atop("\n","" ~ bold("Mean R") ^ bold("2")))) + labs(color = "Disease Group") + scale_x_discrete(breaks=seq(1:length(StandardisationLabel)),labels=StandardisationLabel) + theme_bw() +
  theme(axis.text.x = element_text(size=10, angle=30,face="bold",hjust=1),axis.text.y = element_text(size=10,face="bold"),
        legend.title =element_text(size = 9,face="bold"), legend.text = element_text(size = 9,face="bold"),
        axis.title.y =element_text(size=10,face="bold"),axis.title.x =element_text(size=10)) 

p.immunoassay.OA <- ggplot(data = CorDatP.OA) + geom_line(aes(x=as.character(CorC),y=as.numeric(CorDatY),group=CorDatX,col=CorDatX)) + ylim(-0.75,1) + 
  xlab("") + ylab("Correlation Coefficient\nOA Samples") +labs(color = "Protein") + scale_x_discrete(breaks=seq(1:length(StandardisationLabel)),labels=StandardisationLabel) + scale_colour_manual(values = c("blue","seagreen","green","blueviolet","red","deepskyblue","hotpink","darkgoldenrod4","orange"))  + theme_bw() + 
  theme(axis.text.x = element_text(size=10, angle=30,face="bold",hjust=1),axis.text.y = element_text(size=10,face="bold"),
        legend.title =element_text(size = 9,face="bold"), legend.text = element_text(size = 9,face="bold"),
        axis.title.y =element_text(size=10,face="bold"),axis.title.x =element_text(size=10))

p.immunoassay.INJ <- ggplot(data = CorDatP.INJ) + geom_line(aes(x=as.character(CorC),y=as.numeric(CorDatY),group=CorDatX,col=CorDatX)) + ylim(-0.75,1) + 
  xlab("") + ylab("Correlation Coefficient\nInjury Samples") +labs(color = "Protein") + scale_x_discrete(breaks=seq(1:length(StandardisationLabel)),labels=StandardisationLabel) + scale_colour_manual(values = c("blue","seagreen","green","blueviolet","red","deepskyblue","hotpink","darkgoldenrod4","orange"))  + theme_bw() + 
  theme(axis.text.x = element_text(size=10, angle=30,face="bold",hjust=1),axis.text.y = element_text(size=10,face="bold"),
        legend.title =element_text(size = 9,face="bold"), legend.text = element_text(size = 9,face="bold"),
        axis.title.y =element_text(size=10,face="bold"),axis.title.x =element_text(size=10))

p1  <- ggpubr::ggarrange(p.meanCV,p.meanR2,labels = c("A", "B"),common.legend=TRUE,legend="right")
p2 <- ggpubr::ggarrange(p.immunoassay.OA,p.immunoassay.INJ,labels = c("C", "D"),common.legend=TRUE,legend="right")

### plot Figure1 in QC paper
cowplot::plot_grid(p1, p2,nrow=2,ncol=1,align = "v")

###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
### Investigate the PC1 drivers -- Intracellular Protein Score
###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
##########################################################################################################################################################
### Variation explained by top 10 PCs after standardisation and IPS adjustment
##########################################################################################################################################################
p.variation.PC.before <- ggplot() + geom_bar(aes(x=1:10,y=eig.val.standardised.all[1:10]),stat="identity",fill="steelblue") + scale_x_discrete(limits=paste0("PC",1:10)) +
  xlab("") + ylab("Variation Explained (%)") + 
  geom_text(aes(x=1:10,y=eig.val.standardised.all[1:10]+0.9,label=paste0(signif(eig.val.standardised.all[1:10],2),"%")),size=3,fontface = "bold") + theme_bw() + 
  theme(axis.text.x = element_text(size=7,face="bold"),axis.text.y =element_text(size=10,face="bold"),axis.title.y =element_text(size=12.5,face="bold",vjust=-1),
        panel.grid.major.x = element_blank()) 

p.variation.PC.after  <- ggplot() + geom_bar(aes(x=1:10,y=eig.val.batchDone2[1:10]),stat="identity",fill="steelblue") + scale_x_discrete(limits=paste0("PC",1:10)) +
  xlab("") + ylab("Variation Explained (%)\nafter IPS Adjustment") + theme_bw() + 
  geom_text(aes(x=1:10,y=eig.val.batchDone2[1:10]+0.25,label=paste0(signif(eig.val.batchDone2[1:10],2),"%")),size=3,fontface = "bold") + theme_bw() + 
  theme(axis.text.x = element_text(size=7,face="bold"),axis.text.y =element_text(size=10,face="bold"),axis.title.y =element_text(size=12.5,face="bold",vjust=-1),
        panel.grid.major.x = element_blank()) 

##########################################################################################################################################################
### Visualization of PC1 driver -- protein abundance 
##########################################################################################################################################################
abundance.PC1 <- data.frame("abunadance"=log(corPerPro.BatchCorrected[[2]]),"correlation"=corPerPro.BatchCorrected[[1]])
albminSeq <- rownames(ProMeta)[which(ProMeta$EntrezGeneSymbol=="ALB")]
LDHsea <- rownames(ProMeta)[which(ProMeta$EntrezGeneSymbol=="LDHB")]
p.abundance.PC1 <- ggplot(abundance.PC1,aes(x=abunadance,y=correlation)) + geom_point(color="royalblue4",size=0.5) + geom_smooth(method='lm',color="black") +
  geom_point(aes(x=abundance.PC1[albminSeq,"abunadance"],y=abundance.PC1[albminSeq,"correlation"]),color="#009E73",shape=24,size=3) +
  annotate(geom="text",x=abundance.PC1[albminSeq,"abunadance"]+0.1, y=abundance.PC1[albminSeq,"correlation"]-0.15, label="Albumin",color="#009E73",fontface="bold") +
  geom_point(aes(x=abundance.PC1[LDHsea,"abunadance"],y=abundance.PC1[LDHsea,"correlation"]),color="#D55E00",shape=24,size=3) +
  annotate(geom="text",x=abundance.PC1[LDHsea,"abunadance"]+2, y=abundance.PC1[LDHsea,"correlation"], label="LDH",color="#D55E00",fontface="bold") +
  xlab("Log Mean Protein Abudance") + ylab("Correlation Coefficient with PC1") + theme_bw() + 
  theme(axis.title.x = element_text(size=10,face="bold"),axis.text.x = element_text(size=10,face="bold"),
        axis.title.y =element_text(size=12.5,face="bold",vjust=-1),axis.text.y = element_text(size=10,face="bold"),
        legend.position = c(0.8,0.1),legend.title =element_text(size = 11,face="bold"), legend.text = element_text(size = 11,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

##########################################################################################################################################################
### PC1 of 18 paired spun/unspun samples 
##########################################################################################################################################################
adjustY <- rep(0,36)
adjustY[c(8,17,1,22,36,29,19,26,32,30,25,33,28,20,31)] <- c(4,2,1,-4,-3.5,-2,-3,-1.3,-6,-2,5,7,11,7,5)
p.spin.PC1 <- ggplot(data=spinFrame,aes(x=Spin,y=as.numeric(PC1))) + geom_point(aes(color=Spin),size=2) + geom_line(aes(group=rep(1:18,2),color=rep(spin.col,2))) +
  geom_text(aes(x=c(rep(0.7,18),rep(2.3,18)),y=as.numeric(PC1)-adjustY,label = c(sub("-V.*$","",rownames(spinFrame)))),size=2) + xlab("") + ylab("PC1") + scale_x_discrete(limits=c("Unspun","Spun")) + theme_bw() +
  theme(axis.title.x = element_text(size=10,face="bold"),axis.text.x = element_text(size=10,face="bold"), axis.text.y = element_text(size=10,face="bold"),legend.position = "none",
        axis.title.y =element_text(size=12.5,face="bold",vjust=-1),panel.grid.major.x = element_blank()) 

##########################################################################################################################################################
### Intracellular protein score vs PC1 on the IPS adjusted and non-adjusted data
##########################################################################################################################################################
p.RegProBefore <- ggplot() + geom_point(aes(x=pcDat.batchDone[,"PC1"],y=IPScore)) + geom_smooth(aes(x=pcDat.batchDone[,"PC1"],y=IPScore),method='lm',color="red") + 
  xlab("PC1 before IPS Adjustment") + ylab("Intracellular Protein Score") + annotate(geom="text",x=60, y=35000, label=paste0("Correlation Coefficient: ",signif(cor(pcDat.batchDone[,"PC1"],IPScore),2)),color="red",fontface = "bold",size=3) + theme_bw() + 
  theme(axis.text.x = element_text(size=10,face="bold"),axis.title.x =element_text(size = 10,face="bold"), 
        axis.text.y =element_text(size=10,face="bold"),axis.title.y =element_text(size=12.5,face="bold",vjust=-1)) 

p.RegProAfter <- ggplot() + geom_point(aes(x=pcDat.batchDone2[,"PC1"],y=IPScore)) + geom_smooth(aes(x=pcDat.batchDone2[,"PC1"],y=IPScore),method='lm',color="red") + 
  xlab("PC1 after IPS Adjustment") + ylab("Intracellular Protein Score") + annotate(geom="text",x=-50, y=35000, label=paste0("Correlation Coefficient: ",signif(cor(pcDat.batchDone2[,"PC1"],IPScore),2)),color="red",fontface = "bold",size=3) + theme_bw() + 
  theme(axis.text.x = element_text(size=10,face="bold"),axis.title.x =element_text(size = 10,face="bold"), 
        axis.text.y =element_text(size=10,face="bold"),axis.title.y =element_text(size=12.5,face="bold",vjust=-1)) 

### plot Figure2 in QC paper
cowplot::plot_grid(p.variation.PC.before, p.abundance.PC1,p.spin.PC1,p.RegProBefore,p.RegProAfter,p.variation.PC.after,labels = c("A", "B","C","D","E","F"),nrow=2,ncol=3,align = "v")

##########################################################################################################################################################
### PC1 driver regression model
##########################################################################################################################################################
SubLocation <- read.csv(paste0(pathIn,"online resource/subcellular_location.tsv"),sep="\t") 
LocationList <- diffProCombineOnline(SubLocation,"Main.location")
SubLocationDat = LocationList[[1]]
SubLocationType = LocationList[[2]]
SublocationSet = LocationList[[3]]
CytoplasmL <- read.csv(paste0(pathIn,"online resource/Cytoplasm.txt"))[,1]
NucleusL <- read.csv(paste0(pathIn,"online resource/Nucleus.txt"))[,1]
EndomembraneL <- read.csv(paste0(pathIn,"online resource/Endomembrane.txt"))[,1]
secreted <- SubLocation$Gene.name[which(SubLocation$Extracellular.location=="Predicted to be secreted")]

# for each protein map their detailed cell location to the broad cell location types (Cytoplasm, Nucleus, Endomembrane, and Secreted) 
for (indexCounter in 1:nrow(SubLocationDat)){
  x = SubLocationDat$Main.location[indexCounter]
  pro = SubLocationDat$Gene.name[indexCounter]
  mainLoc = strsplit(x,";")[[1]]
  if(any(CytoplasmL %in% mainLoc) & !(any(pro %in% secreted))){newLocation = "Cytoplasm"
  }else if(any(NucleusL %in% mainLoc) & !(any(pro %in% secreted))){newLocation = "Nucleus"
  }else if(any(EndomembraneL %in% mainLoc) & !(any(pro %in% secreted))){newLocation = "Endomembrane"
  }else if(pro %in% secreted){newLocation = "secreted"
  }else{newLocation = NA}
  SubLocationDat$Broad.location[indexCounter]= newLocation
}

BroadLocationList <- diffProCombineOnline(SubLocationDat,"Broad.location")
SubLocationBroad = BroadLocationList[[1]]
SubLocationBroadType = BroadLocationList[[2]]
SubLocationBroadSet = BroadLocationList[[3]]

keepseq <- rownames(ProMeta)[which(ProMeta$Organism=="Human" & ProMeta$Type=="Protein")]
Nucleus <- sapply(keepseq,function(x) {ifelse(any(ProMeta[x,"EntrezGeneSymbol"] %in% SubLocationBroadSet$Nucleus),1,0)})
Cytoplasm <- sapply(keepseq,function(x) {ifelse(any(ProMeta[x,"EntrezGeneSymbol"] %in% SubLocationBroadSet$Cytoplasm),1,0)})
Endomembrane <- sapply(keepseq,function(x) {ifelse(any(ProMeta[x,"EntrezGeneSymbol"] %in% SubLocationBroadSet$Endomembrane),1,0)})
secretom <- sapply(keepseq,function(x) {ifelse(any(ProMeta[x,"EntrezGeneSymbol"] %in% SubLocationBroadSet$secreted),1,0)})

# for each protein check whether it is a marker for monocyte, neutrophil, macrophage
Tissue.Cell.Raw <- read.table(paste0(pathIn,"online resource/PanglaoDB_markers_27_Mar_2020.tsv"),h=T,sep="\t",quote="")
Tissue.Cell.Raw <- Tissue.Cell.Raw[grep("Hs",Tissue.Cell.Raw$species),] # keep only human information
TissueList <- na.omit(unique(Tissue.Cell.Raw$organ)) # extract all the tissue names

Tissue.Cell.GeneSets1 <- vector("list",length=length(TissueList)) # generate the tissue cell gene sets 
for(tsCount in 1:length(TissueList)){
  Tissue <- TissueList[tsCount]
  TissueFrame <- Tissue.Cell.Raw[which(Tissue.Cell.Raw$organ==Tissue),c("official.gene.symbol","cell.type")]
  Tissue.CellType <- unique(TissueFrame$cell.type)
  Tissue.GeneSets <- sapply(Tissue.CellType,function(x){TissueFrame$official.gene.symbol[which(TissueFrame$cell.type==x)]})
  Tissue.Cell.GeneSets1[[tsCount]] <- Tissue.GeneSets
  names(Tissue.Cell.GeneSets1)[tsCount] <- Tissue
}

MonoCyte <- sapply(keepseq,function(x) {ifelse(any(ProMeta[x,"EntrezGeneSymbol"] %in% Tissue.Cell.GeneSets1$`Immune system`$Monocytes),1,0)})
Neutrophil <- sapply(keepseq,function(x) {ifelse(any(ProMeta[x,"EntrezGeneSymbol"] %in% Tissue.Cell.GeneSets1$`Immune system`$Neutrophils),1,0)})
Macrophage <- sapply(keepseq,function(x) {ifelse(any(ProMeta[x,"EntrezGeneSymbol"] %in% Tissue.Cell.GeneSets1$`Immune system`$Macrophages),1,0)})

# regression model with broad subcellular location as predictors, correlation coefficient between protein and PC1 as response variable, both on batch corrected/non-IPS adjusted and batch corrected/IPS-adjusted data
summary(lm(corPerPro.BatchCorrected[[1]][keepseq] ~ log(corPerPro.BatchCorrected[[2]][keepseq])))
summary(lm(corPerPro.BatchCorrected[[1]][keepseq] ~ as.factor(Nucleus) + log(corPerPro.BatchCorrected[[2]][keepseq])))
summary(lm(corPerPro.BatchCorrected[[1]][keepseq] ~ as.factor(secretom) + log(corPerPro.BatchCorrected[[2]][keepseq])))
summary(lm(corPerPro.BatchCorrected[[1]][keepseq] ~ as.factor(MonoCyte) + log(corPerPro.BatchCorrected[[2]][keepseq])))
summary(lm(corPerPro.BatchCorrected[[1]][keepseq] ~ as.factor(Neutrophil) + log(corPerPro.BatchCorrected[[2]][keepseq])))
summary(lm(corPerPro.BatchCorrected[[1]][keepseq] ~ as.factor(Macrophage) + log(corPerPro.BatchCorrected[[2]][keepseq])))

summary(lm(corPerPro.reg[[1]][keepseq] ~ log(corPerPro.reg[[2]][keepseq])))
summary(lm(corPerPro.reg[[1]][keepseq] ~ as.factor(Nucleus) + log(corPerPro.reg[[2]][keepseq])))
summary(lm(corPerPro.reg[[1]][keepseq] ~ as.factor(secretom) + log(corPerPro.reg[[2]][keepseq])))
summary(lm(corPerPro.reg[[1]][keepseq] ~ as.factor(MonoCyte) + log(corPerPro.reg[[2]][keepseq])))
summary(lm(corPerPro.reg[[1]][keepseq] ~ as.factor(Neutrophil) + log(corPerPro.reg[[2]][keepseq])))
summary(lm(corPerPro.reg[[1]][keepseq] ~ as.factor(Macrophage) + log(corPerPro.reg[[2]][keepseq])))

### TableS5 in QC paper shows the summary of (PC1 driver) results from the above regression model

###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
### Investigate the PC2 drivers -- bimodal signal 
###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
##########################################################################################################################################################
### Histogram along PC2
##########################################################################################################################################################
p.bimodal.PC2.noColor <- ggplot(data=data.frame(pcDat.standardised[,c(1,2)])) + geom_histogram(aes(x=PC2,y=..density..),bins = 100,fill="white",color="steelblue3") + 
  mapply(function(mean, sd, lambda) {
    stat_function(
      fun = function(x) {
        (dnorm(x, mean = mean, sd = sd)) * lambda}, size=1)},
    mean = gmm_fit[["mu"]], 
    sd = gmm_fit[["sigma"]],
    lambda = gmm_fit[["lambda"]]) + ylab("Density") + theme_bw() + 
  theme(axis.title.x = element_text(size=11,face="bold"),axis.text.x = element_text(size=10,face="bold"), 
        axis.title.y =element_text(size=11,face="bold"),axis.text.y = element_text(size=10,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p.bimodal.PC2.Before <- ggplot(data=data.frame(pcDat.standardised[,c(1,2)])) + geom_histogram(aes(x=PC2,y=..density..),bins = 100,fill="white",color="steelblue3") +
  stat_function(fun = function(x){(dnorm(x, mean = gmm_fit[["mu"]][[1]], sd = gmm_fit[["sigma"]][[1]])) * gmm_fit[["lambda"]][[1]]}, colour = "#D55E00", size=1) +
  stat_function(fun = function(x){(dnorm(x, mean = gmm_fit[["mu"]][[2]], sd = gmm_fit[["sigma"]][[2]])) * gmm_fit[["lambda"]][[2]]}, colour = "#009E73", size=1) + 
  xlab("PC2") + ylab("Density") + theme_bw() +
  theme(axis.title.x = element_text(size=11,face="bold"),axis.text.x = element_text(size=10,face="bold"), 
        axis.title.y =element_text(size=11,face="bold"),axis.text.y = element_text(size=10,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p.bimodal.PC2.After <- ggplot(data=data.frame(pcDat.batchDone[,c(1,2)])) + geom_histogram(aes(x=PC2,y=..density..),bins = 100,fill="white",color="steelblue3") + 
  stat_function(fun = function(x) {dnorm(x,mean=mean(pcDat.batchDone[,2]),sd=sd(pcDat.batchDone[,2]))},color="#D55E00",size=1) + 
  stat_function(fun = function(x) {dnorm(x,mean=mean(pcDat.batchDone[,2]),sd=sd(pcDat.batchDone[,2]))},color="#009E73",linetype="dashed",size=1) + 
  xlab("PC2") + ylab("Density") + theme_bw() +
  theme(axis.title.x = element_text(size=11,face="bold"),axis.text.x = element_text(size=10,face="bold"), 
        axis.title.y =element_text(size=11,face="bold"),axis.text.y = element_text(size=10,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##########################################################################################################################################################
### UMAP before and after batch correction
##########################################################################################################################################################
p.umap.BimodalBefore.noColor <- ggplot() + geom_point(aes(x=myUmap1$layout[,1],y=myUmap1$layout[,2]),size=0.3) +
  xlab("Dimension1") + ylab("Dimension2") + theme_bw() + 
  theme(axis.title.x = element_text(size=11,face="bold"),axis.text.x = element_text(size=10,face="bold"),
        axis.title.y =element_text(size=11,face="bold"), axis.text.y = element_text(size=10,face="bold"))

p.umap.BimodalBefore <- ggplot() + geom_point(aes(x=myUmap1$layout[,1],y=myUmap1$layout[,2],color=as.factor(bimodalLabel)),size=0.3) +
  xlab("Dimension1") + ylab("Dimension2") + labs(color="Bimodal Signal Status") + scale_colour_manual(values = c("#D55E00","#009E73"),labels = c("Status1", "Status2")) + 
  guides(color = guide_legend(override.aes = list(size = 3))) + theme_bw() + 
  theme(axis.title.x = element_text(size=11,face="bold"),axis.text.x = element_text(size=10,face="bold"),
        axis.title.y =element_text(size=11,face="bold"), axis.text.y = element_text(size=10,face="bold"),
        legend.position="bottom", legend.title =element_text(size = 9,face="bold"), legend.text = element_text(size = 9,face="bold"),
        legend.margin=margin(0,0,0,0),legend.box.margin=margin(-4,-4,-4,-4))

p.umap.BimodalAfter <- ggplot() + geom_point(aes(x=myUmap2$layout[,1],y=myUmap2$layout[,2],color=as.factor(bimodalLabel)),size=0.3) +
  xlab("Dimension1") + ylab("Dimension2") + labs(color="Bimodal Signal Status") + scale_colour_manual(values = c("#D55E00","#009E73"),labels = c("Status1", "Status2")) + 
  guides(color = guide_legend(override.aes = list(size = 3))) + theme_bw() +
  theme(axis.title.x = element_text(size=11,face="bold"),axis.text.x = element_text(size=10,face="bold"),
        axis.title.y =element_text(size=11,face="bold"), axis.text.y = element_text(size=10,face="bold"),
        legend.position="bottom", legend.title =element_text(size = 9,face="bold"), legend.text = element_text(size = 9,face="bold"),
        legend.margin=margin(0,0,0,0),legend.box.margin=margin(-4,-4,-4,-4))

##########################################################################################################################################################
### TSG101 vs processing order, processing batch; and for reprocessed samples
##########################################################################################################################################################
plotData <- data.frame(ProcessingOrder = CombinedFrameKeep$sf_iknee_proc_order, TSG101 = exprDat_normKeep[,TSGseq], ProcessingBatch = CombinedFrameKeep$sf_iknee_proc_batch)

cols <- rep(RColorBrewer::brewer.pal(9, "Set1"),100)[1:length(levels(plotData$ProcessingBatch))]
names(cols) <- sort(as.numeric(levels(plotData$ProcessingBatch)))

p.TSGseq <- ggplot(plotData,aes(x=ProcessingOrder,y=TSG101,color=ProcessingBatch,group=ProcessingBatch)) + geom_point(size=0.5) + 
  geom_line() + scale_colour_manual(values = cols) +
  xlab("Processing Order") + ylab("TSG101") + theme_bw() + 
  theme(legend.position="none",
        axis.title.x = element_text(size=11,face="bold"),axis.text.x = element_text(size=10,face="bold"),
        axis.title.y =element_text(size=11,face="bold"), axis.text.y = element_text(size=10,face="bold"))

p.reprocess <- ggplot(data=reprocessFrame) + geom_point(aes(x=Sample,y=as.numeric(TSG101),group=Processing,color=Processing)) +
  geom_line(aes(x=Sample,y=as.numeric(TSG101),group=Processing,color=Processing)) + labs(color="") + xlab("") + ylab("TSG101") + theme_bw() + 
  theme(axis.text.x = element_text(size=8,face="bold"),
        axis.title.y =element_text(size=10,face="bold"), axis.text.y = element_text(size=10,face="bold"),
        legend.position="top",legend.title =element_text(size = 8,face="bold"), legend.text = element_text(size = 8,face="bold"),
        legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-10,-10))

##########################################################################################################################################################
### plot Figure3 in QC paper
p3 <- ggpubr::ggarrange(p.bimodal.PC2.noColor, p.umap.BimodalBefore.noColor, p.TSGseq, p.reprocess, labels = c("A", "B","C","D"))
p4  <- ggpubr::ggarrange(p.bimodal.PC2.Before, p.umap.BimodalBefore, p.bimodal.PC2.After, p.umap.BimodalAfter, labels = c("E", "F", "G", "H"),common.legend=TRUE,legend="bottom")

cowplot::plot_grid(p3,p4,nrow=2,ncol=1,align = "v")

###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
### Correlation coefficient between SOMAscan and immunoassay, comparison across raw, standardized, non-IPS adjusted and IPS adjusted data
###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
NormalisationLabel = c("Raw Data","Optimized Standardisation","Processed without IPS Adjustment","Processed with IPS Adjustment")

p.immunoassay.OA.IPS <- ggplot(data = CorDatP.OA.IPS) + geom_line(aes(x=as.character(CorC),y=as.numeric(CorDatY),group=CorDatX,col=CorDatX)) + ylim(-0.5,1) + scale_colour_manual(values = c("blue","seagreen","green","blueviolet","red","deepskyblue","hotpink","darkgoldenrod4","orange")) +
  xlab("") + ylab("Correlation Coefficient for OA Samples") +labs(color = "Protein") + scale_x_discrete(breaks=seq(1:length(NormalisationLabel)),labels=NormalisationLabel) + theme_bw() + 
  theme(axis.text.x = element_text(size=10, angle=20,face="bold",hjust=1),axis.text.y = element_text(size=10,face="bold"),
        legend.title =element_text(size = 10,face="bold"), legend.text = element_text(size = 10,face="bold"),
        axis.title.y =element_text(size=10,face="bold"),axis.title.x =element_text(size=10))

p.immunoassay.INJ.IPS <- ggplot(data = CorDatP.INJ.IPS) + geom_line(aes(x=as.character(CorC),y=as.numeric(CorDatY),group=CorDatX,col=CorDatX)) + ylim(-0.5,1) + scale_colour_manual(values = c("blue","seagreen","green","blueviolet","red","deepskyblue","hotpink","darkgoldenrod4","orange")) +
  xlab("") + ylab("Correlation Coefficient for Injury Samples") +labs(color = "Protein") + scale_x_discrete(breaks=seq(1:length(NormalisationLabel)),labels=NormalisationLabel) + theme_bw() + 
  theme(axis.text.x = element_text(size=10, angle=20,face="bold",hjust=1),axis.text.y = element_text(size=10,face="bold"),
        legend.title =element_text(size = 10,face="bold"), legend.text = element_text(size = 10,face="bold"),
        axis.title.y =element_text(size=10,face="bold"),axis.title.x =element_text(size=10))

### plot Figure4 in QC paper
ggpubr::ggarrange(p.immunoassay.OA.IPS,p.immunoassay.INJ.IPS,labels=c("A","B"),ncol=2,nrow=1,common.legend=TRUE,legend="bottom")

### Correlation coefficient between intracellular protein score and 9 immunoassay protein measures, stratified by OA and injury, as shown in TableS6
KeepOA <- unlist(sapply(sandwich_master_Ben$PIN,function(x){names(IPScore)[grep(x,names(IPScore))]}))
sapply(Test1,function(x) {cor.test(as.numeric(sandwich_master_Ben[,x]),IPScore[KeepOA],use="pairwise.complete.obs")$estimate})
sapply(Test1,function(x) {cor.test(as.numeric(sandwich_master_Ben[,x]),IPScore[KeepOA],use="pairwise.complete.obs")$p.value})

KeepINJ <- unlist(sapply(sandwich_master_Historic$PIN,function(x){names(IPScore)[grep(x,names(IPScore))]}))
sapply(Test1[-6],function(x) {cor.test(as.numeric(sandwich_master_Historic[,x]),IPScore[KeepINJ],use="pairwise.complete.obs")$estimate})
sapply(Test1[-6],function(x) {cor.test(as.numeric(sandwich_master_Historic[,x]),IPScore[KeepINJ],use="pairwise.complete.obs")$p.value})

###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
### Technical confounder visualisation on the most associated two PCs for batch corrected with/without filtering, batch corrected + IPS adjusted with/without filtering
###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
### plot Figure6 and FigureS6 in QC paper
plotTechPC(CombinedFrame,pcDat.batchDone)  ### batch corrected, non-IPS adjusted, non-filtered
plotTechPC(CombinedFrame,pcDat.batchDone2) ### batch corrected, IPS adjusted, non-filtered
plotTechPC(CombinedFrame,pcDat.batchCorrected.Filtered) ### batch corrected, non-IPS adjusted, filtered
plotTechPC(CombinedFrame,pcDat.IPSreg.Filtered) ### batch corrected, IPS adjusted, filtered

###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
### UMAP visualisation on filtered data for non-IPS adjusted and IPS adjusted data
###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
DiseaseSize1=ifelse(myUmap.BC.final.F$DiseaseGroup==4,"Big",ifelse(myUmap.BC.final.F$DiseaseGroup==3,"Mid","Small"))
p.umap.BC.Dis <- ggplot(data=myUmap.BC.final.F) + geom_point(aes(x=D1,y=D2,color=as.character(DiseaseGroup),size=DiseaseSize1)) +
  xlab("D1") + ylab("D2") + scale_colour_manual(name="IPS Adjusted", values = c("#D55E00","#009E73","black","purple"),labels=c('OA', 'Injury', 'Healthy control',"Inflammatory control")) + 
  scale_size_manual (values= c(2,2,1.5)) + guides(size=FALSE,color = guide_legend(override.aes = list(size = 3))) + theme_bw() + 
  theme(axis.title.x = element_text(size=12,face="bold"),axis.text.x = element_text(size=12,face="bold"),
        axis.title.y =element_text(size=12,face="bold"), axis.text.y = element_text(size=12,face="bold"),
        legend.title =element_text(size = 13,face="bold"), legend.text = element_text(size = 12,face="bold"))

DiseaseSize2=ifelse(myUmap.IPS.final.F$DiseaseGroup==4,"Big",ifelse(myUmap.IPS.final.F$DiseaseGroup==3,"Mid","Small"))
p.umap.IPS.Dis <- ggplot(data=myUmap.IPS.final.F) + geom_point(aes(x=D1,y=D2,color=as.character(DiseaseGroup),size=DiseaseSize2)) +
  xlab("D1") + ylab("D2") + scale_colour_manual(name="After IPS Adjustment", values = c("#D55E00","#009E73","black","purple"),labels=c('OA', 'Injury', 'Healthy control',"Inflammatory control")) + 
  scale_size_manual (values= c(2,2,1.5)) + guides(size=FALSE,color = guide_legend(override.aes = list(size = 3))) + theme_bw() + 
  theme(axis.title.x = element_text(size=12,face="bold"),axis.text.x = element_text(size=12,face="bold"),
        axis.title.y =element_text(size=12,face="bold"), axis.text.y = element_text(size=12,face="bold"),
        legend.title =element_text(size = 13,face="bold"), legend.text = element_text(size = 12,face="bold"))

### plot Figure7 in QC paper
cowplot::plot_grid(p.umap.BC.Dis, p.umap.IPS.Dis,labels = c("A", "B"),nrow=1,ncol=2)

###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
### %CV distribution based on pooled sample replicates
###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
temp1.OA <- OA.CV[[1]]
temp3.OA <- OA.CV[[2]]
temp5.OA <- OA.CV[[3]]
OA.CV.Frame <- data.frame(cbind(c(rep("Sample Repeats",7596),rep("Freeze Thaw",7596),rep("Reprocessed",7596)),c(temp1.OA,temp3.OA,temp5.OA)))
colnames(OA.CV.Frame) <- c("Conditions","CV")

p.OA.CV.cond <- ggplot(OA.CV.Frame) + stat_ecdf(aes(x=as.numeric(CV)*100,colour=Conditions),geom = "line") + labs(color="") +
  geom_hline(yintercept=0.8, linetype="dashed", color = "black") +
  geom_vline(xintercept=quantile(temp1.OA,0.8)*100, linetype="dashed", color = "blue") + 
  geom_vline(xintercept=quantile(temp3.OA,0.8)*100, linetype="dashed", color = "red") +
  geom_vline(xintercept=quantile(temp5.OA,0.8)*100, linetype="dashed", color = "green") + 
  scale_x_continuous(name="%CV", limits=c(0,80),breaks=seq(0,80,10)) + scale_y_continuous(name="%Cumulative Distribution for OA", limits=c(0,1),breaks=seq(0,1,0.1)) + theme_bw() + 
  theme(axis.title.x = element_text(size=11,face="bold"),legend.position = c(0.8,0.2),legend.title =element_text(size = 10,face="bold"), legend.text = element_text(size = 10,face="bold"), 
        axis.title.y =element_text(size=11,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

temp1.INJ <- INJ.CV[[1]]
temp3.INJ <- INJ.CV[[2]]
INJ.CV.Frame <- data.frame(cbind(c(rep("Sample Repeats",7596),rep("Freeze Thaw",7596)),c(temp1.INJ,temp3.INJ)))
colnames(INJ.CV.Frame) <- c("Conditions","CV")

p.INJ.CV.cond <- ggplot(INJ.CV.Frame) + stat_ecdf(aes(x=as.numeric(CV)*100,colour=Conditions),geom = "line") + labs(color="") +
  geom_hline(yintercept=0.8, linetype="dashed", color = "black") +
  geom_vline(xintercept=quantile(temp1.INJ,0.8)*100, linetype="dashed", color = "cyan4") + 
  geom_vline(xintercept=quantile(temp3.INJ,0.8)*100, linetype="dashed", color = "red") +
  scale_x_continuous(name="%CV", limits=c(0,80),breaks=seq(0,80,10)) + scale_y_continuous(name="%Cumulative Distribution for Injury", limits=c(0,1),breaks=seq(0,1,0.1)) + theme_bw() + 
  theme(axis.title.x = element_text(size=11,face="bold"),legend.position = c(0.8,0.2),legend.title =element_text(size = 10,face="bold"), legend.text = element_text(size = 10,face="bold"), 
        axis.title.y =element_text(size=11,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### plot FigureS2 in QC paper
cowplot::plot_grid(p.OA.CV.cond, p.INJ.CV.cond, labels = c("A","B"),nrow=1,ncol=2)

###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
### 18 paired spun/unspun samples visualisation on PCA 
###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
## extract out spun/unspun and pool samples
unspun_sams <- sort(grep("-UN",MySoma$STEpUpId,value=T))
spun_sams <- sort(grep("-SP",MySoma$STEpUpId,value=T))

MySoma_spun.before <- exprDat.batchDone[spun_sams,]
MySoma_unspun.before <- exprDat.batchDone[unspun_sams,]
MySoma_spun_unspun.before <- exprDat.batchDone[c(spun_sams,unspun_sams),]

MySoma_spun.after <- exprDat.batchDone2[spun_sams,]
MySoma_unspun.after <- exprDat.batchDone2[unspun_sams,]
MySoma_spun_unspun.after <- exprDat.batchDone2[c(spun_sams,unspun_sams),]

pca_svu.before <- prcomp(MySoma_spun_unspun.before)$x
pca_svu.after <- prcomp(MySoma_spun_unspun.after)$x

### movement on PCA 
p.spun.unspun.PC.before <- ggplot(data=data.frame(pca_svu.before)) + geom_point(aes(x=PC1,y=PC2,color=c(rep("Spun",18),rep("Unspun",18))),size=2) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific,limits=c(-4e5,1.3e6)) + 
  geom_line(aes(x=PC1,y=PC2,group=(rep(1:18,2)))) + labs(color="before IPS Adjustment") + xlab("PC1") + ylab("PC2") + guides(color = guide_legend(override.aes = list(size = 5))) + theme_bw() + 
  theme(axis.title.x = element_text(size=11,face="bold"),axis.title.y =element_text(size=11,face="bold"),
        legend.position="bottom",legend.title =element_text(size = 11,face="bold"),legend.text = element_text(size = 11,face="bold"))

p.spun.unspun.PC.after <- ggplot(data=data.frame(pca_svu.after)) + geom_point(aes(x=PC1,y=PC2,color=c(rep("Spun",18),rep("Unspun",18))),size=2) + scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) + 
  geom_line(aes(x=PC1,y=PC2,group=(rep(1:18,2)))) + labs(color="after IPS Adjustment") + xlab("PC1") + ylab("PC2") + guides(color = guide_legend(override.aes = list(size = 5))) + theme_bw() + 
  theme(axis.title.x = element_text(size=11,face="bold"),axis.title.y =element_text(size=11,face="bold"),
        legend.position="bottom",legend.title =element_text(size = 11,face="bold"), legend.text = element_text(size = 11,face="bold"))

##calculate correlation statistics
spun_unspun_cors.before <- sapply(1:ncol(MySoma_spun.before),function(i) cor(log(MySoma_spun.before[,i]),log(MySoma_unspun.before[,i])))
spun_unspun_cor_p.before <- sapply(which(colnames(MySoma_spun.before)=="seq.10000.28"):ncol(MySoma_spun.before),function(i) cor.test(log(MySoma_spun.before[,i]),log(MySoma_unspun.before[,i]))$p.val)
spun_unspun_cor_ci.before <- sapply(which(colnames(MySoma_spun.before)=="seq.10000.28"):ncol(MySoma_spun.before),function(i) cor.test(log(MySoma_spun.before[,i]),log(MySoma_unspun.before[,i]))$conf.int)

##calculate differential abundance statistics
spun_unspun_de_p.before <- sapply(which(colnames(MySoma_spun.before)=="seq.10000.28"):ncol(MySoma_spun.before),function(i) t.test(log(MySoma_spun.before[,i]),log(MySoma_unspun.before[,i]),paired=T)$p.val) 
spun_unspun_d.before <- sapply(which(colnames(MySoma_spun.before)=="seq.10000.28"):ncol(MySoma_spun.before),function(i) mean(log(MySoma_spun.before[,i]) - log(MySoma_unspun.before[,i]))/sd(log(c(MySoma_spun.before[,i],MySoma_unspun.before[,i]))))
names(spun_unspun_d.before) <- names(spun_unspun_de_p.before) <- colnames(MySoma_spun.before)[which(colnames(MySoma_spun.before)=="seq.10000.28"):ncol(MySoma_spun.before)]

cor_padj.before <- p.adjust(spun_unspun_cor_p.before,method="BH")
de_padj.before <- p.adjust(spun_unspun_de_p.before,method="BH")
padj.Frame.before <- data.frame(t(rbind(spun_unspun_d.before,spun_unspun_cors.before)))
colnames(padj.Frame.before) <- c("spun_unspun_d","spun_unspun_cors")

cols.before <- rep("NA",length(de_padj.before))
cols.before[cor_padj.before < 0.05 & de_padj.before > 0.05] <- "Correlated & the Same Mean"
cols.before[cor_padj.before > 0.05 & de_padj.before > 0.05] <- "Uncorrelated & the Same Mean"
cols.before[cor_padj.before > 0.05 & de_padj.before < 0.05] <- "Uncorrelated & Different Means"
cols.before[cor_padj.before < 0.05 & de_padj.before < 0.05] <- "Correlated & Different Means"

p.d.cor.before <- ggplot(padj.Frame.before) + geom_point(aes(x=spun_unspun_d,y=spun_unspun_cors,color=cols.before),size=0.3) + xlim(-1.5,1.25) + guides(color = guide_legend(override.aes = list(size = 5))) + 
  labs(color="before IPS Adjustment") + xlab("Spun vs Unspun Difference (Cohen's d)") + ylab("Spun and Unspun Correlation (rho)") + theme_bw() + 
  theme(axis.title.x = element_text(size=11,face="bold"),legend.title =element_text(size = 10,face="bold"), legend.text = element_text(size = 10,face="bold"),
        axis.title.y =element_text(size=11,face="bold"))

spun_unspun_cors.after <- sapply(1:ncol(MySoma_spun.after),function(i) cor(log(MySoma_spun.after[,i]),log(MySoma_unspun.after[,i])))
spun_unspun_cor_p.after<- sapply(which(colnames(MySoma_spun.after)=="seq.10000.28"):ncol(MySoma_spun.after),function(i) cor.test(log(MySoma_spun.after[,i]),log(MySoma_unspun.after[,i]))$p.val)
spun_unspun_cor_ci.after <- sapply(which(colnames(MySoma_spun.after)=="seq.10000.28"):ncol(MySoma_spun.after),function(i) cor.test(log(MySoma_spun.after[,i]),log(MySoma_unspun.after[,i]))$conf.int)

##calculate differential abundance statistics
spun_unspun_de_p.after <- sapply(which(colnames(MySoma_spun.after)=="seq.10000.28"):ncol(MySoma_spun.after),function(i) t.test(log(MySoma_spun.after[,i]),log(MySoma_unspun.after[,i]),paired=T)$p.val)
spun_unspun_d.after <- sapply(which(colnames(MySoma_spun.after)=="seq.10000.28"):ncol(MySoma_spun.after),function(i) mean(log(MySoma_spun.after[,i]) - log(MySoma_unspun.after[,i]))/sd(log(c(MySoma_spun.after[,i],MySoma_unspun.after[,i]))))
names(spun_unspun_d.after) <- names(spun_unspun_de_p.after) <- colnames(MySoma_spun.after)[which(colnames(MySoma_spun.after)=="seq.10000.28"):ncol(MySoma_spun.after)]

cor_padj.after <- p.adjust(spun_unspun_cor_p.after,method="BH")
de_padj.after <- p.adjust(spun_unspun_de_p.after,method="BH")
padj.Frame.after <- data.frame(t(rbind(spun_unspun_d.after,spun_unspun_cors.after)))
colnames(padj.Frame.after) <- c("spun_unspun_d","spun_unspun_cors")

cols.after <- rep("NA",length(de_padj.after))
cols.after[cor_padj.after < 0.05 & de_padj.after > 0.05] <- "Correlated & the Same Mean"
cols.after[cor_padj.after > 0.05 & de_padj.after > 0.05] <- "Uncorrelated & the Same Mean"
cols.after[cor_padj.after > 0.05 & de_padj.after < 0.05] <- "Uncorrelated & Different Means"
cols.after[cor_padj.after < 0.05 & de_padj.after < 0.05] <- "Correlated & Different Means"

p.d.cor.after <- ggplot(padj.Frame.after) + geom_point(aes(x=spun_unspun_d,y=spun_unspun_cors,color=cols.after),size=0.3) + xlim(-1.5,1.28) + guides(color = guide_legend(override.aes = list(size = 5))) + 
  labs(color="after IPS Regession") + xlab("Spun vs Unspun Difference (Cohen's d)") + ylab("Spun and Unspun Correlation (rho)") + theme_bw() + 
  theme(axis.title.x = element_text(size=11,face="bold"),legend.title =element_text(size = 10,face="bold"), legend.text = element_text(size = 10,face="bold"),
        axis.title.y =element_text(size=11,face="bold"))

### plot FigureS3 in QC paper
cowplot::plot_grid(p.spun.unspun.PC.before, p.d.cor.before, p.spun.unspun.PC.after, p.d.cor.after,
                   labels = c("A","C","B", "D"),nrow=2,ncol=2,rel_widths = c(2,3))

###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
### Blood staining
###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
CombinedFrame$sf_iknee_bloodstaining <- as.numeric(CombinedFrame$sf_iknee_bloodstaining)
BMarker = "HBA"
CommanName = "Hemoglobin (HBA)"
seqName <- "seq.17137.160"

keepID = which(!is.na(CombinedFrame$sf_iknee_bloodstaining))

CombinedHere <- cbind(CombinedFrame[,1:(which(colnames(CombinedFrame)=="seq.10000.28")-1)],exprDat.batchDone)

CorStr1.before = cor.test(log(exprDat.batchDone[keepID,seqName]),CombinedHere$sf_iknee_bloodstaining[keepID], method ="spearman")
blood.all.Frame <- CombinedHere[keepID,c(seqName,"sf_iknee_bloodstaining")]
p.blood.all.before <- ggplot(blood.all.Frame) + geom_boxplot(aes(x=as.character(sf_iknee_bloodstaining),y=log(get(seqName)),color=as.factor(sf_iknee_bloodstaining))) +
  scale_x_discrete(name = "", labels = c("no visual blood\ngrade1", "mildly blood stained\ngrade2\n(red tinge)","moderate blood staining\ngrade3\n(reduced light passage)","severe blood staining\ngrade4\n(opaque)")) +
  scale_color_discrete("Sample Size", labels=c(paste0("grade1: ",length(which(blood.all.Frame$sf_iknee_bloodstaining==1))),paste0("grade2: ",length(which(blood.all.Frame$sf_iknee_bloodstaining==2))),paste0("grade3: ",length(which(blood.all.Frame$sf_iknee_bloodstaining==3))),paste0("grade4: ",length(which(blood.all.Frame$sf_iknee_bloodstaining==4))))) + 
  ylab("log(HBA)\nbefore IPS Adjustment") + ggtitle(paste0("Correlation Coefficient ",signif(CorStr1.before$estimate,2))) + theme_bw() + 
  theme(plot.title=element_text(size = 9,face="bold",hjust=0.5),
        axis.text.x = element_text(size=7,face="bold"), axis.text.y = element_text(size=9,face="bold"),
        axis.title.x = element_text(size=9,face="bold"),axis.title.y =element_text(size=9,face="bold"),
        legend.position="none",legend.title =element_text(size = 9,face="bold"), legend.text = element_text(size = 8,face="bold"))

CorStrOA.before = cor.test(log(exprDat.batchDone[intersect(which(CombinedHere$sf_iknee_qc_group==0),keepID),seqName]),as.numeric(CombinedHere$sf_iknee_bloodstaining[intersect(which(CombinedHere$sf_iknee_qc_group==0),keepID)]),method ="spearman")
blood.all.OA.Frame <- CombinedHere[intersect(which(CombinedHere$sf_iknee_qc_group==0),keepID),c(seqName,"sf_iknee_bloodstaining")]
p.blood.OA.before <- ggplot(blood.all.OA.Frame) + geom_boxplot(aes(x=as.character(sf_iknee_bloodstaining),y=log(get(seqName)),color=as.factor(sf_iknee_bloodstaining))) +
  scale_x_discrete(name = "", labels = c("no visual blood\ngrade1", "mildly blood stained\ngrade2\n(red tinge)","moderate blood staining\ngrade3\n(reduced light passage)","severe blood staining\ngrade4\n(opaque)")) +
  scale_color_discrete("Sample Size", labels=c(paste0("grade1: ",length(which(blood.all.OA.Frame$sf_iknee_bloodstaining==1))),paste0("grade2: ",length(which(blood.all.OA.Frame$sf_iknee_bloodstaining==2))),paste0("grade3: ",length(which(blood.all.OA.Frame$sf_iknee_bloodstaining==3))),paste0("grade4: ",length(which(blood.all.OA.Frame$sf_iknee_bloodstaining==4))))) + 
  ylab("log(HBA) OA Group\nbefore IPS Adjustment") + ggtitle(paste0("Correlation Coefficient ",signif(CorStrOA.before$estimate,2))) + theme_bw() + 
  theme(plot.title=element_text(size = 9,face="bold",hjust=0.5),
        axis.text.x = element_text(size=7,face="bold"), axis.text.y = element_text(size=9,face="bold"),
        axis.title.x = element_text(size=9,face="bold"),axis.title.y =element_text(size=9,face="bold"),
        legend.position="none",legend.title =element_text(size = 9,face="bold"), legend.text = element_text(size = 8,face="bold"))

CorStrINJ.before = cor.test(log(exprDat.batchDone[intersect(which(CombinedHere$sf_iknee_qc_group==1),keepID),seqName]),as.numeric(CombinedHere$sf_iknee_bloodstaining[intersect(which(CombinedHere$sf_iknee_qc_group==1),keepID)]),method ="spearman")
blood.all.INJ.Frame <- CombinedHere[intersect(which(CombinedHere$sf_iknee_qc_group==1),keepID),c(seqName,"sf_iknee_bloodstaining")]
p.blood.INJ.before <- ggplot(blood.all.INJ.Frame) + geom_boxplot(aes(x=as.character(sf_iknee_bloodstaining),y=log(get(seqName)),color=as.factor(sf_iknee_bloodstaining))) +
  scale_x_discrete(name = "", labels = c("no visual blood\ngrade1", "mildly blood stained\ngrade2\n(red tinge)","moderate blood staining\ngrade3\n(reduced light passage)","severe blood staining\ngrade4\n(opaque)")) +
  scale_color_discrete("Sample Size", labels=c(paste0("grade1: ",length(which(blood.all.INJ.Frame$sf_iknee_bloodstaining==1))),paste0("grade2: ",length(which(blood.all.INJ.Frame$sf_iknee_bloodstaining==2))),paste0("grade3: ",length(which(blood.all.INJ.Frame$sf_iknee_bloodstaining==3))),paste0("grade4: ",length(which(blood.all.INJ.Frame$sf_iknee_bloodstaining==4))))) + 
  ylab("log(HBA) Injury Group\nbefore IPS Adjustment") + ggtitle(paste0("Correlation Coefficient ",signif(CorStrINJ.before$estimate,2))) + theme_bw() + 
  theme(plot.title=element_text(size = 9,face="bold",hjust=0.5),
        axis.text.x = element_text(size=7,face="bold"), axis.text.y = element_text(size=9,face="bold"),
        axis.title.x = element_text(size=9,face="bold"),axis.title.y =element_text(size=9,face="bold"),
        legend.position="none",legend.title =element_text(size = 9,face="bold"), legend.text = element_text(size = 8,face="bold"))

CombinedHere2 <- cbind(CombinedFrame[,1:(which(colnames(CombinedFrame)=="seq.10000.28")-1)],exprDat.batchDone2)

CorStr1.after = cor.test(log(exprDat.batchDone2[keepID,seqName]),CombinedHere2$sf_iknee_bloodstaining[keepID], method ="spearman")
blood.all.Frame <- CombinedHere2[keepID,c(seqName,"sf_iknee_bloodstaining")]
p.blood.all.after <- ggplot(blood.all.Frame) + geom_boxplot(aes(x=as.character(sf_iknee_bloodstaining),y=log(get(seqName)),color=as.factor(sf_iknee_bloodstaining))) +
  scale_x_discrete(name = "", labels = c("no visual blood\ngrade1", "mildly blood stained\ngrade2\n(red tinge)","moderate blood staining\ngrade3\n(reduced light passage)","severe blood staining\ngrade4\n(opaque)")) +
  scale_color_discrete("Sample Size", labels=c(paste0("grade1: ",length(which(blood.all.Frame$sf_iknee_bloodstaining==1))),paste0("grade2: ",length(which(blood.all.Frame$sf_iknee_bloodstaining==2))),paste0("grade3: ",length(which(blood.all.Frame$sf_iknee_bloodstaining==3))),paste0("grade4: ",length(which(blood.all.Frame$sf_iknee_bloodstaining==4))))) + 
  ylab("log(HBA)\nafter IPS Adjustment") + ggtitle(paste0("Correlation Coefficient ",signif(CorStr1.after$estimate,2))) + theme_bw() + 
  theme(plot.title=element_text(size = 9,face="bold",hjust=0.5),
        axis.text.x = element_text(size=7,face="bold"), axis.text.y = element_text(size=9,face="bold"),
        axis.title.x = element_text(size=9,face="bold"),axis.title.y =element_text(size=9,face="bold"),
        legend.position="right",legend.title =element_text(size = 9,face="bold"), legend.text = element_text(size = 8,face="bold"))

CorStrOA.after = cor.test(log(exprDat.batchDone2[intersect(which(CombinedHere2$sf_iknee_qc_group==0),keepID),seqName]),as.numeric(CombinedHere2$sf_iknee_bloodstaining[intersect(which(CombinedHere2$sf_iknee_qc_group==0),keepID)]),method ="spearman")
blood.all.OA.Frame <- CombinedHere2[intersect(which(CombinedHere2$sf_iknee_qc_group==0),keepID),c(seqName,"sf_iknee_bloodstaining")]
p.blood.OA.after <- ggplot(blood.all.OA.Frame) + ggtitle(paste0("Correlation Coefficient ",signif(CorStrOA.after$estimate,2))) + geom_boxplot(aes(x=as.character(sf_iknee_bloodstaining),y=log(get(seqName)),color=as.factor(sf_iknee_bloodstaining))) +
  scale_x_discrete(name = "", labels = c("no visual blood\ngrade1", "mildly blood stained\ngrade2\n(red tinge)","moderate blood staining\ngrade3\n(reduced light passage)","severe blood staining\ngrade4\n(opaque)")) +
  scale_color_discrete("Sample Size", labels=c(paste0("grade1: ",length(which(blood.all.OA.Frame$sf_iknee_bloodstaining==1))),paste0("grade2: ",length(which(blood.all.OA.Frame$sf_iknee_bloodstaining==2))),paste0("grade3: ",length(which(blood.all.OA.Frame$sf_iknee_bloodstaining==3))),paste0("grade4: ",length(which(blood.all.OA.Frame$sf_iknee_bloodstaining==4))))) + 
  ylab("log(HBA) OA Group\nafter IPS Adjustment") + theme_bw() + 
  theme(plot.title=element_text(size = 9,face="bold",hjust=0.5),
        axis.text.x = element_text(size=7,face="bold"), axis.text.y = element_text(size=9,face="bold"),
        axis.title.x = element_text(size=9,face="bold"),axis.title.y =element_text(size=9,face="bold"),
        legend.position="right",legend.title =element_text(size = 9,face="bold"), legend.text = element_text(size = 8,face="bold"))

CorStrINJ.after = cor.test(log(exprDat.batchDone2[intersect(which(CombinedHere2$sf_iknee_qc_group==1),keepID),seqName]),as.numeric(CombinedHere2$sf_iknee_bloodstaining[intersect(which(CombinedHere2$sf_iknee_qc_group==1),keepID)]),method ="spearman")
blood.all.INJ.Frame <- CombinedHere2[intersect(which(CombinedHere2$sf_iknee_qc_group==1),keepID),c(seqName,"sf_iknee_bloodstaining")]
p.blood.INJ.after <- ggplot(blood.all.INJ.Frame) + geom_boxplot(aes(x=as.character(sf_iknee_bloodstaining),y=log(get(seqName)),color=as.factor(sf_iknee_bloodstaining))) +
  scale_x_discrete(name = "", labels = c("no visual blood\ngrade1", "mildly blood stained\ngrade2\n(red tinge)","moderate blood staining\ngrade3\n(reduced light passage)","severe blood staining\ngrade4\n(opaque)")) +
  scale_color_discrete("Sample Size", labels=c(paste0("grade1: ",length(which(blood.all.INJ.Frame$sf_iknee_bloodstaining==1))),paste0("grade2: ",length(which(blood.all.INJ.Frame$sf_iknee_bloodstaining==2))),paste0("grade3: ",length(which(blood.all.INJ.Frame$sf_iknee_bloodstaining==3))),paste0("grade4: ",length(which(blood.all.INJ.Frame$sf_iknee_bloodstaining==4))))) + 
  ylab("log(HBA) Injury Group\nafter IPS Adjustment") + ggtitle(paste0("Correlation Coefficient ",signif(CorStrINJ.after$estimate,2))) + theme_bw() + 
  theme(plot.title=element_text(size = 9,face="bold",hjust=0.5),
        axis.text.x = element_text(size=7,face="bold"), axis.text.y = element_text(size=9,face="bold"),
        axis.title.x = element_text(size=9,face="bold"),axis.title.y =element_text(size=9,face="bold"),
        legend.position="right",legend.title =element_text(size = 9,face="bold"), legend.text = element_text(size = 8,face="bold"))

p.blood1  <- ggpubr::ggarrange(p.blood.all.before,p.blood.all.after,labels = c("A", "B"),common.legend=TRUE,legend="right")
p.blood2  <- ggpubr::ggarrange(p.blood.OA.before,p.blood.OA.after,labels = c("C", "D"),common.legend=TRUE,legend="right")
p.blood3  <- ggpubr::ggarrange(p.blood.INJ.before,p.blood.INJ.after,labels = c("E", "F"),common.legend=TRUE,legend="right")

### plot FigureS4 in QC paper
cowplot::plot_grid(p.blood1,p.blood2,p.blood3,nrow=3,ncol=1,align = "v")

###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
### %CV and R2 distribution based on OA and injury pooled sample replicates after QC
###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
exprDat <- MySoma2022[grep("STEP",MySoma2022$SampleId),which(colnames(MySoma2022)=="seq.10000.28"):ncol(MySoma2022)]
rownames(exprDat) <- MySoma2022[grep("STEP",MySoma2022$SampleId),"SampleId"]

tempOA.cv.final <- CalibratorCheck2(MySoma2022,"OA",exprDat,"OA")
tempINJ.cv.final <- CalibratorCheck2(MySoma2022,"INJ",exprDat,"Injury")

OA.INJ.CV.Frame <- data.frame(cbind(c(rep("OA Sample Repeats",7596),rep("Injury Sample Repeats",7596)),c(tempOA.cv.final,tempINJ.cv.final)))
colnames(OA.INJ.CV.Frame) <- c("Conditions","CV")

p.OA.INJ.CV <- ggplot(OA.INJ.CV.Frame) + stat_ecdf(aes(x=as.numeric(CV)*100,colour=Conditions),geom = "line") + labs(color="") +
  geom_hline(yintercept=0.8, linetype="dashed", color = "black") +
  geom_vline(xintercept=quantile(tempOA.cv.final,0.8)*100, linetype="dashed", color = "blue") + 
  geom_vline(xintercept=quantile(tempINJ.cv.final,0.8)*100, linetype="dashed", color = "red") +
  scale_x_continuous(name="%CV", limits=c(0,80),breaks=seq(0,80,10)) + scale_y_continuous(name="%Cumulative Distribution", limits=c(0,1),breaks=seq(0,1,0.1)) + theme_bw() + 
  theme(axis.title.x = element_text(size=11,face="bold"),legend.position = c(0.75,0.2),legend.title =element_text(size = 10,face="bold"), legend.text = element_text(size = 10,face="bold"),
        axis.title.y =element_text(size=11,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  annotate(geom="text",x=13, y=0.82, label=as.character(signif(quantile(tempOA.cv.final,seq(0,1,length.out=100))[80],4)),color="blue") +
  annotate(geom="text",x=21, y=0.78, label=as.character(signif(quantile(tempINJ.cv.final,seq(0,1,length.out=100))[80],4)),color="red")

R2.OA.final <- VarExp2(MySoma2022,"OA",exprDat,"")
R2.INJ.final <- VarExp2(MySoma2022,"INJ",exprDat,"")

OA.INJ.R2.Frame <- data.frame(cbind(c(rep("OA Sample Repeats",7596),rep("Injury Sample Repeats",7596)),c(R2.OA.final,R2.INJ.final)))
colnames(OA.INJ.R2.Frame) <- c("Conditions","R2")

p.OA.INJ.R2 <- ggplot(OA.INJ.R2.Frame) + stat_ecdf(aes(x=as.numeric(R2),colour=Conditions),geom = "line") + labs(color="") + 
  geom_hline(yintercept=0.8, linetype="dashed", color = "black") +
  geom_vline(xintercept=quantile(R2.OA.final,0.2), linetype="dashed", color = "blue") +
  geom_vline(xintercept=quantile(R2.INJ.final,0.2), linetype="dashed", color = "red") +
  scale_x_reverse(name=bquote(bold("Mean R")^bold("2")),limits=c(1,0),breaks=seq(1,0,-0.1)) + scale_y_continuous(name="%Cumulative Distribution", limits=c(0,1),breaks=seq(0,1,0.1)) + theme_bw() + 
  theme(axis.title.x = element_text(size=11,face="bold"),legend.position = c(0.75,0.2),legend.title =element_text(size = 10,face="bold"), legend.text = element_text(size = 10,face="bold"), 
        axis.title.y =element_text(size=11,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate(geom="text",x=0.92, y=0.82, label=as.character(signif(quantile(R2.OA.final,1 - seq(0,1,length.out=100))[80],4)),color="blue") +
  annotate(geom="text",x=0.8, y=0.78, label=as.character(signif(quantile(R2.INJ.final,1 - seq(0,1,length.out=100))[80],4)),color="red")

### plot in FigureS5 in QC paper 
cowplot::plot_grid(p.OA.INJ.CV, p.OA.INJ.R2, labels = c("A","B"),nrow=1,ncol=2)

###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
### Bimodal signal visualisation on paired top 5 PCs 
###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
p.PCA.BimodalBefore <- ggpairs(data.frame(pcDat.standardised[,1:5]), columns=1:5, aes(color= as.factor(bimodalLabel)), 
                               diag=list(continuous=wrap("densityDiag",alpha=0.4)),
                               lower=list(continuous = wrap("points",alpha=0.9,size=0.5)),
                               upper = list(continuous = "blank"),
                               legend = c(3,1)) + labs(color="before Batch Correction") + theme_bw() +
  theme(axis.title.x = element_text(size=11,face="bold"),axis.text.x = element_text(size=10,face="bold"),
        axis.title.y =element_text(size=11,face="bold"), axis.text.y = element_text(size=10,face="bold"),
        legend.title =element_text(size = 11,face="bold"), legend.text = element_text(size = 10,face="bold"))



p.PCA.BimodalAfter <- ggpairs(data.frame(pcDat.batchDone[,1:5]), columns=1:5, aes(color= as.factor(bimodalLabel)), 
                              diag=list(continuous=wrap("densityDiag",alpha=0.4)),
                              lower=list(continuous = wrap("points",alpha=0.9,size=0.5)),
                              upper = list(continuous = "blank"),
                              legend= c(3,1)) + labs(color="after Batch Correction") + theme_bw() + 
  theme(axis.title.x = element_text(size=11,face="bold"),axis.text.x = element_text(size=10,face="bold"),
        axis.title.y =element_text(size=11,face="bold"), axis.text.y = element_text(size=10,face="bold"),
        legend.title =element_text(size = 11,face="bold"), legend.text = element_text(size = 10,face="bold"))

### plot FigureS7 in QC paper
p.PCA.BimodalBefore
p.PCA.BimodalAfter

###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
### Filter summary for batch corrected non-IPS adjusted and IPS adjusted data
###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
### Summary of filters for samples and proteins, as shown in TableS4 in QC paper
getRemoveTable(NonHuman,RemoveS1,RemoveS2,RemoveS3,RemoveS4,RemovePro1,RemovePro2,bimodalP1,ConfounderTable1.batchDone[,c(7,10)],exprDat.batchDone)[[1]]
getRemoveTable(NonHuman,RemoveS1,RemoveS2,RemoveS3,RemoveS4,RemovePro1,RemovePro2,bimodalP2,ConfounderTable1.batchDone2[,c(7,10)],exprDat.batchDone2)[[1]]

###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
### Summary of associations between technical confounders and top 10 PCs
###_______________________________________________________________________________________________________________________________________________________
###_______________________________________________________________________________________________________________________________________________________
### for standardised data
View(signif(ConfounderTable2.Standerdised[1:10,],3))
ConfounderTable2.Standerdised.padj <- sapply(1:ncol(ConfounderTable2.Standerdised),function(x){p.adjust(ConfounderTable2.Standerdised[,x],method="BH")})
signif(eig.val.standardised.all[1:10],2) ### variation expained per top 10 PC

### for batch corrected, non-IPS adjusted, non-filtered data
View(signif(ConfounderTable2.batchCorrected[1:10,],3))
ConfounderTable2.batchCorrected.padj <- sapply(1:ncol(ConfounderTable2.batchCorrected),function(x){p.adjust(ConfounderTable2.batchCorrected[,x],method="BH")})
signif(eig.val.batchDone[1:10],2) ### variation expained per top 10 PC

### for batch corrected, IPS adjusted, non-filtered data
View(signif(ConfounderTable2.IPSreg[1:10,],3))
ConfounderTable2.IPSreg.padj <- sapply(1:ncol(ConfounderTable2.IPSreg),function(x){p.adjust(ConfounderTable2.IPSreg[,x],method="BH")})
signif(eig.val.batchDone2[1:10],2) ### variation expained per top 10 PC

### for batch corrected, non-IPS adjusted, filtered data
View(signif(ConfounderTable2.batchCorrected.Filtered[1:10,],3)) 
ConfounderTable2.batchCorrected.Filtered.padj <- sapply(1:ncol(ConfounderTable2.batchCorrected.Filtered),function(x){p.adjust(ConfounderTable2.batchCorrected.Filtered[,x],method="BH")})
signif(eig.val.batchCorrected.Filtered[1:10],2) ### variation expained per top 10 PC

### for batch corrected, IPS adjusted, filtered data
View(signif(ConfounderTable2.IPSreg.Filtered[1:10,],3))
ConfounderTable2.IPSreg.Filtered.padj <- sapply(1:ncol(ConfounderTable2.IPSreg.Filtered),function(x){p.adjust(ConfounderTable2.IPSreg.Filtered[,x],method="BH")})
signif(eig.val.IPSreg.Filtered[1:10],2) ### variation expained per top 10 PC

### TableS7 in QC paper summaries the above five data frames