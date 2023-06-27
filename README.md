---
<img src="STEPUPlogo.png" width="30%" />

### About STEpUp_QC_Paper repository
R scripts in this repository constructed processing and quality control pipeline for proteomic data of Synovial Fluid measured on SOMAscan platform. Outputs from the codes in this repository contributed to the paper "Methodological development of molecular endotype discovery from synovial fluid of individuals with knee osteoarthritis: the STEpUP OA Consortium".

### Breif direcotry structure of the input/output data referred in the codes shows as follows


```
##                           levelName
## 1  STEpUp_QC_Paper                 
## 2   ¦--STEPUP_DAG_rel002_discovery1
## 3   ¦   ¦--online resource         
## 4   ¦   ¦--somascan                
## 5   ¦   °--clinical                
## 6   ¦--Robj.Paper                  
## 7   ¦--QC.paper.202306.R           
## 8   ¦--QCnorm.Paper.202306.R       
## 9   ¦--QCassess.Paper.202306.R     
## 10  °--QC.plot.202306.R
```

### About input/output folder
* online resource: subcellular locations of proteins: https://www.proteinatlas.org/humanproteome/subcellular
                   cell markers: https://panglaodb.se/ 
* somascan: Adat file including protein expression profiles from SOMAscan platform
* clinical: clinical matching variables used for quality control 
* Robj.Paper: stored R objects, for convenient referencing to make plots and tables  

### About each R script
* QC.paper.202306.R: major code to transform data and perform series of data quality assessments.
* QCnorm.Paper.202306.R: functions of standardization steps, called by QC.paper.202306.R 
* QCassess.Paper.202306.R: functions to perform quality evaluation, and use plots and tables for visualisation, called by QC.paper.202306.R.
* QC.plot.202306.R: summurised code to make plots and tables for QC paper. Need to run QC.paper.202306.R first.

