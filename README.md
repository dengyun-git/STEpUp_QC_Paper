---
output: 
  html_document: 
    keep_md: yes
---
<img src="STEPUPlogo.png" width="30%" />

Breif direcotry structure of the codes and input/output data shows as follows


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
Brief introduction of each R script and input/output folder.

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

