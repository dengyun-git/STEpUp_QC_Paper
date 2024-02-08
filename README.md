<img src="STEPUPlogo.png" width="30%" />

Date: Friday 30th June 2023

Project: STEpUP OA

Author: Dr Yun Deng     | Email:  yun.deng @ kennedy.ox.ac.uk\
Author: Dr.Luke jostins | Email:  luke.jostins @ kennedy.ox.ac.uk\
Author: Dr.Thomas Perry | Email:  thomas.perry @ ndorms.ox.ac.uk

Github repository access: https://github.com/dengyun-git/STEpUp_QC_Paper (administrator: Y.Deng)

Summary: this repository holds scripts & data required to generate figures/plots/tables for our 'QC Manuscript'

### STEpUp_QC_Paper repository - overview
This repository comprises scripts and data required to process SomaLogic proteomic data used in the STEpUP OA project, and to generate figures, plots and results for our 'QC Manuscript' entitled:
'Methodological development of molecular endotype discovery from synovial fluid of individuals with knee osteoarthritis: the STEpUP OA Consortium'

### Repository structure


```
##                                                               levelName
## 1 STEpUP_QC_Paper-main                                                 
## 2  ¦--Scripts                                                          
## 3  ¦   °--4 R scripts used to prepare and generate figures plots tables
## 4  ¦--Resource File                                                    
## 5  ¦   °--11 files used by the various R scripts                       
## 6  ¦--Robj.Paper                                                       
## 7  °--minimal datasets
```

### Repository scripts and data (as stored in the GitHub folder: )
* "Scripts" folder includes x4 R-scripts which processed the proteomic data and clinical data are described below;
  1) QC.paper.202306.R --------> main script used to transform Somalogic synovial fluid (SF) proteomic data and to perform a series of data quality assessments (i.e. our 'QC pipeline').
   2) QC.plot.202306.R ---------> script used to make plots and tables as shown in our 'QC manuscript'. 
   3) QCnorm.Paper.202306.R ----> functions used to perform each of our standardization steps. This function is called by QC.paper.202306.R 
   4) QCassess.Paper.202306.R --> functions used to perform quality evaluation, and to generate plots and tables. This function is called by QC.paper.202306.R.\
These scripts are stored in the Github folder called 'Scripts'.

* Full dataset required by R code in "Scripts" folder includes x4 .xlsx/csv files, 3x .txt files and 4x .tsv/.adat files which are required by the scripts: 
   1) adat_to_redcap_SIN_map.csv 
   2) discovery_QApheno_1_120522.csv 
   3) discovery_QApheno_2_270522.xlsx 
   4) Masterlist.xlsx 
   5) Cytoplasm.txt 
   6) Endomembrane.txt 
   7) Nucleus.txt 
   8) SS-228545_v4.1.20220224.adat 
   9) SS-228545_v4.1.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.20220224.adat 
   10) subcellular_location.tsv (online resources downloaded from: https://www.proteinatlas.org/humanproteome/subcellular)
   11) PanglaoDB_markers_27_Mar_2020.tsv (online resources downloaded from: https://panglaodb.se) \
Currently only four files without privacy information about patients are stored in "Resource File". The full STEpUP OA dataset may be made available by application to the Data Access and Publication Group of STEpUP OA (stepupoa@kennedy.ox.ac.uk). 

However for the purpose of users to replicate results in the paper, we prepared the folder "minimal datasets", which contains the 'minimal' dataset and R code to generate figures in our Quality-Control (QC) Manuscript.
   
### General notes and comments for use
* Once the user has downloaded the 'Resource File' folder and the 'Scripts' folder from GitHub to their 'Downloads' folder on their local machine, the following user input is required:
  1) Script: 'QC.paper.202306.R'
   --> on line 15 of the code, the user will need to change the path name: userPath <-       "/Users/ydeng/Documents/QCpaper.Code/" ---> CHANGE TO user local path <--- e.g. "C:/Users/tperry/Downloads/"

  2) Script: 'QC.plot.202306.R'
   --> on line 14 of the code, the user will need to change the path name: userPath <- "/Users/ydeng/Documents/QCpaper.Code/" --->CHANGE TO user local path <--- e.g. "C:/Users/tperry/Downloads/"
   --> the user will then need to save this script

* The following scripts need to be altered/run in the following order:
   Order: 1) QC.paper.202306.R first & 2) QC.plot.202306.R
   
* The script called 'QC.paper.202306' will generate a new folder called 'Robj.Paper' in the users personal Downloads folder (if the user is working out of this folder location) and will store R objects for plots and tables.  

* To use the minimal dataset to replicate figures in the paper please refer to the detailed guide described in "vignette-for-replication.html" under the "minimal datasets" folder.

### Other
* The user will need to install the R packages:
readxl, SomaDataIO, nlme, sva, factoextra, umap, ggplot2, GGally, mixtools, data.table, openxlsx, writexl, ggpubr, gdata cowplot, ggpubr, ggforce, scales.\
  An example to install package "limma" with commond as follows:\
  if (!requireNamespace("BiocManager", quietly = TRUE))\
  install.packages("BiocManager")\
  BiocManager::install("limma")

