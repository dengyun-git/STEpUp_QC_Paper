---
title: "STEpUP OA: Script used to process 'minimal' dataset required to generate figures in our Quality-Control (QC) Manuscript"
output:
  html_document:
    keep_md: yes
    toc: yes
    toc_depth: 3
    fig_width: 12
    fig_height: 7
  pdf_document:
    toc: yes
    toc_depth: '3'
---



## Purpose of this directory:

The full STEpUP OA dataset may be made available by application to the Data Access and Publication Group of STEpUP OA (stepupoa@kennedy.ox.ac.uk) once the primary analysis manuscript is published, in accordance with what is stipulated in our Consortium Agreement. The minimal datasets necessary for replicating figures along with the required R code are provided here. 

## Directory structure:

```
##                                      levelName
## 1 minimal datasets                            
## 2  ¦--PC1 Driver - Standardisation            
## 3  ¦--PC1 Driver - Intracellular Protein Score
## 4  ¦--PC2 Driver - Bimodal Signal             
## 5  ¦--Compare to Immunoassay                  
## 6  ¦--Disease Group After Filtering           
## 7  ¦--vignette for replication.Rmd            
## 8  °--vignette-for-replication.html
```

## The minimal datasets to replicate figures in the paper:
Minimal data to replicate main figures in the paper are stored in the folders: "PC1 Driver - Standardisation", "PC1 Driver - Intracellular Protein Score", "PC2 Driver - Bimodal Signal", "Compare to Immunoassay", and "Disease Group After Filtering".

## R code required to process the minimal data to replicate figures in the paper
A guide for users to utilize R for processing the minimal data to replicate figures in the paper is presented in the 'vignette-for-replication.Rmd' file (to be opened in the R editor) and 'vignette-for-replication.html' file (to be opened directly for convenient viewing).


