Description:
===
This document provides the analysis code for all the data and images mentioned in this article

Requirements:
===
R 4.0.1

Packages needed to be installed: MetaCycle,pheatmap,VennDiagram,ggplot2,cowplot,ggplotify,dplyr,forcats

Usage:
===
Cycling gene calculation: 
---
All species expect yeast:

meta2d(infile=input_file, filestyle=“csv”, outdir=output_file,timepoints=“Line1”,outRawData=TRUE,minper=20,maxper=28, cycMethod = c(“ARS”, “JTK”, “LS”))
Yeast high glucose:

meta2d(infile=input_file, filestyle=“csv”, outdir=output_file,timepoints=“Line1”,outRawData=TRUE,minper=4,maxper=14, cycMethod = c(“ARS”, “JTK”, “LS”),ARSdefaultPer=120,ARSmle=“nomle”)
Yeast low glucose:

 meta2d(infile=input_file, filestyle=“csv”, outdir=output_file,timepoints=“Line1”,outRawData=TRUE,minper=5,maxper=15, cycMethod = c(“ARS”, “JTK”, “LS”),ARSdefaultPer=360,ARSmle=“nomle”)
 
 Data contained in this paper:
 ---
 All data in this paper can be found at 'ContainedData' file
 
 Code in this paper:
 ---
 All code in this paper can be found in all .md files
