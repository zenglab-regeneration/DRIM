# DRIM
DRIM:Deconvolution followed by Region-growing,Interpolation and iterative-Mapping  
website url:https://cell.ownbox.cn/

## Overview
## Installation
### Setup
DRIM is available from GitHub with:

```
#If you don't have devtools installed, please install it first

install.packages("devtools")

devtools::install_github("zenglab-regeneration/DRIM")

```

### Depend

Some programs of our project require a python environment, so if you don't have a python environment, please follow the steps below.  
* install [anaconda](https://www.anaconda.com/ "anaconda")
* Create a python3.9 environment  

When you have anaconda installed, you need to create a python conda.
```
conda create -n testconda python = 3.9
```

## Examples from paper
### Dataset 
- Single-cell transcriptome gene expression matrix
- Spatial transcriptome gene expression matrix
- Single cell type data
### Environment settings


```
library('DRIM')
help(package = 'DRIM')

#set the py conda
env_python_set("D:/anaconda/envs/testconda")

#Check the dependent environment for the program to run, and automatically install the missing python package
env_test()
```
### Datadeal
**sc_rds**_:*Single-cell transcriptome data*  
**st_rds**_:Spatial transcriptome data    
**plot_data**_:Cell type file divided by deconvolution  
If you don't have plot data, don't worry, you can use **data_plot** to generate it.

```
sc_rds_<-readRDS("D:/code/data/sct_data_1.5d.rds")
st_rds_<-readRDS("D:/code/data/spatial_obj.rds")
plot_data_=read.csv("D:/code/data/deconvolution.csv", row.names = 1,header = TRUE)
```
**sc_exp_data**：Single-cell transcriptome gene expression matrix  
**st_exp_data**：Spatial transcriptome gene expression matrix  
**sc_celltype_data**：meta.data  
**loc_data **：  
**plot_data**：Cell type file divided by deconvolution  
```
DRIM::data_deal(sc_exp_data = sc_rds_@assays$SCT@data,  
                st_exp_data = st_rds_$timing_0h@assays$Spatial@data,
                sc_celltype_data = sc_rds_@meta.data,
                loc_data = st_rds_$timing_0h@images$slice1@coordinates,
                plot_data = plot_data_)
```
### Run
**resolution**：magnification  
**colname**：selected data column  
```
result <- DRIM::planarian_main(resolution = 4,colname="final_celltype")
```
return is a S4 object , you could use result@Amplification.

