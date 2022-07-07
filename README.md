# DRIM
Deconvolution followed by Region-growing,Interpolation and iterative-Mapping  
https://cell.ownbox.cn/
## Installation
DRIM is available from GitHub with:
```
# install.packages("devtools")
devtools::install_github("zenglab-regeneration/DRIM")
```
**Request**
If you have already installed conda environment, please skip this step.
-  **install **[anaconda](https://www.anaconda.com/ "anaconda")
- **Create a python conda environment**
**example: conda create -n test python=3.9**
```
library(DRIM)
help(package="DRIM")
```
Set the path of python conda
```
DRIM::env_python_set("D:/anaconda/envs/py")
```

Test if environment dependencies are satisfied
```
DRIM::env_test()
```

**sc_rds**_:*Single-cell transcriptome data*  
**st_rds**_:Spatial transcriptome data    
**plot_data**_:Cell type file divided by deconvolution    
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
DRIM::data_deal(sc_exp_data = sc_rds_@assays$SCT@data,st_exp_data = st_rds_$timing_0h@assays$Spatial@data,sc_celltype_data = sc_rds_@meta.data,
                loc_data = st_rds_$timing_0h@images$slice1@coordinates,plot_data = plot_data_)
```
**resolution**：magnification  
**colname**：selected data column  
```
result <- DRIM::planarian_main(resolution = 4,colname="final_celltype")
```

```
