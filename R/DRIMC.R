#' @aliases DRIM
#' @author sh hay
packages_path <- function(){
  packages_dir <- system.file(package = "DRIM")
  return(packages_dir)
}
#' @title pythonEnvSet
#' @description
#' Set the path of the conda environment you use
#' Many python programs are called in our program,
#' so you need to set the path of the python interpreter,
#' and you can set it directly to the conda environment here.
#' @param py_path path of the conda environment
#' @import reticulate
#' @export
pythonEnvSet <- function(py_path){
  # library(reticulate)
  use_condaenv(py_path)
  #use_python(py_path)
}
#' @title testEnv
#' Detect environment dependencies of python
#' @description
#' You can use this function to detect if a package is missing from a dependent python environment.
#' @return bool TRUE or FALSE
#' @import reticulate
#' @export
testEnv <- function(){
  library(reticulate)
  package_flag <- TRUE
  package_detect <- c('anndata',
                    'collections',
                    'copy',
                    'datetime',
                    'math',
                    'matplotlib',
                    'numba',
                    'numpy==1.20.3',
                    'operator',
                    'pandas',
                    'pathlib',
                    'random',
                    'scanpy',
                    'seaborn',
                    'sklearn',
                    'time',
                    'rich',
                    'plotly',
                    'kaleido==0.1.0')
  for(package_name in package_detect){
    if(!py_module_available(package_name)){
      py_install(package_name, pip = T)
    }
  }
  for(package_name in package_detect){
    if(!py_module_available(package_name)){
      error_message <- paste(package_name,"is not ready",sep=" ")
      print(error_message)
    }
  }
  return(package_flag)
}
#' @title dataProcessing
#' data processing
#' @description
#' Preprocess the data to get the data we need later.
#' @param sc_exp_data Single Cell Transcriptome Expression Matrix
#' @param st_exp_data Spatial transcriptome expression matrix
#' @param sc_celltype_data Cell type mete data
#' @param loc_matrix Spatial transcriptome coordinates information
#' @param deconvolution_matrix Spatial spot deconvolution results
#' @import data.table
#' @export
dataProcessing <- function(sc_exp_data,st_exp_data,sc_celltype_data,loc_matrix,deconvolution_matrix){
  library(data.table)
  data_path <- paste(packages_path(),'/data',sep = "")
  if(!dir.exists(data_path)){
    dir.create(data_path)
  }
  fwrite(as.data.frame(deconvolution_matrix),file = paste(data_path,"/deconvolution.csv",sep = ""))
  fwrite(as.data.frame(sc_exp_data),file = paste(data_path,"/sc_exp.csv",sep = ""),row.names = TRUE)
  fwrite(as.data.frame(st_exp_data),file = paste(data_path,"/st_exp.csv",sep = ""),row.names = TRUE)
  fwrite(as.data.frame(loc_matrix),file = paste(data_path,'/st_loc.csv',sep = ""))
  fwrite(as.data.frame(sc_celltype_data),file = paste(data_path,'/sc_celltype.csv',sep = ""))
}

#' @title setParameters
#' Hyperparameter setting
#' @description
#' Set Resolution and Cell Columns,resolution is the multiple of program amplification,
#' Cell Columns is the specified column name
#' @param resolution program magnification
#' @param colname selected column name
setParameters <- function(resolution=4,thread=7,colname){
  if(thread==7){
    parameter_settings_csv <- c(resolution,colname)
  }
  else {
     parameter_settings_csv <-c (resolution,colname,thread)
  }
  dir=packages_path()
  parameter_settings_path = paste(dir,"/data/parameter_settings.csv",sep="")
  write.table(parameter_settings_csv,file = parameter_settings_path,row.names = FALSE,col.names = 'parameter')
}

#' @import reticulate
#' @export
call_python_program <- function(pyname){
  package_path_dir <- packages_path()
  os <- import('os')
  os$chdir(package_path_dir)
  print(pyname)
  py_dir = paste(package_path_dir,'/code/',pyname,'.py',sep="")
  source_python(py_dir)
}

Planarian_run <- function(){
  call_python_program('sc_st_gene_charge')
  call_python_program('spot_pre')
  call_python_program('GS_mapping_HVG_gene')
  call_python_program('newRegineGrowing_use_RG')
  call_python_program('comb_mapping_spot')
  call_python_program('it_final_celltype')
  print("end")
}


#' @title planarian_main
#' main
#' @description
#' If the environment is correct, you can run this program directly to get the result
#' @param conda_path Location of the conda
#' @param resolution The default program magnification is four times
#' @param thread The number of running threads, the default is 7
#' @param colname selected column name
#' @import data.table
#' @export
drim <- function(resolution=4,thread=7,colname){
  library(data.table)
  data_path <- packages_path()
  data_dir <- paste(data_path,'/data',sep="")
  if(!dir.exists(data_dir)){
    dir.create(data_dir)
  }
  setParameters(resolution,thread,colname)
  Planarian_run()
  iterative_mapping_result_celltype_it_dir=paste(data_path,'/data/',resolution,'/mapping_result.csv',sep='')
  iterative_mapping_result_celltype_it=fread(input = iterative_mapping_result_celltype_it_dir)
  return(iterative_mapping_result_celltype_it)
  print("over")
}



uni_name <- function(names){
  if(length(names) == length(unique(names))){
    message('no names repate!')
  }
  new_names <- sapply(1:length(names), function(n){
    loc <- which(names == names[n])
    id <- match(n, loc)
    new_name <- paste0(names[n],'_', id)
    return(new_name)
  })
  return(new_names)
}

#' @import Seurat
get_map_count <- function(sc_dat, map_dat){
  id_map <- match(map_dat$single_cell_name, colnames(sc_dat))
  counts <- GetAssayData(sc_dat[,id_map],slot='count')
  colnames(counts) <- uni_name(colnames(counts))
  return(counts)
}

#' @import Seurat
get_map_visiumV1 <- function(map_dat,cell_name){
  cell_coords <- map_dat[,c('tissue','row','col','imagerow','imagecol')]
  tissue_lowres_image <- matrix(1, max(cell_coords$row), max(cell_coords$col))
  tissue_positions_list <- data.frame(row.names = cell_name,
                                      tissue = cell_coords$tissue,
                                      row = cell_coords$row, col = cell_coords$col,
                                      imagerow = cell_coords$imagerow, imagecol = cell_coords$imagecol)
  scalefactors_json <- list(fiducial_diameter_fullres = 1,
                            tissue_hires_scalef = 1,
                            tissue_lowres_scalef = 1)
  generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE){
    if (filter.matrix) {
      tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
    }
    unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
    spot.radius <- unnormalized.radius / max(dim(x = image))
    return(new(Class = 'VisiumV1', 
               image = image, 
               scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, 
                                            fiducial = scale.factors$fiducial_diameter_fullres,
                                            hires = scale.factors$tissue_hires_scalef, 
                                            lowres = scale.factors$tissue_lowres_scalef), 
               coordinates = tissue.positions, 
               spot.radius = spot.radius))
  }
  spatialObj <- generate_spatialObj(image = tissue_lowres_image, 
                                    scale.factors = scalefactors_json, 
                                    tissue.positions = tissue_positions_list)
  return(spatialObj)
}

#' @import Seurat
map_toseurat <- function(sc_dat, map_dat, project = 'project'){
  counts <- get_map_count(sc_dat, map_dat)
  spatialObj <- get_map_visiumV1(map_dat, cell_name = colnames(counts))
  meta <- sc_dat@meta.data[map_dat$single_cell_name,-1:-3]
  rownames(meta) <- colnames(counts)
  meta$spatial_name <- map_dat$spatial_name
  obj <- CreateSeuratObject(counts = counts, project = project, assay = 'Spatial', meta.data = meta)
  DefaultAssay(spatialObj) <- 'Spatial'
  obj[['slice1']] <- spatialObj
  obj <- subset(obj, subset = nCount_Spatial > 0)
  return(obj)
}

#' @title getSeuratObject
#' @description Returns a seurat object, in testing, not recommended
#' @import data.table
#' @export
getSeuratObject <- function(sc_rds,resolution){
  library(data.table)
  data_path <- packages_path()
  iterative_mapping_result_celltype_it_dir=paste(data_path,'/data/',resolution,'/mapping_result.csv',sep='')
  map_dat <- fread(input = iterative_mapping_result_celltype_it_dir)
  obj <- map_toseurat(sc_dat = sc_rds,map_dat = map_dat)
  return(obj)
}
#' @title drimPlot
#' @description Returns the result of a simple draw
#' @import png
#' @export
drimPlot <- function(){
  library(png)
  call_python_program('draw')
  dir=packages_path()
  parameter_settings_path = paste(dir,"/data/parameter_settings.csv",sep="")
  resolution <- read.csv(parameter_settings_path)[1,1]
  simple_draw_pic_path <- paste(dir,'/data/',resolution,'/mapping_result.png',sep = '')
  simple_draw_pic <- readPNG(simple_draw_pic_path)
  return(simple_draw_pic)
}
