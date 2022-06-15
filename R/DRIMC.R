packages_path<-function(){
  packages_dir<-system.file(package = "DRIM")
  return(packages_dir)
}
#'@title env_python_set
#'@description
#'Set the path of the conda environment you use
#'Many python programs are called in our program,
#'so you need to set the path of the python interpreter,
#' and you can set it directly to the conda environment here.
#'@param py_path path of the conda environment
#'@import reticulate
#'@export
env_python_set<-function(py_path){
  # library(reticulate)
  use_condaenv(py_path)
  #use_python(py_path)
}
#'@title env_test
#'Detect environment dependencies of python
#'@description
#'You can use this function to detect if a package is missing from a dependent python environment.
#'@return bool TRUE or FALSE
#'@import reticulate
#'@export
env_test<-function(){
  library(reticulate)
  package_flag<-TRUE
  if(!py_available()){
    print("The conda is not ready")
    package_flag<-FALSE
  }
  package_detect<-c('anndata',
                    'collections',
                    'copy',
                    'datetime',
                    'math',
                    'matplotlib',
                    'numba',
                    'numpy',
                    'operator',
                    'pandas',
                    'pathlib',
                    'random',
                    'scanpy',
                    'seaborn',
                    'sklearn',
                    'time')
  for(package_name in package_detect){
    if(!py_module_available(package_name)){
      package_flag<-FALSE
      error_message<-paste(package_name,"is not ready",sep=" ")
      print(error_message)
    }
  }
  return(package_flag)
}
#'@title data_deal
#'data processing
#'@description
#'Preprocess the data to get the data we need later.
#'@param sc_rds Single-cell transcriptome data
#'@param st_rds Spatial transcriptome data
#'@param plot_data convolution data
#'@export
data_deal<-function(sc_rds,st_rds,plot_data){
  data_path<-packages_path()
  data_dir<-paste(data_path,'/data',sep="")
  if(!dir.exists(data_dir)){
    dir.create(data_dir)
  }
  sp_loc = st_rds@images$anterior1@coordinates
  # plot_data <- read.csv("/home/sunhang/data/New_Mouse_brain/deconvolution.csv", row.names = 1,header = TRUE)
  sc_exp = sc_rds@assays[["RNA"]]@data

  row.names(sc_exp)
  st_exp = st_rds@assays[["Spatial"]]@data
  row.names(st_exp)
  row.names(sc_exp)


  plot_data_dir=paste(data_path,"/data/deconvolution.csv",sep="")
  sc_exp_dir=paste(data_path,"/data/sc_exp.csv",sep="")
  st_exp_dir=paste(data_path,"/data/st_exp.csv",sep="")
  sp_loc_dir=paste(data_path,'/data/sp_loc.csv',sep="")
  sc_celltype_dir=paste(data_path,'/data/sc_celltype.csv',sep="")
  plot_data_=plot_data
  write.csv(plot_data_,plot_data_dir)
  write.csv(sc_exp,sc_exp_dir)
  write.csv(st_exp,st_exp_dir)
  write.csv(sp_loc,sp_loc_dir)
  sc_celltype = sc_rds@meta.data
  write.csv(sc_celltype,sc_celltype_dir)

}

#'@title parameter_setting
#'Hyperparameter setting
#'@description
#'Set Resolution and Cell Columns,resolution is the multiple of program amplification,
#'Cell Columns is the specified column name
#'@param resolution program magnification
#'@param colname selected column name
Parameter_settings<-function(resolution=4,colname){
  parameter_settings_csv<-c(resolution,colname)
  dir=packages_path()
  parameter_settings_path = paste(dir,"/data/parameter_settings.csv",sep="")
  write.table(parameter_settings_csv,file = parameter_settings_path,row.names = FALSE,col.names = 'parameter')
}
Planarian_run<-function(){
  now_path = getwd()
  package_path<-packages_path()
  setwd(package_path)
  print("start sc_st_gene_charge")
  source_python("sc_st_gene_charge.py")
  print("start spot_pre")
  source_python("spot_pre.py")
  print("start GS_mapping_HVG_gene")
  source_python("GS_mapping_HVG_gene.py")
  print("start newRegineGrowing_use_RG")
  source_python("newRegineGrowing_use_RG.py")
  print("start comb_mapping_spot")
  source_python("comb_mapping_spot.py")
  print("start it_final_celltype")
  source_python("it_final_celltype.py")
  setwd(now_path)
  print("end")
}
#'@title planarian_main
#'main
#'@description
#'If the environment is correct, you can run this program directly to get the result
#'@param conda_path Location of the conda
#'@param resolution The default program magnification is four times
#'@param colname selected column name
#'@export
planarian_main<-function(conda_path,resolution=4,colname){
  data_path<-packages_path()
  data_dir<-paste(data_path,'/data',sep="")
  if(!dir.exists(data_dir)){
    dir.create(data_dir)
  }
  env_python_set(py_path=conda_path)
  if(!env_test()){
    stop("The conda is not ready")
  }
  Parameter_settings(resolution,colname)
  Planarian_run()
  print("over")
}
