#'@aliases DRIM




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
      py_install(package_name, pip = T)
    }
  }
  for(package_name in package_detect){
    if(!py_module_available(package_name)){
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
#'@param sc_exp_data Single Cell Transcriptome Expression Matrix
#'@param st_exp_data Spatial transcriptome expression matrix
#'@param sc_celltype_data Cell type mete data
#'@param loc_data loc_data
#'@param plot_data convolution data
#'@export
data_deal<-function(sc_exp_data,st_exp_data,sc_celltype_data,loc_data,plot_data){
  data_path<-packages_path()
  data_dir<-paste(data_path,'/data',sep="")
  if(!dir.exists(data_dir)){
    dir.create(data_dir)
  }
  sp_loc = loc_data
  # plot_data <- read.csv("/home/sunhang/data/New_Mouse_brain/deconvolution.csv", row.names = 1,header = TRUE)
  sc_exp = sc_exp_data

  row.names(sc_exp)
  st_exp = st_exp_data
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
  sc_celltype = sc_celltype_data
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
  package_path<-packages_path()
  os<-import('os')
  os$chdir(package_path)
  print("start sc_st_gene_charge")
  sc_st_gene_charge_dir=package_path
  sc_st_gene_charge_dir=paste(package_path,"/code/sc_st_gene_charge.py",sep="")
  source_python(sc_st_gene_charge_dir)
  print("start spot_pre")
  spot_pre_dir=paste(package_path,"/code/spot_pre.py",sep="")
  source_python(spot_pre_dir)
  print("start GS_mapping_HVG_gene")
  GS_mapping_HVG_gene_pre_dir=paste(package_path,"/code/GS_mapping_HVG_gene.py",sep="")
  source_python(GS_mapping_HVG_gene_pre_dir)
  print("start newRegineGrowing_use_RG")
  newRegineGrowing_use_RG_dir=paste(package_path,"/code/newRegineGrowing_use_RG.py",sep="")
  source_python(newRegineGrowing_use_RG_dir)
  print("start comb_mapping_spot")
  comb_mapping_spot_dir=paste(package_path,"/code/comb_mapping_spot.py",sep="")
  source_python(comb_mapping_spot_dir)
  print("start it_final_celltype")
  it_final_celltype_dir=paste(package_path,"/code/it_final_celltype.py",sep="")
  source_python(it_final_celltype_dir)
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
  #读取最终的结果
  print("over")
}
