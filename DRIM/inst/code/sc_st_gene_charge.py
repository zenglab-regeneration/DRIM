#from math import comb
import numpy as np
import pandas as pd
import sys
import os
# time_point = sys.argv[1]
# print(time_point)
#num = 0
dir_root=os.getcwd()
pl_sc_file = dir_root+"/data/sc_exp.csv"
#pl_celltype_file = "D://SingleCellDate//planarian//serve-result//"+ time_point + "//"+ time_point + "_sc_celltype.csv"
pl_st_file = dir_root+"/data/st_exp.csv"
#cell_type_gene_file = "/home/sunhang/data/planarian/9_lineage_charge_markers.csv"
pl_sc_charge_file = dir_root+"/data/sc_charge_exp.csv"
pl_st_charge_file = dir_root+"/data/st_charge_exp.csv"
spot_ratio_file = dir_root+"/data/deconvolution.csv"
position_file = dir_root+"/data/sp_loc.csv"
sc_celltype_file = dir_root+"/data/sc_celltype.csv"
sc_celltype_csv = pd.read_csv(sc_celltype_file,index_col=0)
position = pd.read_csv(position_file,index_col=0)

pl_sc_csv = pd.read_csv(pl_sc_file,index_col=0)
pl_sc_csv.columns = sc_celltype_csv.index 
pl_st_csv = pd.read_csv(pl_st_file,index_col=0)

pl_st_csv.columns = position.index 
#取相交的gene
comb_gene = [val for val in pl_sc_csv.index if val in pl_st_csv.index]
sc_charge = pl_sc_csv.loc[comb_gene]
st_charge = pl_st_csv.loc[comb_gene]
sc_charge.to_csv(pl_sc_charge_file)
st_charge.to_csv(pl_st_charge_file)