
import pandas as pd
import os
import warnings
dir_root = os.getcwd()
warnings.filterwarnings("ignore")
pl_sc_file = dir_root+"/data/sc_exp.csv"
pl_st_file = dir_root+"/data/st_exp.csv"
pl_sc_charge_file = dir_root+"/data/sc_charge_exp.csv"
pl_st_charge_file = dir_root+"/data/st_charge_exp.csv"
position_file = dir_root+"/data/st_loc.csv"
sc_celltype_file = dir_root+"/data/sc_celltype.csv"
sc_celltype_csv = pd.read_csv(sc_celltype_file,index_col=0)
position = pd.read_csv(position_file,index_col=0)
pl_sc_csv = pd.read_csv(pl_sc_file,index_col=0)
pl_sc_csv.columns = sc_celltype_csv.index 
pl_st_csv = pd.read_csv(pl_st_file,index_col=0)
pl_st_csv.columns = position.index 
comb_gene = [val for val in pl_sc_csv.index if val in pl_st_csv.index]
sc_charge = pl_sc_csv.loc[comb_gene]
st_charge = pl_st_csv.loc[comb_gene]
sc_charge.to_csv(pl_sc_charge_file)
st_charge.to_csv(pl_st_charge_file)