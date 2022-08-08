import numpy as np
from datetime import datetime 
import pandas as pd
import scanpy as sc
import anndata as ad
import sys
#import util_code
from sklearn.utils import shuffle
from numba import jit
from rich.progress import track
import os
import warnings
warnings.filterwarnings("ignore")
dir_root = os.getcwd()
parameter_settings_f=pd.read_csv(dir_root+'/data/parameter_settings.csv')
resolution = int(parameter_settings_f['parameter'][0])
celltype_column_name = str(parameter_settings_f['parameter'][1])


# resolution = sys.argv[1]
# resolution = int(resolution)
# celltype_column_name = sys.argv[2]
@jit(nopython=True)
def GS_matching(MP,WP):
    m = len(MP)
    n = len(WP)
    isManFree = [True]*m
    isWomenFree = [True]*n
    isManProposed = [[False for i in range(n)]for j in range(m)]
    match = [(-1,-1)]*m
    while(True in isManFree):
        indexM = isManFree.index(True)
        if(False in isManProposed[indexM]):
            indexW = -1  
            for i in range(len(MP[indexM])):
                w = MP[indexM][i]
                if(not isManProposed[indexM][w]):
                    indexW = w
                    break
            isManProposed[indexM][indexW] = True
            if(isWomenFree[indexW]):
                isWomenAccept = False
                for i in range(len(WP[indexW])):
                    if(WP[indexW][i] == indexM):
                        isWomenAccept = True
                if(isWomenAccept):
                    isWomenFree[indexW] = False
                    isManFree[indexM] = False
                    match[indexM] = (indexM,indexW)
            else:
                indexM1 = -1  
                isWomenAccept = False
                for i in range(len(WP[indexW])):
                    if (WP[indexW][i] == indexM):
                        isWomenAccept = True
                if(isWomenAccept):
                    for j in range(m):
                        if(match[j][1] == indexW):
                            indexM1 = j
                            break

                    if(np.where(WP[indexW] == indexM)[0][0]<np.where(WP[indexW] == indexM1)[0][0]):
                        isManFree[indexM1] = True
                        isManFree[indexM] = False
                        match[indexM] = (indexM,indexW)

    for i in range(len(match) - 1, -1, -1):
        if (match[i][0] == -1 or match[i][1] == -1):
                match.pop(i)
    return match
@jit(nopython=True) 
def distance_comp(before_gene_choose_np,after_gene_choose_np,distance_np):
    for i in range(before_gene_choose_np.shape[1]):
        for j in range(after_gene_choose_np.shape[1]):
            distance_np[i][j] = np.linalg.norm(before_gene_choose_np[:,i] - after_gene_choose_np[:,j])
    return distance_np
def comp_mapping(before_pd,after_pd ,myframe,change = 0):
    before_gene_choose_np = before_pd.values
    after_gene_choose_np = after_pd.values
    distance_np = np.zeros([len(before_pd.columns),len(after_pd.columns)])
    distance_np = distance_comp(before_gene_choose_np,after_gene_choose_np,distance_np)
    before_to_after_sort = distance_np.argsort()
    after_to_before_sort = distance_np.T.argsort()
    mapping_result = GS_matching(before_to_after_sort,after_to_before_sort)
    frame_in = pd.DataFrame(columns=["single_cell_name","single_cell_num","spatial_num","spatial_name"])
    before_columns= before_pd.columns
    after_columns = after_pd.columns
    if change == 0:
        for i in mapping_result:
            frame_in.loc[i[0]] = [before_columns[i[0]],i[0],i[1],after_columns[i[1]]]
    else :
        for i in mapping_result:
            frame_in.loc[i[0]] = [after_columns[i[1]],i[1],i[0],before_columns[i[0]]]
    myframe = myframe.append(frame_in)
    return myframe
def find_around_cell(cell,mapping_csv,single_cell_gene_exp_csv,layer = 2):
    layer_point_spot = []
    layer_point_cell = []
    for i in range(layer + 1):
        interval = (1/resolution) * i  + 0.0001
        center_col = mapping_csv.loc[cell,"col"]
        center_row = mapping_csv.loc[cell,"row"]
        around_cell_csv = mapping_csv.loc[(mapping_csv["row"] <= center_row + interval) & 
                            (mapping_csv["row"] >= center_row - interval) & 
                            (mapping_csv["col"] >= center_col - interval) &
                            (mapping_csv["col"] <= center_col + interval)
                       ]
        if(i > 0 ):
            layer_cell_csv = around_cell_csv.drop(layer_point_spot[i - 1])
            around_cell = layer_cell_csv["single_cell_name"].values
            around_spot = layer_cell_csv.index
        else:
            around_cell = around_cell_csv["single_cell_name"].values
            around_spot = around_cell_csv.index
        layer_point_cell.append(around_cell)
        layer_point_spot.append(around_spot)
    all_num = (len(layer_point_cell) + 1 ) * len(layer_point_cell) /2
    center_cell_gene_exp = single_cell_gene_exp_csv.loc[:,mapping_csv.loc[cell,"single_cell_name"]]
    around_cell_gene_exp = pd.Series(index = single_cell_gene_exp_csv.index,dtype='float64')
    around_cell_gene_exp = around_cell_gene_exp.fillna(0)
    center_cell_gene_exp = single_cell_gene_exp_csv.loc[:,mapping_csv.loc[cell,"single_cell_name"]]
    for i in range(len(layer_point_cell)):
        layer_gene_exp = single_cell_gene_exp_csv.loc[:,layer_point_cell[i]].mean(axis=1) * (len(layer_point_cell) - i)
        around_cell_gene_exp = around_cell_gene_exp + layer_gene_exp
    around_cell_gene_exp = around_cell_gene_exp 
    center_cell_fix_gene_exp = 0.5 * center_cell_gene_exp + 0.5 * around_cell_gene_exp
    return center_cell_fix_gene_exp
def find_celltype(column_name, celltype_name):
    if column_name.split("__")[2]  ==  celltype_name:
        return True
    else:
        return False
def find_celltype_name_id(celltype_name , columns):
    have_celltype_name = []
    n = 0
    for i in columns:
        if find_celltype(i,celltype_name):
            have_celltype_name.append(i)
            
        n = n + 1
    return have_celltype_name
single_cell_file =dir_root+"/data/sc_charge_exp.csv"
result_mapping_file = dir_root+"/data/" + str(resolution) +"/mapping_result.csv"
mapping_celltype_pic = dir_root+"/data/" + str(resolution) +"/mapping_result.pdf"
sc_celltype_file = dir_root+"/data/sc_celltype.csv"
result_file = dir_root+"/data/" + str(resolution) +"/intermediate_result/comb_position_mapping.csv"
sc_celltype_csv = pd.read_csv(sc_celltype_file,index_col=0)
single_cell_csv = pd.read_csv(single_cell_file, index_col=0)

single_cell_adata = ad.AnnData(single_cell_csv.values.T)
sc.pp.highly_variable_genes(single_cell_adata, flavor="seurat", n_top_genes=3000)
before_mapping_result = pd.read_csv(result_file,index_col=0)
pro = 0
iterative_num = 0
mapping_csv = before_mapping_result
hvg_gene_name = []
num = 0
for cell_serial_tag in single_cell_adata.var.highly_variable:
    if cell_serial_tag == True:
        hvg_gene_name.append(single_cell_csv.index[num])
    num = num + 1
#single_cell_csv = single_cell_csv.loc[hvg_gene_name]
mapping_csv = before_mapping_result
cluster_name = set()
for i in sc_celltype_csv.index:
    cluster_name.add(sc_celltype_csv.loc[i,celltype_column_name])
single_cell_celltype_column_num = {}
for in_celltype in cluster_name:
    single_cell_celltype_column_num[in_celltype] = sc_celltype_csv.loc[sc_celltype_csv[celltype_column_name] == in_celltype].index
spot_cell_celltype_index_num = {}
for in_celltype in cluster_name:
    in_celltype = str(in_celltype)
    spot_cell_celltype_index_num[in_celltype] = find_celltype_name_id(in_celltype,mapping_csv.index)
for step in track(sequence = range(20),transient = True):
    frame = pd.DataFrame(columns=["single_cell_name","single_cell_num","spatial_num","spatial_name"])
    before_time = datetime.now()
    for in_celltype in cluster_name:
        spatial_mapping_gene_csv = pd.DataFrame(index = single_cell_csv.index)
        single_cell_celltype = single_cell_csv.loc[:,single_cell_celltype_column_num[in_celltype]]
        in_celltype = str(in_celltype)
        in_celltype_mapping_csv = mapping_csv.loc[spot_cell_celltype_index_num[in_celltype]]
        shuffle_index = shuffle(in_celltype_mapping_csv.index)
        in_celltype_mapping_csv = in_celltype_mapping_csv.reindex(shuffle_index)
        for cell in shuffle_index:
            cell_fix_gene_exp = find_around_cell(cell,mapping_csv= mapping_csv,single_cell_gene_exp_csv= single_cell_csv,layer=2)
            spatial_mapping_gene_csv[cell] = cell_fix_gene_exp

        if len(single_cell_celltype.columns) >= len(spatial_mapping_gene_csv.columns) :
            change = 1
            frame = comp_mapping(spatial_mapping_gene_csv, single_cell_celltype, frame,change=1)
        else :
            while(len(single_cell_celltype.columns) <= len(spatial_mapping_gene_csv.columns)):

                frame = comp_mapping(single_cell_celltype, spatial_mapping_gene_csv,frame,change=0)
                already_mapping = []
                for col in frame["spatial_name"]:
                    if col in spatial_mapping_gene_csv.columns :
                        already_mapping.append(col)
                spatial_mapping_gene_csv = spatial_mapping_gene_csv.drop(labels=already_mapping, axis=1)
                spatial_num = len(spatial_mapping_gene_csv.columns)
            frame = comp_mapping(spatial_mapping_gene_csv, single_cell_celltype, frame,change=1)

    after_time = datetime.now()
    frame.index  = frame["spatial_name"]
    frame = frame.reindex(mapping_csv.index)
    single_cell_change_list = mapping_csv["single_cell_name"] == frame["single_cell_name"]
    num = 0
    for i in single_cell_change_list:
        if( i == True):
            num = num + 1
    mapping_csv["single_cell_name"] = frame["single_cell_name"].values
    pro  = num/len(single_cell_change_list)
    if(pro >= 0.95):
        break
    
#util_code.drawPicture(mapping_csv,row_name= "row",col_name="col",colorattribute="celltype",save_file=mapping_celltype_pic,is_save=True,is_show= False,save_type="pdf")
mapping_csv.to_csv(result_mapping_file)



