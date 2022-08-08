import numpy as np
import pandas as pd
import pandas as pd 
import numpy as np
import scanpy as sc
import anndata as ad
from numba import jit
import sys
import warnings
import os
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

def find_celltype(column_name, celltype_name):
    if column_name.find(celltype_name) > -1:
        return True
    else:
        return False
def ed(m, n):
    return np.sqrt(np.sum((m - n) ** 2))
def find_celltype_name_id(celltype_name , columns):
    have_celltype_name = []
    n = 0
    for i in columns:
        if find_celltype(i,celltype_name):
            have_celltype_name.append(i)
        n = n + 1
    return have_celltype_name
single_cell_file = dir_root+"/data/sc_charge_exp.csv"
spatial_file = dir_root+"/data/st_charge_exp.csv"
spot_cell_num_file = dir_root+"/data/" + str(resolution) +"/intermediate_result/position_comb_celltype_num.csv"
gm_mapping_result = dir_root+"/data/" + str(resolution) +"/intermediate_result/first_mapping_result.csv"
sc_celltype_file = dir_root+"/data/sc_celltype.csv"
single_cell_csv = pd.read_csv(single_cell_file, index_col=0)

spot_cell_num_csv = pd.read_csv(spot_cell_num_file,index_col=0)
sc_celltype_csv = pd.read_csv(sc_celltype_file,index_col=0)
spatial_csv = pd.read_csv(spatial_file,index_col=0)
specificity_gene = set()
cluster_name = set()
for i in sc_celltype_csv.index:
    cluster_name.add(sc_celltype_csv.loc[i,celltype_column_name])
frame = pd.DataFrame(columns=["single_cell_name","single_cell_num","spatial_num","spatial_name"])
single_cell_adata = ad.AnnData(single_cell_csv.values.T)
sc.pp.highly_variable_genes(single_cell_adata, flavor="seurat", n_top_genes=3000)


for cluster in cluster_name:
    specificity_gene = []
    find_columns = sc_celltype_csv.loc[sc_celltype_csv[celltype_column_name] == cluster].index
    single_cell_celltype_csv = single_cell_csv.loc[:,find_columns]
    single_cell_celltype_adata = ad.AnnData(single_cell_celltype_csv.values.T)
    try:
        sc.pp.highly_variable_genes(single_cell_celltype_adata, flavor="seurat", n_top_genes=50)
    except:
        specificity_time = 0
    else:
        num = 0
        for cell_serial_tag in single_cell_celltype_adata.var.highly_variable:
            if cell_serial_tag == True and single_cell_celltype_csv.index[num] in spatial_csv.index:
                specificity_gene.append(single_cell_celltype_csv.index[num])
            num = num + 1
        specificity_time = (3000/len(specificity_gene))**0.5
    hvg_gene_name = []
    num = 0
    for cell_serial_tag in single_cell_adata.var.highly_variable:

        if cell_serial_tag == True and single_cell_csv.index[num] in spatial_csv.index :

            hvg_gene_name.append(single_cell_csv.index[num])
        num = num + 1
    single_cell_cut_csv = single_cell_csv.loc[specificity_gene, find_columns]
    single_cell_cut_csv =  pd.DataFrame(data = single_cell_cut_csv.values * specificity_time * 2000,columns= single_cell_cut_csv.columns, index = single_cell_cut_csv.index)
    single_cell_hvg_csv = single_cell_csv.loc[hvg_gene_name, find_columns]
    single_cell_cut_csv = single_cell_cut_csv.append(single_cell_hvg_csv)
   
    find_spatial_colomn = []
    spatial_cut_csv = pd.DataFrame(index = specificity_gene)
    emp_serial = spot_cell_num_csv.columns.get_loc("empty")
    for spatial in spot_cell_num_csv.index:
        if (spot_cell_num_csv.loc[spatial,cluster] > 0):
            spatial_cut_csv[spatial] = spatial_csv.loc[specificity_gene, spatial]/spot_cell_num_csv.loc[spatial,cluster]
            find_spatial_colomn.append(spatial)
    spatial_cut_csv =  pd.DataFrame(data = spatial_cut_csv.values * specificity_time * 2000,columns= spatial_cut_csv.columns, index = spatial_cut_csv.index)  
    spatial_hvg_csv = spatial_csv.loc[hvg_gene_name, find_spatial_colomn]
    spatial_cut_csv = spatial_cut_csv.append(spatial_hvg_csv)
    spatial_add_celltype_column = []
    for column in spatial_cut_csv.columns:
        spatial_add_celltype_column.append(column+" "+cluster)
    spatial_cut_csv.columns = spatial_add_celltype_column
    single_cell_num = len(single_cell_cut_csv.columns)
    spatial_num = len(spatial_cut_csv.columns)
    change = 0
    if single_cell_num >= spatial_num :
        change = 1
        frame = comp_mapping(spatial_cut_csv, single_cell_cut_csv, frame,change=1)
    else:
        while(single_cell_num <= spatial_num):
            frame = comp_mapping(single_cell_cut_csv, spatial_cut_csv, frame)
            already_mapping = []
            for col in frame["spatial_name"]:
                if col in spatial_cut_csv.columns :
                    already_mapping.append(col)
            spatial_cut_csv = spatial_cut_csv.drop(labels=already_mapping, axis=1)
            spatial_num = len(spatial_cut_csv.columns)
        frame = comp_mapping(spatial_cut_csv,single_cell_cut_csv,frame,change=1)
    print(cluster)
frame.to_csv(gm_mapping_result)
