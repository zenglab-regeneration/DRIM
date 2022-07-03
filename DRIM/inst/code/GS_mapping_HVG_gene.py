import numpy as np
import pandas as pd
import pandas as pd 
import numpy as np
import scanpy as sc
import anndata as ad
import random
from numba import jit
from collections import defaultdict
import sys
import os
from sklearn.utils import shuffle
#from sympy import re 
#所有的gene找HVG gene用的是3000
#在细胞类型内部找的HVG gene用的是50
# time_point = sys.argv[1]
#resolution = sys.argv[1]
dir_root=os.getcwd()
parameter_settings_f=pd.read_csv(dir_root+'/data/parameter_settings.csv')
#resolution = int(resolution)
resolution = int(parameter_settings_f['parameter'][0])
celltype_column_name = "final_celltype"
@jit(nopython=True) # jit，numba装饰器中的一种
def GS_matching(MP,WP):
    #MP是男士的择偶排序的集合 WP是女士的
    m = len(MP)
    n = len(WP)
    #给出男士和女士是否单身的数组用以评价
    isManFree = [True]*m
    isWomenFree = [True]*n
    #男士是否向女士求过婚的表格
    isManProposed = [[False for i in range(n)]for j in range(m)]
    #print(isManProposed)
    #最后匹配得出的组合 返回结果
    match = [(-1,-1)]*m

    while(True in isManFree):
        #找到第一个单身男士的索引值
        indexM = isManFree.index(True)
        #对每个女生求婚  找到男士优先列表中还没找到对象的女士
        if(False in isManProposed[indexM]):
            indexW = -1  #找到还没被求婚的排名靠前的女士的索引
            for i in range(len(MP[indexM])):
                w = MP[indexM][i]
                if(not isManProposed[indexM][w]):
                    indexW = w
                    break
            isManProposed[indexM][indexW] = True
            if(isWomenFree[indexW]):#女士单身且她愿意接受这个男士
                isWomenAccept = False
                for i in range(len(WP[indexW])):
                    if(WP[indexW][i] == indexM):
                        isWomenAccept = True
                if(isWomenAccept):
                    isWomenFree[indexW] = False
                    isManFree[indexM] = False
                    match[indexM] = (indexM,indexW)
            else:
                indexM1 = -1   #与当前女士已匹配的男士的索引
                isWomenAccept = False
                for i in range(len(WP[indexW])):#找到这个人是否能被这个女士接受
                    if (WP[indexW][i] == indexM):
                        isWomenAccept = True
                if(isWomenAccept):
                    for j in range(m):
                        if(match[j][1] == indexW):
                            indexM1 = j
                            break
                    #print(np.where(WP[indexW] == indexM)[0][0])
                    #print(np.where(WP[indexW] == indexM1)[0][0])
                    if(np.where(WP[indexW] == indexM)[0][0]<np.where(WP[indexW] == indexM1)[0][0]):
                        isManFree[indexM1] = True
                        isManFree[indexM] = False
                        match[indexM] = (indexM,indexW)
    # 删除没有进行配对的值  即任一方值为-1
    for i in range(len(match) - 1, -1, -1):
        # 倒序循环，从最后一个元素循环到第一个元素。不能用正序循环，因为正序循环删除元素后后续的列表的长度和元素下标同时也跟着变了，len(list)是动态的。
        if (match[i][0] == -1 or match[i][1] == -1):
                match.pop(i)
    #print(match)
    return match
@jit(nopython=True) # jit，numba装饰器中的一种
def distance_comp(before_gene_choose_np,after_gene_choose_np,distance_np):
    for i in range(before_gene_choose_np.shape[1]):
        for j in range(after_gene_choose_np.shape[1]):
            distance_np[i][j] = np.linalg.norm(before_gene_choose_np[:,i] - after_gene_choose_np[:,j])
    return distance_np
def comp_mapping(before_pd,after_pd ,myframe,change = 0):
    #后面的比前面的多，而且single_cell 在前面的话，change = 0
    #myframe = pd.DataFrame(columns=["single_cell_name","single_cell_num","spatial_num","spatial_name"])
    #change = 0
    #single_data = single_pd.values
    #spatial_data = spatial_pd.values
    #distance_csv = pd.DataFrame(data=distance_np , index = spatial_pd.columns,columns=single_pd.columns)
    before_gene_choose_np = before_pd.values
    after_gene_choose_np = after_pd.values
    #print(len(before_pd.columns))
    #print(len(after_pd.columns))
    distance_np = np.zeros([len(before_pd.columns),len(after_pd.columns)])
    distance_np = distance_comp(before_gene_choose_np,after_gene_choose_np,distance_np)
    #转化为整数
    #print(len(before_pd.columns))
    #print(len(after_pd.columns))
    before_to_after_sort = distance_np.argsort()
    #转化为整数
    after_to_before_sort = distance_np.T.argsort()
    #print(single_cell_list)
    #before_time = time.time()
    mapping_result = GS_matching(before_to_after_sort,after_to_before_sort)
    #print(mapping_result)
    #after_time = time.time()
    #print(after_time - before_time)
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
# 判断一个str-column_name里面是否有celltype_name
# 返回值为True False
def find_celltype(column_name, celltype_name):
    #print(column_name)
    if column_name.find(celltype_name) > -1:
        return True
    else:
        return False
# 返回值为一个列表名字 [Epicardium A0_P1_c1318,....,]
def ed(m, n):
    return np.sqrt(np.sum((m - n) ** 2))
def find_celltype_name_id(celltype_name , columns):
    #print(columns)
    have_celltype_name = []
    n = 0
    for i in columns:
        if find_celltype(i,celltype_name):
            have_celltype_name.append(i)
            #print(i)
        n = n + 1
    return have_celltype_name
single_cell_file = dir_root+"/data/sc_charge_exp.csv"
#spatical_seg_file = "D:\\SingleCellDate\\heart\\st\\" + time + "_expand_spot.csv"
spatial_file = dir_root+"/data/st_charge_exp.csv"
#cell_type_gene_file = "/home/sunhang/data/pl_new/9_lineage_charge_markers.csv"
spot_cell_num_file = dir_root+"/data/" + str(resolution) +"/number.csv"
gm_mapping_result = dir_root+"/data/" + str(resolution) +"/mapping_result.csv"
sc_celltype_file = dir_root+"/data/sc_celltype.csv"
#cell_type_gene = pd.read_csv(cell_type_gene_file,index_col=1)
single_cell_csv = pd.read_csv(single_cell_file, index_col=0)
#spatical_seg_csv = pd.read_csv(spatical_seg_file,index_col=0)
spot_cell_num_csv = pd.read_csv(spot_cell_num_file,index_col=0)
sc_celltype_csv = pd.read_csv(sc_celltype_file,index_col=0)
#print(spatical_seg_csv)
spatial_csv = pd.read_csv(spatial_file,index_col=0)
specificity_gene = set()
cluster_name = set()
for i in sc_celltype_csv.index:
    cluster_name.add(sc_celltype_csv.loc[i,celltype_column_name])
frame = pd.DataFrame(columns=["single_cell_name","single_cell_num","spatial_num","spatial_name"])
single_cell_adata = ad.AnnData(single_cell_csv.values.T)
sc.pp.highly_variable_genes(single_cell_adata, flavor="seurat", n_top_genes=3000)


for cluster in cluster_name:
    #print(cluster)
    #cluster = "tgs1+_Neoblast"
    specificity_gene = []
    find_columns = sc_celltype_csv.loc[sc_celltype_csv[celltype_column_name] == cluster].index
    single_cell_celltype_csv = single_cell_csv.loc[:,find_columns]
    single_cell_celltype_adata = ad.AnnData(single_cell_celltype_csv.values.T)
    #sc.pp.highly_variable_genes(single_cell_celltype_adata, flavor="seurat", n_top_genes=50)
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
    #hvg_gene_num = []
    hvg_gene_name = []
    num = 0
    for cell_serial_tag in single_cell_adata.var.highly_variable:
        #print(i)
        if cell_serial_tag == True and single_cell_csv.index[num] in spatial_csv.index :
            #num = num +1
            #print(num)
            #print(two_sc_csv.index[num])
            #num = num +1
            #hvg_gene_num.append(num)
            hvg_gene_name.append(single_cell_csv.index[num])
        num = num + 1
        #print(num)
    single_cell_cut_csv = single_cell_csv.loc[specificity_gene, find_columns]
    single_cell_cut_csv =  pd.DataFrame(data = single_cell_cut_csv.values * specificity_time * 2000,columns= single_cell_cut_csv.columns, index = single_cell_cut_csv.index)
    single_cell_hvg_csv = single_cell_csv.loc[hvg_gene_name, find_columns]
    single_cell_cut_csv = single_cell_cut_csv.append(single_cell_hvg_csv)
    #find_spatial_colomn = []
    
    find_spatial_colomn = []
    spatial_cut_csv = pd.DataFrame(index = specificity_gene)
    for spatial in spot_cell_num_csv.index:
        if (spot_cell_num_csv.loc[spatial,cluster] >= 0):
            spatial_cut_csv[spatial] = spatial_csv.loc[specificity_gene, spatial]/spot_cell_num_csv.loc[spatial,cluster]
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
            # 后面的必须比前面的大
            frame = comp_mapping(single_cell_cut_csv, spatial_cut_csv, frame)
            already_mapping = []
            #num = 0
            for col in frame["spatial_name"]:
                #print(frame["single_cell_name"][num])
                if col in spatial_cut_csv.columns :
                    already_mapping.append(col)
                #num = num + 1
            spatial_cut_csv = spatial_cut_csv.drop(labels=already_mapping, axis=1)
            spatial_num = len(spatial_cut_csv.columns)
        frame = comp_mapping(spatial_cut_csv,single_cell_cut_csv,frame,change=1)
    print(cluster)
frame.to_csv(gm_mapping_result)
print(frame)