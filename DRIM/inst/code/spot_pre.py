import numpy as np
import pandas as pd
from pathlib import Path
import sys
import os
#time_point = sys.argv[1]
#resolution = sys.argv[1]
dir_root=os.getcwd()
parameter_settings_f=pd.read_csv(dir_root+'/data/parameter_settings.csv')
#resolution = int(resolution)
resolution=int(parameter_settings_f['parameter'][0])
#print(time_point)
spot_ratio_file = dir_root+"/data/deconvolution.csv"
position_file = dir_root+"/data/sp_loc.csv"
#生成文件
TILE_PATH = Path(dir_root+"/data/" + str(resolution))
TILE_PATH.mkdir(parents=True, exist_ok=True)
st_num_file = dir_root+"/data/" + str(resolution) +"/number.csv"
position_celltype_num_file = dir_root+"/data/" + str(resolution) +"/position_comb_celltype_num.csv"
spot_ratio_csv = pd.read_csv(spot_ratio_file,index_col=0)
#print(spot_ratio_csv.index)

data = spot_ratio_csv.values
round_data = np.zeros((data.shape[0],data.shape[1]))
for i in range(data.shape[0]):
    #print(i)
    for j in range(data.shape[1]):
        #print(i)
        #print(j)
        #print(data[i,j])
        "change value 0.04"
        if data[i][j]<0.04:
            #print("i ="+str(i))
            #print(j)
            round_data[i][j] = 0
        else:
            round_data[i][j] = int(round(data[i][j] * resolution * resolution))
#print(round_data[799])
for i in range(round_data.shape[0]):
    #print(round_data[i])
    #print(np.sum(round_data[i]))
    #print(np.max(round_data[i]))
    if np.sum(round_data[i]) <= resolution * resolution :
        round_data[i,np.argmax(round_data[i])] = round_data[i,np.argmax(round_data[i])] + (resolution * resolution - int(np.sum(round_data[i])))
    if np.sum(round_data[i]) >= resolution * resolution - 1:
        round_data[i,np.argmax(round_data[i])] = round_data[i,np.argmax(round_data[i])] - (int(np.sum(round_data[i]))- resolution * resolution)
for i in range(round_data.shape[0]):
    #print(np.sum(round_data[i]))
    if np.sum(round_data[i]) != resolution * resolution:
        print("坏事")
charge_column = []
for i in spot_ratio_csv.columns:
    print(i)
print("*******************")
#print(spot_ratio_csv.index)

spot_round_csv = pd.DataFrame(data=round_data,index=spot_ratio_csv.index,columns=spot_ratio_csv.columns,dtype=int)
spot_round_csv["empty"] = 0
spot_round_csv.to_csv(st_num_file)
#把每行相加看看等于多少
#小于9的把最大的加一个1
position = pd.read_csv(position_file,index_col=0)

position_celltype_num_merge = pd.merge(spot_round_csv,position,left_index=True,right_index=True)
position_celltype_num_merge["location"] = 0
position_celltype_num_merge.to_csv(position_celltype_num_file)