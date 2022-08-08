import numpy as np
import pandas as pd
from pathlib import Path
import sys
#import util_code
import os
dir_root = os.getcwd()
parameter_settings_f=pd.read_csv(dir_root+'/data/parameter_settings.csv')
# cell_type_choose= sys.argv[1]
resolution = int(parameter_settings_f['parameter'][0])
# resolution = sys.argv[1]
# resolution = int(resolution)

spot_ratio_file = dir_root+"/data/deconvolution.csv"
position_file = dir_root+"/data/st_loc.csv"
TILE_PATH = Path(dir_root+"/data/" + str(resolution) + "/intermediate_result")
TILE_PATH.mkdir(parents=True, exist_ok=True)
position_celltype_num_file = dir_root+"/data/" + str(resolution) +"/intermediate_result/position_comb_celltype_num.csv"
before_pic_file =  dir_root+"/data/" + str(resolution) +"/before_pic.pdf"
spot_ratio_csv = pd.read_csv(spot_ratio_file,index_col=0)
data = spot_ratio_csv.values
round_data = np.zeros((data.shape[0],data.shape[1]))
for i in range(data.shape[0]):
    for j in range(data.shape[1]):
        if data[i][j]<0.04:

            round_data[i][j] = 0
        else:
            round_data[i][j] = int(round(data[i][j] * resolution * resolution))

for i in range(round_data.shape[0]):

    if np.sum(round_data[i]) <= resolution * resolution :
        round_data[i,np.argmax(round_data[i])] = round_data[i,np.argmax(round_data[i])] + (resolution * resolution - int(np.sum(round_data[i])))
    if np.sum(round_data[i]) >= resolution * resolution - 1:
        round_data[i,np.argmax(round_data[i])] = round_data[i,np.argmax(round_data[i])] - (int(np.sum(round_data[i]))- resolution * resolution)

spot_round_csv = pd.DataFrame(data=round_data,index=spot_ratio_csv.index,columns=spot_ratio_csv.columns,dtype=int)
spot_round_csv["empty"] = 0
position = pd.read_csv(position_file,index_col=0)
position_celltype_num_merge = pd.merge(spot_round_csv,position,left_index=True,right_index=True)
position_celltype_num_merge["location"] = 0
position_celltype_num_merge.to_csv(position_celltype_num_file)
max_celltype = []
for i in np.argmax(spot_round_csv.values,axis=1):
    max_celltype.append(spot_ratio_csv.columns[i])
position_celltype_num_merge["celltype"] = max_celltype
#util_code.drawPicture(position_celltype_num_merge,row_name= "row",col_name="col",colorattribute="celltype",save_file=before_pic_file,is_save=True,is_show= False,save_type="pdf",point_size=8)