import pandas as pd
import sys
import warnings
import os
warnings.filterwarnings("ignore")
dir_root = os.getcwd()
parameter_settings_f=pd.read_csv(dir_root+'/data/parameter_settings.csv')
resolution = int(parameter_settings_f['parameter'][0])

# resolution = sys.argv[1]
# resolution = int(resolution)

region_growing_file =dir_root+"/data/" + str(resolution) +"/intermediate_result/region_growing_result.csv"

position_file = dir_root+"/data/st_loc.csv"

mapping_result_file = dir_root+"/data/" + str(resolution) +"/intermediate_result/first_mapping_result.csv"

result_file = dir_root+"/data/" + str(resolution) +"/intermediate_result/comb_position_mapping.csv"

result_csv = pd.read_csv(region_growing_file)
position_csv = pd.read_csv(position_file,index_col=0)
mapping_result_csv = pd.read_csv(mapping_result_file,index_col=0)
result_csv.columns=["spot","row","col","celltype"]
spot_result_csv = pd.DataFrame(columns=["spot","row","col","celltype"])

for index in result_csv.index:
    aa = result_csv.loc[index]
    if (aa["spot"] != "none" ):
        spot_result_csv = spot_result_csv.append(aa,ignore_index = True)
single_cell = []
spatial_cell = []
for num in spot_result_csv.index:
    tt = spot_result_csv.loc[num]

    spatial_name = tt["spot"]+ " "+str(tt["celltype"])

    ss = mapping_result_csv.loc[mapping_result_csv["spatial_name"] == spatial_name]["single_cell_name"]
    try:
        single_cell.append(ss[ss.index[0]])
        spatial_cell.append(spatial_name)
    except:
        print("ERROR")

spot_result_csv["single_cell_name"] = single_cell
spot_result_csv["spatial_cell_name"] = spatial_cell
index_name = []
for i in spot_result_csv.index:
    inname = str(spot_result_csv.loc[i,"row"]) + "__" + str(spot_result_csv.loc[i,"col"]) + "__" + str(spot_result_csv.loc[i,"celltype"])
    index_name.append(inname)
spot_result_csv.index = index_name
spot_result_csv.to_csv(result_file)

