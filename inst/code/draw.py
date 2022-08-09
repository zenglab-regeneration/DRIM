import util_code
import pandas as pd
dir_root = os.getcwd()
parameter_settings_f=pd.read_csv(dir_root+'/data/parameter_settings.csv')
resolution = int(parameter_settings_f['parameter'][0])
mapping_result = pd.read_csv(dir_root+"/data/" + str(resolution) +"/intermediate_result/comb_position_mapping.csv")
mapping_celltype_pic = dir_root+"/data/" + str(resolution) +"/mapping_result.png"
util_code.drawPicture(mapping_csv,row_name= "row",col_name="col",colorattribute="celltype",save_file=mapping_celltype_pic,is_save=True,is_show= False,save_type="png")
