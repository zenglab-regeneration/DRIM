import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")
from copy import deepcopy
import math
import sys
import os
dir_root = os.getcwd()
parameter_settings_f=pd.read_csv(dir_root+'/data/parameter_settings.csv')
resolution = int(parameter_settings_f['parameter'][0])

# resolution = sys.argv[1]

# if(len(sys.argv) == 2):   
#     expansion_threshold = 7
# else:
#     expansion_threshold = int(sys.argv[2])

if(len(parameter_settings_f['parameter']) == 2):   
    expansion_threshold = 7
else:
    expansion_threshold = int(parameter_settings_f['parameter'][2])

resolution = int(resolution)
position_celltype_num_merge_file =  dir_root+"/data/" + str(resolution) +"/intermediate_result/position_comb_celltype_num.csv"
position_celltype_num_merge = pd.read_csv(position_celltype_num_merge_file,index_col=0)

data = position_celltype_num_merge.values
data = data.astype(np.int16)
names = position_celltype_num_merge.index
row_serial = position_celltype_num_merge.columns.get_loc("row")
col_serial = position_celltype_num_merge.columns.get_loc("col")
emp_serial = position_celltype_num_merge.columns.get_loc("empty")
two_range = 1
point_range = two_range / resolution
min_row = np.min(data[:,row_serial]) - two_range * 10
max_row = np.max(data[:,row_serial]) + two_range * 10
min_col = np.min(data[:,col_serial]) - two_range * 10
max_col = np.max(data[:,col_serial]) + two_range * 10
if(min_row <= 0):
    length_row = abs(min_row) + max_row + 1
else:
    length_row = max_row + 1

if(min_col <= 0):
    length_col = abs(min_col) + max_col + 1
else:
    length_col = max_col + 1


all_celltype_num = emp_serial + 1
celltype_num_celltype = {
}
for nn in range(all_celltype_num):
    celltype_num_celltype[nn] = position_celltype_num_merge.columns[nn]
iteration_times = 0

csv_file =  dir_root+"/data/" + str(resolution) +"/intermediate_result/region_growing_result.csv"
class Point :
    def __init__(self, row, col,father_spot,name = "none",celltype= emp_serial,lock= True,):
        self.row = row
        self.col = col
        self.celltype = celltype
        self.lock = lock
        self.name = name
        self.father_spot = father_spot
    def set_lock_True(self):
        self.lock = True
    def set_lock_False(self):
        self.lock = False
    def set_cell_type(self,celltype):
        self.celltype = celltype
        self.set_lock_True()
    def result_get(self):
        result_out = np.zeros(3)
        result_out[0] = self.row
        result_out[1] = self.col
        result_out[2] = self.celltype
        return result_out, self.name
    def get_Spot(self):
        return self.father_spot
    def result_get(self):
        result_out = np.zeros(3)
        result_out[0] = self.row
        result_out[1] = self.col

        result_out[2] = self.celltype
        return result_out , self.name


class Spot:
        def __init__(self, center_row, center_col,name = "none", celltype_np=np.zeros([all_celltype_num]),resolution  = 4,lock_all=True):
            self.center_row = center_row
            self.center_col = center_col
            self.celltype_np = celltype_np
            self.before_celltype_np = celltype_np
            self.resolution = resolution
            self.lock_all = lock_all
            self.name = name
            self.point_all = []
            interval = 1/resolution
            start_row = self.center_row - 1/(resolution * 2) * (resolution + 1)
            start_col = self.center_col - 1/(resolution * 2) * (resolution + 1)
            for i in range(resolution):
                self.point_all.append([])
                start_row = start_row + interval
                start_col_in = start_col
                for j in range(resolution):
                    start_col_in = start_col_in + interval
                    self.point_all[i].append(Point(col = start_col_in, row = start_row, father_spot = self))
                    
            self.lock_num = resolution * resolution
        def all_celltype_set(self, celltype):
            for i in range(self.resolution):
                for j in range(self.resolution):
                    point = self.point_all[i][j]
                    point.set_cell_type(celltype)
        def activate(self, name,celltype_np):
            self.name = name
            self.before_celltype_np = deepcopy(celltype_np)
            for i in range(self.resolution):
                for j in range(self.resolution):
                    point = self.point_all[i][j]
                    point.name = name
            if len(np.nonzero(celltype_np)[0]) == 1:
                self.all_celltype_set(np.nonzero(celltype_np)[0][0])
                self.celltype_np = celltype_np
                self.celltype_np[np.nonzero(celltype_np)[0]] -= resolution * resolution
            else:
                self.celltype_np = celltype_np
                self.lock_all = False
                self.lock_num = 0
                for i in range(self.resolution):
                    for j in range(self.resolution):
                        point = self.point_all[i][j]
                        point.set_lock_False()
        def get_spot(self,in_row_num,in_col_num):
            return(self.point_all[in_row_num][in_col_num])
        def change_lock_num(self):
            self.lock_num = 0
            for i in range(self.resolution):
                for j in range(self.resolution):
                    point = self.point_all[i][j]
                    if point.lock == True:
                        self.lock_num = self.lock_num + 1

            self.change_lock_all()
        def change_lock_all(self):
            if self.lock_num == self.resolution * self.resolution:
                self.lock_all = True

        def set_point_celltype(self,in_row_num,in_col_num, celltype):
            if(self.point_all[in_row_num][in_col_num].lock == False):
                self.point_all[in_row_num][in_col_num].set_cell_type(celltype)
                self.change_lock_num()
                self.change_celltype_np(celltype)
        def change_celltype_np(self, celltype):
            self.celltype_np[celltype] -= 1
            if(self.celltype_np[celltype]<0):
                print("ERROR")
        def spot_get_result(self):
            result_out = np.zeros([self.resolution * self.resolution, 3])
            names = []
            num = 0
            for i in range(self.resolution):
                for j in range(self.resolution):

                    point = self.point_all[i][j]

                    result_out[num] = point.result_get()[0]
                    names.append(point.result_get()[1])
                    num = num + 1
            return result_out,names
def find_left_spot(spot):
    if spot.center_col <= min_col + 0 * two_range:
        return None
    return spots[spot.center_row][spot.center_col - 1 * two_range]
def find_left_left_spot(spot):
    if spot.center_col <= min_col + 2 * two_range:
        return None
    return spots[spot.center_row][spot.center_col - 2 * two_range]
def find_right_spot(spot):
    if spot.center_col >= max_col - 1 * two_range:
        return None
    return spots[spot.center_row][spot.center_col + 1 * two_range]
def find_right_right_spot(spot):
    if spot.center_col >= max_col - 3 * two_range:
        return None
    return spots[spot.center_row][spot.center_col + 2 * two_range]
def find_up_spot(spot):
    if spot.center_row <= min_row + 0 * two_range:
        return None
    return spots[spot.center_row - 1 * two_range][spot.center_col]
def find_up_up_spot(spot):
    if spot.center_row <= min_row + 2 * two_range:
        return None
    return spots[spot.center_row - 2 * two_range][spot.center_col]
def find_down_spot(spot):
    if spot.center_row >= max_row - 1 * two_range:
        return None
    return spots[spot.center_row + 1 * two_range][spot.center_col]
def find_down_down_spot(spot):
    if spot.center_row >= max_row - 3 * two_range:
        return None
    return spots[spot.center_row + 2 * two_range][spot.center_col]
def find_left_up_spot(spot):
    if spot.center_col <= min_col + 0 * two_range or  spot.center_row <= min_row + 0 * two_range:
        return None
    return spots[spot.center_row - 1 * two_range][spot.center_col - 1 * two_range]
def find_right_up_spot(spot):
    if spot.center_col >=max_col - 2 * two_range or  spot.center_row <= min_row + 0 * two_range:
        return None
    return spots[spot.center_row - 1 * two_range][spot.center_col + 1 * two_range]
def find_left_down_spot(spot):
    if spot.center_col <= min_col + 0 * two_range or  spot.center_row >= max_row - 2 * two_range :
        return None
    return spots[spot.center_row + 1 * two_range][spot.center_col - 1 * two_range]
def find_right_down_spot(spot):
    if spot.center_col >= max_col - 1 * two_range or  spot.center_row >= max_row -1 * two_range:
        return None
    return spots[spot.center_row + 1 * two_range][spot.center_col + 1 * two_range]
def find_left_spot(spot):
    if spot.center_col <= min_col + 0 * two_range:
        return None
    return spots[spot.center_row][spot.center_col - 1 * two_range]
def find_left_left_spot(spot):
    if spot.center_col <= min_col + 2 * two_range:
        return None
    return spots[spot.center_row][spot.center_col - 2 * two_range]
def find_right_spot(spot):
    if spot.center_col >= max_col - 1 * two_range:
        return None
    return spots[spot.center_row][spot.center_col + 1 * two_range]
def find_right_right_spot(spot):
    if spot.center_col >= max_col - 3 * two_range:
        return None
    return spots[spot.center_row][spot.center_col + 2 * two_range]
def find_up_spot(spot):
    if spot.center_row <= min_row + 0 * two_range:
        return None
    return spots[spot.center_row - 1 * two_range][spot.center_col]
def find_up_up_spot(spot):
    if spot.center_row <= min_row + 2 * two_range:
        return None
    return spots[spot.center_row - 2 * two_range][spot.center_col]
def find_down_spot(spot):
    if spot.center_row >= max_row - 1 * two_range:
        return None
    return spots[spot.center_row + 1 * two_range][spot.center_col]
def find_down_down_spot(spot):
    if spot.center_row >= max_row - 3 * two_range:
        return None
    return spots[spot.center_row + 2 * two_range][spot.center_col]
def find_left_up_spot(spot):
    if spot.center_col <= min_col + 0 * two_range or  spot.center_row <= min_row + 0 * two_range:
        return None
    return spots[spot.center_row - 1 * two_range][spot.center_col - 1 * two_range]
def find_right_up_spot(spot):
    if spot.center_col >=max_col - 2 * two_range or  spot.center_row <= min_row + 0 * two_range:
        return None
    return spots[spot.center_row - 1 * two_range][spot.center_col + 1 * two_range]
def find_left_down_spot(spot):
    if spot.center_col <= min_col + 0 * two_range or  spot.center_row >= max_row - 2 * two_range :
        return None
    return spots[spot.center_row + 1 * two_range][spot.center_col - 1 * two_range]
def find_right_down_spot(spot):
    if spot.center_col >= max_col - 1 * two_range or  spot.center_row >= max_row -1 * two_range:
        return None
    return spots[spot.center_row + 1 * two_range][spot.center_col + 1 * two_range]

def distance_np(m,n):
    return np.linalg.norm(m - n)
def comp_Point_Spot(point,p_spot, celltype):
    value = 0
    if(point.lock == True):
        return value 
    father_spot = point.get_Spot()
    point_np = np.array([point.row , point.col])
    spot_np = np.array([p_spot.center_row , p_spot.center_col])
    distance = distance_np(spot_np,point_np)
    value = (p_spot.before_celltype_np[celltype] * father_spot.celltype_np[celltype])/distance
    return value
def sum_around_spot(point, around_spot,celltype):
    sum_value = 0
    for i in around_spot:
        sum_value = sum_value + comp_Point_Spot(point,i,celltype)
    return sum_value
def confirm_point_celltype(point,around_spot):
    celltype_list = []
    value_list = []
    for i in np.nonzero(point.father_spot.celltype_np)[0]:
        celltype_list.append(i)
        value_list.append(sum_around_spot(point,around_spot,i))
    return [celltype_list[value_list.index(max(value_list))] ,max(value_list)]
def confirm_point_celltype_add_region_grow(point,around_spot,center_celltype):
    celltype_list = []
    value_list = []
    for i in np.nonzero(point.father_spot.celltype_np)[0]:
        celltype_list.append(i)
        if( i == center_celltype):
            value_list.append(sum_around_spot(point,around_spot,i) * 2)
        else:
            value_list.append(sum_around_spot(point,around_spot,i))
    return [celltype_list[value_list.index(max(value_list))] ,max(value_list)]
def find_around(lock_spot):
    around_spot = []
    left_up_spot = find_left_up_spot(spot = lock_spot)
    if(left_up_spot != None):
        around_spot.append(left_up_spot)
    right_up_spot = find_right_up_spot(spot =lock_spot)
    if(right_up_spot != None):
        around_spot.append(right_up_spot)


    left_down_spot = find_left_down_spot(spot=lock_spot)
    if(left_down_spot != None):
        around_spot.append(left_down_spot)


    right_down_spot = find_right_down_spot(spot=lock_spot)
    if(right_down_spot != None):
        around_spot.append(right_down_spot)


    up_up_spot = find_up_up_spot(spot= lock_spot)
    if(up_up_spot != None):
        around_spot.append(up_up_spot)

    down_down_spot = find_down_down_spot(spot= lock_spot)
    if(down_down_spot != None):
        around_spot.append(down_down_spot)


    left_left_spot = find_left_left_spot(spot= lock_spot)    
    if(left_left_spot != None):
        around_spot.append(left_left_spot)

    right_right_spot = find_right_right_spot(spot = lock_spot)
    if(right_right_spot != None):
        around_spot.append(right_right_spot)
    return around_spot
def find_around_center_spot(lock_spot):
    around_spot = []
    around_spot.append(lock_spot)
    left_spot = find_left_spot(spot = lock_spot)
    if(left_spot != None):
        around_spot.append(left_spot)
    right_spot = find_right_spot(spot =lock_spot)
    if(right_spot != None):
        around_spot.append(right_spot)


    down_spot = find_down_spot(spot=lock_spot)
    if(down_spot != None):
        around_spot.append(down_spot)

    up_spot = find_up_spot(spot= lock_spot)
    if(up_spot != None):
        around_spot.append(up_spot)
    left_up_spot = find_left_up_spot(spot = lock_spot)
    if(left_up_spot != None):
        around_spot.append(left_up_spot)
    right_up_spot = find_right_up_spot(spot =lock_spot)
    if(right_up_spot != None):
        around_spot.append(right_up_spot)


    left_down_spot = find_left_down_spot(spot=lock_spot)
    if(left_down_spot != None):
        around_spot.append(left_down_spot)


    right_down_spot = find_right_down_spot(spot=lock_spot)
    if(right_down_spot != None):
        around_spot.append(right_down_spot)
    return around_spot  
def get_spot_layer_point(spot,layer):
    point_set = set()
    for i in range(resolution - layer * 2 - 1):
        point_set.add(spot.point_all[layer][layer  + i])
        point_set.add(spot.point_all[layer + i][resolution -1  - layer])
        point_set.add(spot.point_all[resolution - 1  - layer - i][layer])
        point_set.add(spot.point_all[resolution -1 - layer][resolution- 1 -layer - i])
    return point_set 
def set_point_celltype_by_around_point(around_and_center_spots,center_point,expansion_threshold = 8):
    center_row = center_point.row
    center_col = center_point.col
    cell_type_num_np = np.zeros([all_celltype_num])
    point_name_list = []
    for i in range(all_celltype_num):
        point_name_list.append([])
    for spot in around_and_center_spots:
        for i in range(resolution):
            for j in range(resolution):
                if(spot != None):
                    choose_point =  spot.point_all[i][j]
                    choose_point_row = choose_point.row
                    choose_point_col = choose_point.col
                    if( (abs(choose_point_row - center_row) <= point_range * 1.1) and abs(choose_point_col - center_col) <= point_range *1.1):
                        choose_point_celltype = choose_point.celltype
                        cell_type_num_np[choose_point_celltype] = cell_type_num_np[choose_point_celltype] + 1
                        choose_point_name = choose_point.name
                        point_name_list[choose_point_celltype].append(choose_point_name)

    except_last = cell_type_num_np[0:all_celltype_num - 1]
    #if(np.all(except_last == 0) or cell_type_num_np[all_celltype_num - 1] >= resolution*resolution/2):
    #9:7
    #16:7
    #25:8
    if(np.all(except_last == 0) or cell_type_num_np[all_celltype_num - 1] >= expansion_threshold ):
        cell_num = cell_type_num_np.argmax()
    else:
        cell_num = except_last.argmax()
        

    return cell_num, point_name_list[cell_num][0]



spots = [[0] * length_col for _ in range(length_row)]

for i in range(min_row,max_row,two_range):
    for j in range(min_col,max_col,two_range):
        spots[i][j] = Spot(center_row = i,center_col = j,resolution = resolution)
for num in range(data.shape[0]):
    if (np.sum(data[num, 0:row_serial-1]) != resolution * resolution):
        print("ERROR")
    spots[data[num ,row_serial]][data[num , col_serial]].activate(names[num],data[num, 0:all_celltype_num])
lock_spots = []
for i in range(min_row,max_row,two_range):
    for j in range(min_col,max_col,two_range):
        if spots[i][j].lock_all == True and spots[i][j].point_all[int(resolution/2)][int(resolution/2)].celltype != emp_serial:
            lock_spots.append(spots[i][j])
before_num = 0
while(before_num < len(lock_spots)):
    iteration_times = iteration_times + 1
    before_num = len(lock_spots)

    for lock_spot in lock_spots:
        spot_max_cell_type = lock_spot.before_celltype_np.argmax()

        in_un_lock_spot = []
        left_up_spot = find_left_up_spot(spot=lock_spot)
        if(left_up_spot != None):
            in_un_lock_spot.append(left_up_spot)
        right_up_spot = find_right_up_spot(spot=lock_spot)
        if(right_up_spot != None):
            in_un_lock_spot.append(right_up_spot)
        left_down_spot = find_left_down_spot(spot=lock_spot)
        if(left_down_spot != None):
            in_un_lock_spot.append(left_down_spot)
        right_down_spot = find_right_down_spot(spot=lock_spot)
        if(right_down_spot != None):
            in_un_lock_spot.append(right_down_spot)
        for un_lock_spot in in_un_lock_spot:
            while(un_lock_spot.lock_all  != True):
                un_lock_num = []
                for i in range(resolution):
                    for j in range(resolution): 
                        if(un_lock_spot.point_all[i][j].lock == False):
                            un_lock_num.append((i,j))
                around_spots = find_around(un_lock_spot)
                point_get_celltype = np.zeros([resolution,resolution])
                point_get_value = np.zeros([resolution,resolution])
                for mm in range(len(un_lock_num)):
                    point_num = un_lock_num[mm]
                    row_num = point_num[0]
                    col_num = point_num[1]
                    unlock_point = un_lock_spot.point_all[row_num][col_num]
                    point_get_celltype[row_num][col_num],point_get_value[row_num][col_num] = confirm_point_celltype_add_region_grow(unlock_point,around_spots,spot_max_cell_type)
                if(point_get_value.max() != 0):
                    center_celltype = int(point_get_celltype[np.unravel_index(point_get_value.argmax(), point_get_value.shape)])
                    center_point_num = np.unravel_index(point_get_value.argmax(), point_get_value.shape)
                    un_lock_spot.set_point_celltype(in_row_num=center_point_num[0] ,in_col_num = center_point_num[1],celltype=center_celltype)
                else:
                    center_celltype = un_lock_spot.celltype_np.nonzero()[0][0]
                    center_point_num = un_lock_num[0]
                    un_lock_spot.set_point_celltype(in_row_num=center_point_num[0] ,in_col_num = center_point_num[1],celltype=center_celltype)
    lock_spots = []
    for i in range(min_row,max_row,two_range): 
            for j in range(min_col,max_col,two_range):
                if spots[i][j].lock_all == True and spots[i][j].point_all[int(resolution/2)][int(resolution/2)].celltype != emp_serial:
                    lock_spots.append(spots[i][j])
un_lock_spots = []
for i in range(min_row,max_row,two_range):

    for j in range(min_col,max_col,two_range):
        if spots[i][j].lock_all != True:
            un_lock_spots.append(spots[i][j])
for un_lock_spot in un_lock_spots:
    while(un_lock_spot.lock_all  != True):
        un_lock_num = []
        for i in range(resolution):
            for j in range(resolution): 
                if(un_lock_spot.point_all[i][j].lock == False):
                    un_lock_num.append((i,j))
        around_spots = find_around(un_lock_spot)
        point_get_celltype = np.zeros([resolution,resolution])
        point_get_value = np.zeros([resolution,resolution])
        for mm in range(len(un_lock_num)):
            point_num = un_lock_num[mm]
            row_num = point_num[0]
            col_num = point_num[1]
            unlock_point = un_lock_spot.point_all[row_num][col_num]
            point_get_celltype[row_num][col_num],point_get_value[row_num][col_num] = confirm_point_celltype(unlock_point,around_spots)
        if(point_get_value.max() != 0):
            center_celltype = int(point_get_celltype[np.unravel_index(point_get_value.argmax(), point_get_value.shape)])

            center_point_num = np.unravel_index(point_get_value.argmax(), point_get_value.shape)
            un_lock_spot.set_point_celltype(in_row_num=center_point_num[0] ,in_col_num = center_point_num[1],celltype=center_celltype)
        else:
            center_celltype = un_lock_spot.celltype_np.nonzero()[0][0]
            center_point_num = un_lock_num[0]
            un_lock_spot.set_point_celltype(in_row_num=center_point_num[0] ,in_col_num = center_point_num[1],celltype=center_celltype)
lock_spots = []
for i in range(min_row,max_row,two_range):
    for j in range(min_col,max_col,two_range):
        if spots[i][j].lock_all == True and spots[i][j].point_all[int(resolution/2)][int(resolution/2)].celltype == emp_serial:
            lock_spots.append(spots[i][j])
for lock_spot in lock_spots:
    for i in range(resolution):
        for j in range(resolution):
            lock_spot.point_all[i][j].lock = False
for lock_spot in lock_spots:
    around_spots = find_around_center_spot(lock_spot)
    for layer in range(math.ceil(resolution/2)):
        point_list = list(get_spot_layer_point(lock_spot,layer))
        for point in point_list:
            in_celltype,point_name = set_point_celltype_by_around_point(around_and_center_spots= around_spots,center_point=point,expansion_threshold= expansion_threshold)
            point.set_cell_type(in_celltype)
            point.name = point_name
result_np = np.array([0, 0, 0])
in_names = []
for i in range(min_row,max_row,two_range):

    for j in range(min_col,max_col,two_range):
        result_np = np.vstack([result_np, spots[i][j].spot_get_result()[0]])
        in_names = in_names + spots[i][j].spot_get_result()[1]
result_np = np.delete(result_np, 0, 0)
result_pd = pd.DataFrame(result_np, index=in_names, columns=["row", "col", "celltype"])
celltype = []
for re in range(result_np.shape[0]):
    celltype.append(celltype_num_celltype[result_np[re][2]])
result_pd["celltype"] = celltype
result_pd = result_pd.loc[result_pd["celltype"]!= "empty"]
result_pd.to_csv(csv_file)
