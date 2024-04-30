import scipy.io as sio
import numpy as np
import pandas as pd
import src.helper_functions as helper
from src.mooring import Mooring
from src.location import Location
    
# ===== load_and_sort_moorings ======
file_path = "/media/sf_VM_Folder/data/moorings/supermooring.mat"
struct = sio.loadmat(file_path)
data = struct["mooring"]
data.dtype

#load necessary data into usable variables
lat_list = np.concatenate(np.squeeze(data["LAT"])).ravel() #flattens the list of arrays into one numpy array
lon_list = np.concatenate(np.squeeze(data["LON"])).ravel() #flattens the list of arrays into one numpy array
depth_list = np.squeeze(data["DEPTH"])
pressure_list = np.squeeze(data["PRES"])

u_list = np.squeeze(data["UC"])
v_list = np.squeeze(data["VC"])
SA_list = np.squeeze(data["SA"])
Temp_list = np.squeeze(data["TEMP"])

#Convert matlab time arrays to arrays filled with python datetime objects
time_list = []
for i in range(lat_list.size):
    mtime = np.squeeze(data["DATETIME"])[i]
    time_list.append([helper.Conversions.matlab2datetime(tval) for tval in mtime.flatten().tolist()])
#time_list = np.asarray(time_list)

# Create dictionary with the depth of the deepest mooring per longitude value
max_depth_lon_dict = {}
for i in range(lat_list.size):
    _depth = np.mean(np.squeeze(depth_list[i]))
    _lon = float(np.squeeze(lon_list[i]))
    
    try:
        if _depth > max_depth_lon_dict[_lon]:
            max_depth_lon_dict[_lon] = _depth
    except KeyError:
        max_depth_lon_dict[_lon] = _depth


#create list of all locations of the time series
list_of_locations = list(zip(lat_list.flatten(), lon_list.flatten()))
#select all unique pairs of lat,lon tuples
unique_locations = list(dict.fromkeys(list_of_locations)) 
#sort after their longitude
unique_locations.sort(key=lambda x: x[1])
print(f"\n\n{unique_locations = }")
unique_locations = np.array(unique_locations)


#create as many mooring objects as unqiue locations and assign them a location
list_of_moorings = [Mooring() for _ in unique_locations]
for (mooring,(lat,lon)) in zip(list_of_moorings,unique_locations):
    mooring.location = Location(lat = lat, lon = lon)
       
#get maximum number of rows (length of time series) needed for each mooring
for mooring in list_of_moorings:
    mooring.max_length = 0
    for i in range(lon_list.size):
        #check if location is the same
        measurement_location = Location(lat = float(np.squeeze(lat_list[i])), lon = float(np.squeeze(lon_list[i])))
        if not mooring.location ==  measurement_location: continue
        
        #check existence of data
        u = np.squeeze(u_list[i])
        u_condition = not np.all(np.isnan(u))
        if not u_condition: continue

        #set new max_length 
        if mooring.max_length < len(u):
            mooring.max_length = len(u)


print("\n\n== Loaded velocity measurements ==")
#use synchronization variable i to iterate over all time series u,v,lat,lon,depth,time
for i in range(lon_list.size):

    u = np.squeeze(u_list[i])
    #Condition if the time series contains actual velocity data
    u_condition = not np.all(np.isnan(u))
    if not u_condition: 
        #print(f"Time Series {i:02} contains no velocity data")
        continue
    
    if np.sum(~np.isnan(u)) < 100:
        print(f"Time Series {i:2}, Mooring {mooring_index} at {mooring.location}, skipped as only {np.sum(~np.isnan(u))} not-null data points") 
        continue
    
    #get data from the lists   
    v = np.squeeze(v_list[i])
    time = pd.Series(time_list[i])
    time_delta = time.diff()
    mean_time_delta = time_delta.mean(numeric_only = False)
    mean_depth = np.nanmean(np.squeeze(depth_list[i]))
    measurement_location = Location(lat = float(np.squeeze(lat_list[i])), lon = float(np.squeeze(lon_list[i])))
    
    #measurement should be deeper than 10% of the water column
    if mean_depth/max_depth_lon_dict[measurement_location.lon] < 0.1: 
        print(f"Time Series {i:2}, Mooring {mooring_index} at {mooring.location}, skipped as {mean_depth:.0f}m depth, {mean_depth/max_depth_lon_dict[measurement_location.lon]:.1%} is too shallow") 
        continue
            
    #select correct mooring from the lat,lon location
    for mooring_index,mooring in enumerate(list_of_moorings):
        if mooring.location ==  measurement_location:
            break
    
    print(f"Time Series {i:2}, Mooring {mooring_index} at {mooring.location}, {time[0]:%H:%M %d.%m.%Y} -- {time.iloc[-1]:%H:%M %d.%m.%Y}") 
    
    #check if time delta is the same (difference of 1/1000 second is allowed)
    try:
        assert abs((mooring.time_delta - mean_time_delta).total_seconds()) < 1e-3, f"Time Difference of {(mooring.time_delta - mean_time_delta).total_seconds()} seconds"
        
    #if time delta is not defined yet, create it
    except AttributeError:
        mooring.time_delta = mean_time_delta
    
    #convert velocity from cm/s to m/s
    u = u / 100
    v = v / 100 

    
    #pad to max_length
    pad_length = mooring.max_length - len(u)
    u = np.pad(u, (0, pad_length), 'constant', constant_values=np.nan )
    v = np.pad(v, (0, pad_length), 'constant', constant_values=np.nan )
    assert len(u) == mooring.max_length

    #interpolate NaNs if there is a gap between valid data points and that gap is smaller than 20 values
    u = pd.Series(u).interpolate(limit = 20, limit_area = "inside").to_numpy()
    v = pd.Series(v).interpolate(limit = 20, limit_area = "inside").to_numpy()
                
    #add new column titled with the average depth     
    mooring[f"{mean_depth:.0f}"] = u+1j*v
         
    #Time axis of time series from each location can have dissaligned time stamps, 
    #up to 20 min in the case for Mooring 2 at (-63.66,-50.81). Therefore the time axis here consists of NaNs
    #TODO Beware! This strategy is only valid, if the values aren"t synchronized in follow up steps
    if mooring_index == 2: 
        if "time" not in mooring.columns:
            mooring["time"] = np.zeros(mooring.max_length).fill(np.datetime64('NaT'))
        continue
     
            
    #if the time axis contains any nans, skip this time data
    if np.sum(pd.isnull(time)) > 0: 
        continue
    try:
        #By appending rows in the previous padding step, NaT can be intoduced in the index/time column. 
        #if there are any, replace by the time series contaning more time stamps
        if pd.isnull(mooring["time"]).sum() > 0 and len(time) == len(mooring["time"]):
            mooring["time"] = time
            mooring.set_index("time")
        else:    
            min_length = min(len(mooring["time"]),len(time))
            assert np.all(mooring["time"][:min_length] == time[:min_length])
    
    #if there is no time axis yet, create one
    except KeyError:
        mooring["time"] = time
        mooring.set_index("time")


#reaarange dataframe columns so time is the first column
for mooring in list_of_moorings:
    if "time" in mooring.columns:
        mooring.insert(0, 'time', mooring.pop('time'))


#print summary of the sortation
print("\n\n\n== Summary ==")
for mooring in list_of_moorings:
    print(mooring)
    print(f"{mooring.location = }, {mooring.time_delta = }")
    print("\n")

print(f"-------------------------------------------------")
#save data as pickle object    
name = "../../data/mooring/list_of_moorings.pkl"
helper.IO.save_as_pickle(list_of_moorings, name = name)
print(f"List of moorings saved in \"{name}\"")

