


CTDs = src.read_CTDs.load_Joinville_transect_CTDs()


# Get each mooring location
#list_of_moorings = helper.IO.load_pickle(name="/home/ole/Desktop/Mooring_Analysis/energy_levels/data/list_of_moorings")
mooring_locations = [mooring.location for mooring in list_of_moorings]
# create dictionary of the mooring coordinates with empty values
dict_of_closest_ctds = {str(mooring):[] for mooring in mooring_locations}

MAXIMUM_DISTANCE_FROM_MOORING = 20 #in kilometers
NUMBER_OF_CTDs_PER_MOORING = 26

#Look only at casts and their geographic location, the cast data itself will be used later
stats = CTDs.groupby(['Event','Latitude','Longitude']).size().reset_index().rename(columns={0:'count'})
stats = stats.loc[(stats["count"]!=0)]
#drop casts that contain not enough data (alternatively looking at the max depth of each cast would also be possible)
stats.drop(stats[stats["count"] < 200].index, inplace = True) 
stats.reset_index(inplace = True, drop=True)

#Calculate Distances to the mooring locations
for i,mooring_location in enumerate(mooring_locations):
    location_tuple = (mooring_location.lat, mooring_location.lon)
    stats[f"Dist. to {mooring_location}"] = stats.apply(lambda row : geopy_distance(location_tuple, (row.Latitude, row.Longitude)).km, axis = 1)

# select the 3 closest CTD profiles
# Appending to a dataframe is bad form, but I still do it here because teh data frame remains very small

# append names of closest cast to the respective mooring key
for mooring,column_mooring_name in zip(mooring_locations,(list(stats.columns)[4:])):
    # check current mooring against the column name to see if they refer to the same mooring
    #print(str(mooring),column_mooring_name[9:])
    assert str(mooring) == column_mooring_name[9:]

    # sort the dataframe after the distances to the current mooring
    # and get the first 3 closest casts 
    temp = stats.sort_values(by = column_mooring_name)

    list_of_events = temp["Event"].iloc[0:NUMBER_OF_CTDs_PER_MOORING].to_list()
    list_of_distances = temp[column_mooring_name].iloc[0:NUMBER_OF_CTDs_PER_MOORING].to_numpy()

    #replace all cast names with NaNs that are further away from the mooring location than allowed
    is_close = np.less_equal(list_of_distances,MAXIMUM_DISTANCE_FROM_MOORING)
    list_of_events = [event if close else np.nan for (event, close) in zip(list_of_events, is_close) ]

    #for sanity checks, take a look again at the new list_of_distances
    #list_of_distances = [np.round(distance,2) if close else np.nan for (distance, close) in zip(list_of_distances, is_close) ]


    print(mooring)
    print("   ",list_of_events,"\n   ",list_of_distances)

    dict_of_closest_ctds[str(mooring)].extend(list_of_events)


np.save('./method_results/dict_of_closest_ctds.npy', dict_of_closest_ctds) 

def calculate_N_with_uncertainties(smoothing_window_size):
    f,axis = plt.subplots(ncols = len(mooring_locations), sharey = True, sharex = True, figsize = (14,6))

    comparison_z_axis = np.arange(0,600,1)

    CTDs_grouped = CTDs.groupby("Event")
    result_dict = {"mab":comparison_z_axis}
    error_dict = {"mab":comparison_z_axis}
    
    for ax, (mooring, closest_ctd_casts) in zip(axis,dict_of_closest_ctds.items()):
        ax.set_title(mooring)
        ax.axvline(0, color = "k", ls = "--", zorder = 5)
        #print(closest_ctd_casts, pd.isna(closest_ctd_casts))
        number_of_valid_profiles = (~pd.isna(closest_ctd_casts)).sum()
        ax.text(0.9,0.95, s= f"{number_of_valid_profiles} casts", ha = "right", transform=ax.transAxes, bbox=dict(facecolor='white', alpha = 0.6, edgecolor='black', boxstyle='round'))

        # for calculating mean N2 profile for each mooring
        N2_array = []
        for cast_name in closest_ctd_casts:
            if pd.isna(cast_name): continue       

            #retrieve data 
            cast = CTDs_grouped.get_group(cast_name)

            #calculate square of buoyancy frequency    
            N2, N2pressure = gsw.Nsquared(
                SA = cast["Absolute Salinity"],
                CT = cast["Conservative Temperature"],
                p = cast["Press [dbar]"],
                lat = cast["Latitude"])
            depth = -1*gsw.z_from_p(p = N2pressure, lat = cast["Latitude"].mean()) 

            #convert N2 from (rad/s)^2 to 1/s^2
            #N2 = N2_in_radians / (2*np.pi)**2
            N2 = smooth(N2,smoothing_window_size)

            # change vertical coordinate 
            # to "meters above ground" or "distance from seafloor"
            mab = max(depth) - depth
            ax.plot(np.sqrt(N2),mab, c = "lightgrey", zorder = -1)

            # interpolate to common grid
            new_N2 = np.interp(
                x = comparison_z_axis,
                xp = mab[::-1], #mab and N2 are reversed so that mab monotonically increases
                fp = N2[::-1]
                )
            N2_array.append(new_N2)

        mean_N2 = np.nanmean(N2_array, axis = 0)  #calculating the average
        median_N2 = np.nanmedian(N2_array, axis = 0)  #calculating the average
        std_N2 = np.nanstd(N2_array, axis = 0)
        d = np.abs(N2_array - np.nanmedian(N2_array, axis = 0))
        median_absolute_distance = np.nanmedian(d, axis = 0)
        
        ax.plot(np.sqrt(median_N2), comparison_z_axis, "k", lw = 2, zorder = 10) #plot median N2
        ax.plot(np.sqrt(mean_N2), comparison_z_axis, "g--", lw = 2, alpha = 0.7) #plot mean N2

        #error of N from gaussian error propagation
        delta_N_std = 1/2 * (mean_N2)**(-1/2) * std_N2
        delta_N_mad = 1/2 * (median_N2)**(-1/2) * median_absolute_distance
        
        #ax.fill_betweenx(y = comparison_z_axis, x1 = np.sqrt(mean_N2) - delta_N_std, x2 = np.sqrt(mean_N2) + delta_N_std, alpha = 0.5) #plot mean N2
        ax.fill_betweenx(y = comparison_z_axis, x1 = np.sqrt(median_N2) - delta_N_mad, x2 = np.sqrt(median_N2) + delta_N_mad, alpha = 1) #plot mean N2
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    
        result_dict[mooring] = np.sqrt(median_N2)

        error_dict[mooring] = delta_N_mad

        ax.set_xlabel("$N$ / (rad/s)")

    axis[0].set_ylabel("Meters above ground")
    ax.set_ylim(-20,550)
    ax.set_xlim(-3e-4,3e-3)
    f.tight_layout()
    
    return result_dict, error_dict

result_dict, error_dict = calculate_N_with_uncertainties(smoothing_window_size = 4)
pd.DataFrame(result_dict).to_pickle("./method_data/N_values.pkl")  
pd.DataFrame(error_dict).to_pickle("./method_data/N_std.pkl")  