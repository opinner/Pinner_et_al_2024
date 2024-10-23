
# Set up conversion from ID number to cast name
def load_stations(path):
    transect_names, transect_numbers = np.loadtxt(path, dtype = (str), delimiter = "\t", unpack = True)
    translate = dict(zip(transect_names, transect_numbers))
    # print(translate)
    return translate

from scipy.stats import binned_statistic
def bin_to_10m_resolution(x, values, bin_center):
    bin_edges = list(bin_center - np.mean(np.diff(bin_center) ) /2)
    bin_edges.append(bin_center[-1] + np.mean(np.diff(bin_center) ) /2)
    bin_means = binned_statistic(x = x, values = values, bins= bin_edges, statistic = "mean")[0]
    return bin_means

def get_PS129_data():

    # get location of CTD files from the LADCP data files
    LADCP_DIRECTORY = f"/media/sf_VM_Folder/PS129_Plots/ladcp_profiles/"

    ladcp_paths = sorted(helper.IO.get_filepaths_from_directory(LADCP_DIRECTORY, inclusive = (".mat",)))

    # all CTD cast can be identified by a number
    ladcp_cast_numbers = [path.split('/')[-1][7:-4] for path in ladcp_paths]

    # create as many CTDCast objects as casts itself
    list_of_LADCP_cast s= [CTDCast() for _ in ladcp_cast_numbers]

    name_to_number_dict = load_stations("/media/sf_VM_Folder/PS129_Plots/conversion.txt")
    number_to_name_dict = {int(v): k for k, v in name_to_number_dict.items()}

    # for every cast object set a location and a name
    for i ,cast in enumerate(list_of_LADCP_casts):
        cast_number = ladcp_cast_numbers[i]
        cast_name = number_to_name_dict[int(cast_number)]
        cast.name = cast_name

        # iterate over LADCP data
    for path, cast in zip(ladcp_paths ,list_of_LADCP_casts):
        # print(sio.whosmat(path))
        data_struct = sio.loadmat(path)
        data = data_struct["dr"][0][0]
        lat = np.round(np.squeeze(data["lat"]) ,3)
        lon = np.round(np.squeeze(data["lon"]) ,3)
        cast["u"] = np.squeeze(data["u"]).astype("double")
        cast["v"] = np.squeeze(data["v"]).astype("double")
        cast["depth"] = np.squeeze(data["z"]).astype("double")
        spacing = np.mean(np.diff(cast["depth"]))
        cast["uz"] = np.gradient(cast["u"] ,spacing)
        cast["vz"] = np.gradient(cast["v"] ,spacing)

        cast.location = Location(lat = lat, lon = lon)
        cast.date = datetime.datetime(*map(int, np.squeeze(data["date"]))) # convert list of values to datetime object


    # create as many CTDCast objects as casts itself
    list_of_CTD_cast s= [CTDCast() for _ in ladcp_cast_numbers]

    # iterate over CTD data
    for LADCP_cast, CTD_cast in zip(list_of_LADCP_casts, list_of_CTD_casts):
        CTD_cast.name = LADCP_cast.name
        CTD_cast.location = LADCP_cast.location
        CTD_cast.date = LADCP_cast.date

        # load CTD data to that LADCP Cast name
        try:
            path = f"/media/sf_VM_Folder/PS129_Plots/ctd_profiles/dps129_{LADCP_cast.name}.cnv"

            # SP = Practical Salinity [PSU]
            # Temperature [ITS-90, °C]

            # add the cast data to the data dictionary
            pressure ,in_situ_temperature ,practical_salinity = np.genfromtxt(path, skip_header = 323, usecols = (0 ,1 ,5), unpack = True)
            CTD_depth = -1 * gsw.z_from_p(p = pressure, lat = LADCP_cast.location.lat)

            LADCP_cast["t"] = bin_to_10m_resolution(
                x = CTD_depth,
                values = in_situ_temperature,
                bin_center = LADCP_cast["depth"].to_list(),
            )

            LADCP_cast["SP"] = bin_to_10m_resolution(
                x = CTD_depth,
                values = practical_salinity,
                bin_center = LADCP_cast["depth"].to_list(),
            )

            CTD_cast["depth"] = CTD_depth
            CTD_cast["t"] = in_situ_temperature
            CTD_cast["SP"] = practical_salinity



        except ValueError as e:
            raise e
            print(f"Not able to load profile PS129_{LADCP_cast.name} at {LADCP_cast.location}")
            continue
        # print(cast.name)

    return [list_of_LADCP_casts, list_of_CTD_casts]


list_of_LADCP_casts, list_of_CTD_casts = get_PS129_data()
# remove all casts further east than 45° W
list_of_LADCP_casts[:] = [cast for cast in list_of_LADCP_casts if cast.location.lon < -46]
list_of_CTD_casts[:] = [cast for cast in list_of_CTD_casts if cast.location.lon < -46]

list_of_LADCP_casts[:] = [cast for cast in list_of_LADCP_casts if cast.location.lon > -54]
list_of_CTD_casts[:] = [cast for cast in list_of_CTD_casts if cast.location.lon > -54]

# check order
for LADCP_cast, CTD_cast in zip(list_of_LADCP_casts, list_of_CTD_casts):
    if LADCP_cast.name != CTD_cast.name:
        print(LADCP_cast.name,CTD_cast.name)
        print("Wrong order")
        raise AssertionError