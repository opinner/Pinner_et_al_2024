import numpy as np
import matplotlib.pyplot as plt
import gsw
import pandas as pd
import pathlib

# data handling
import scipy.io as sio
import datetime

# import my self written functions
from src.location import Location
from src.ctd_cast import CTDCast
import src.helper as helper

RAW_DATA_DIR = " "


# Set up conversion from ID number to cast name
def load_stations(path):
    transect_names, transect_numbers = np.loadtxt(path, dtype=str, delimiter="\t", unpack=True)
    translate = dict(zip(transect_names, transect_numbers))
    # print(translate)
    return translate


def get_PS129_CTD_data():
    # get location of CTD files from the LADCP data files
    LADCP_DIRECTORY = f"/media/sf_VM_Folder/figures/PS129_Plots/ladcp_profiles/"
    ladcp_paths = sorted(helper.IO.get_filepaths_from_directory(LADCP_DIRECTORY, inclusive=(".mat",)))
    assert len(ladcp_paths)!= 0

    # all CTD cast can be identified by a number
    ladcp_cast_numbers = [path.split('/')[-1][7:-4] for path in ladcp_paths]

    # create list of all CTD cast locations
    ctd_locations = []
    ctd_timestamps = []
    for path in ladcp_paths:
        data_struct = sio.loadmat(path)
        data = data_struct["dr"][0][0]
        lat = np.round(np.squeeze(data["lat"]), 3)
        lon = np.round(np.squeeze(data["lon"]), 3)
        ctd_locations.append(Location(lat=lat, lon=lon))
        time_stamp = datetime.datetime(*map(int, np.squeeze(data["date"])))  # convert list of values to datetime object
        ctd_timestamps.append(time_stamp)

    name_to_number_dict = load_stations("/media/sf_VM_Folder/figures/PS129_Plots/conversion.txt")
    number_to_name_dict = {int(v): k for k, v in name_to_number_dict.items()}

    # create as many CTDCast objects as casts itself
    list_of_PS129_casts = [CTDCast() for _ in ladcp_cast_numbers]

    # for every cast object set a location and a name
    for i, cast in enumerate(list_of_PS129_casts):
        cast_number = ladcp_cast_numbers[i]
        cast_name = number_to_name_dict[int(cast_number)]
        cast.name = cast_name
        cast.location = ctd_locations[i]
        cast.date = ctd_timestamps[i]

    column_names = ['Event', "Latitude", "Longitude", "Press [dbar]", "Sal", "Temp [°C]", "Absolute Salinity",
                    "Conservative Temperature", "Date/Time", "Depth water [m]", "Expedition"]
    data_dict = {name: [] for name in column_names}

    for cast in list_of_PS129_casts:
        # load actual data to that Cast name
        try:
            path = f"/media/sf_VM_Folder/figures/PS129_Plots/ctd_profiles/dps129_{cast.name}.cnv"

            # SP = Practical Salinity [PSU]
            # Temperature [ITS-90, °C]

            # add the cast data to the data dictionary
            pressure, in_situ_temperature, practical_salinity = np.genfromtxt(path, skip_header=322, usecols=(0, 1, 5),
                                                                              unpack=True)
            data_dict["Event"].extend([f"PS129_{cast.name}"] * len(pressure))
            data_dict["Latitude"].extend([cast.location.lat] * len(pressure))
            data_dict["Longitude"].extend([cast.location.lon] * len(pressure))
            data_dict["Press [dbar]"].extend(pressure)
            data_dict["Temp [°C]"].extend(in_situ_temperature)
            data_dict["Sal"].extend(practical_salinity)
            data_dict["Date/Time"].extend([cast.date] * len(pressure))
            data_dict["Expedition"].extend(["PS129"] * len(pressure))

            # add new attributes
            SA = gsw.SA_from_SP(SP=practical_salinity, p=pressure, lon=cast.location.lon, lat=cast.location.lat)
            data_dict["Absolute Salinity"].extend(SA)
            data_dict["Conservative Temperature"].extend(gsw.CT_from_t(SA=SA, t=in_situ_temperature, p=pressure))
            data_dict["Depth water [m]"].extend(np.abs(gsw.z_from_p(p=pressure, lat=cast.location.lat)))

        except ValueError as e:
            print("ValueError at ", cast.name, cast.location)
            print(e)
            continue
        # print(cast.name)
    # print(data_dict)
    for (k, v) in data_dict.items():
        print(k, len(v))
    return pd.DataFrame(data=data_dict)


def find_line_number(filename, target_string):
    with open(filename, 'r') as file:
        for line_number, line in enumerate(file):
            if line.startswith(target_string):
                return line_number
    return None


def save_Joinville_transect_CTDs_to_csv():
    columns = ['Event', 'Date/Time', 'Latitude', 'Longitude',
               'Depth water [m]', 'Press [dbar]', 'Temp [°C]', 'Sal', 'Expedition']
    CTDs = pd.DataFrame(columns=columns)

    # definition of box around the transect
    m = (63.21 - 64.22) / (53.8 - 47)
    b = 63.21 - m * 53.8
    shift = 0.14  #0.12

    PS129_CTDs = get_PS129_CTD_data()
    CTDs = CTDs.merge(PS129_CTDs, on=columns, how="outer")
    # Drop all rows that are not even close to the moorings
    CTDs.drop(CTDs[CTDs.Latitude < -64.5].index, inplace=True)
    CTDs.drop(CTDs[CTDs.Latitude > -63].index, inplace=True)
    CTDs.drop(CTDs[CTDs.Longitude < -54].index, inplace=True)
    CTDs.drop(CTDs[CTDs.Longitude > -47].index, inplace=True)

    CTDs.drop(CTDs[
                  m * CTDs.Longitude - b + shift < CTDs.Latitude
                  ].index, inplace=True)

    CTDs.drop(CTDs[
                  m * CTDs.Longitude - b - shift > CTDs.Latitude
                  ].index, inplace=True)

    CTDs.reset_index(inplace=True, drop=True)

    data_paths = helper.IO.get_filepaths_from_directory(directory="/media/sf_VM_Folder/data/CTD", inclusive=".tab",
                                                        exclusive=())
    for p in data_paths: print(p)

    for i, path in enumerate(data_paths):
        target_string = "Event\tDate/Time\tLatitude\tLongitude\tElevation [m]"
        skiprows = find_line_number(path, target_string)
        if skiprows is None:  #if skiprows failed
            # try different target string
            target_string = "Event	Type	Date/Time	Longitude"
            skiprows = find_line_number(path, target_string)
        if path == '/media/sf_VM_Folder/data/CTD/ANT-XXIV_3_phys_oce.tab':
            skiprows = 235  # test

        data = pd.read_csv(path, delimiter="\t", skiprows=skiprows)

        data = data[data.columns.intersection(columns)]

        # Drop all rows that are not even close to the moorings
        data.drop(data[data.Latitude < -64.5].index, inplace=True)
        data.drop(data[data.Latitude > -63].index, inplace=True)
        data.drop(data[data.Longitude < -54].index, inplace=True)
        data.drop(data[data.Longitude > -47].index, inplace=True)

        # skip to next file if no data points lose to teh transect remain
        if data.empty:
            print(f"no data close to the joinville transect in {path}")
            continue

        plt.plot(data.Longitude, data.Latitude, ".", color="lightgrey")

        data.drop(data[
                      m * data.Longitude - b + shift < data.Latitude
                      ].index, inplace=True)

        data.drop(data[
                      m * data.Longitude - b - shift > data.Latitude
                      ].index, inplace=True)

        data.reset_index(inplace=True, drop=True)

        data['Date/Time'] = pd.to_datetime(data['Date/Time'])  # , format='%d%b%Y:%H:%M:%S.%f')
        data['Event'] = data['Event'].astype('category')

        try:
            if data['Event'].iloc[0][4] != "/":
                current_expedition = data['Event'].iloc[0][0:5]
            else:
                current_expedition = data['Event'].iloc[0][0:4]

        except IndexError as e:
            assert data.empty
            print(f"No data on defined transect {path}, {skiprows = }")
            continue

        print(current_expedition, path)
        data['Expedition'] = current_expedition
        plt.plot(data.Longitude, data.Latitude, ".", label=path.split("/")[-1])

        print("\t", path[-25:-8], data['Event'].nunique())

        CTDs = CTDs.merge(data, on=columns, how="outer")

    CTDs['Event'] = CTDs['Event'].astype('category')
    CTDs['Expedition'] = CTDs['Expedition'].astype('category')
    CTDs = CTDs.copy()  # may prevent memory fragmentation created by merging data frames
    x = np.linspace(-54, -47, 4)
    plt.plot(x, m * x - b + shift, "--")
    plt.plot(x, m * x - b - shift)

    # try out different steps up the folder structure
    path = "data/CTD/joinville_transect_ctds.csv"

    for _ in range(5):  #try to find the right folder only five times to avoid an endless loop
        file = pathlib.Path(path)
        # Check if the file exists
        if file.exists():
            # Ask the user if they want to continue
            user_input = input(f"The file {file} already exists. Do you want to overwrite it? (y/n): ").strip().lower()

            # Only continue if the user agrees
            if user_input != 'y':
                print("Operation cancelled by the user.")
                return None

            else:
                print(f"{file} will be overwritten")

        try:
            CTDs.to_csv(path)
            print(f"saved CTDs at {path}")
        except OSError:
            path = "../" + path
            continue
        else:
            break

    return None


def load_Joinville_transect_CTDs():
    # try the file with neutral densities (matlab output)
    path = "data/CTD/joinville_transect_ctds_incl_neutral_density.csv"
    for _ in range(6):  # try to find the right folder only five times to avoid an endless loop
        try:
            CTDs = pd.read_csv(path)
        except FileNotFoundError as error:
            path = "../" + path
            continue

        else:
            print(f"loading of {path} was successful")
            print("renaming of matlab style columns")
            new_columns = ['index', 'Event', 'Latitude', 'Longitude', 'Press [dbar]', 'Sal', 'Temp [°C]', 'Absolute Salinity', 'Conservative Temperature', 'Date/Time',
               'Depth water [m]',   'Expedition', "Neutral density [kg m^-3]"]
            for i,j in zip(CTDs.columns, new_columns):
                print(i,"\t",j)
            CTDs.columns = new_columns
            CTDs.set_index("index", inplace=True)
            return CTDs


    print("CTD data including neutral densities could not be found")

    # try the file without neutral density
    path = "data/CTD/joinville_transect_ctds.csv"
    for _ in range(6):  # try to find the right folder only five times to avoid an endless loop
        try:
            CTDs = pd.read_csv(path)
        except FileNotFoundError as error:
            path = "../" + path
            continue

        else:
            print(f"loading of {path} successfull")
            if "neutral_density" not in CTDs.columns:
                print("!!! Warning: !!!\n"
                      "neutral_density is missing and has to be computed with the corresponding matlab toolbox\n"
                      "This variable may be required for some code to run successfully"
                      "")

            return CTDs

    print("joinville_transect_ctds could not be found")
    print(f"last checked at {path}")
    print(
        "Did you preprocess the raw CTD profiles first by running save_Joinville_transect_CTDs_to_csv() ?")
    raise FileNotFoundError



def bin_to_10m_resolution(x, values, bin_center):
    from scipy.stats import binned_statistic
    bin_edges = list(bin_center - np.mean(np.diff(bin_center)) / 2)
    bin_edges.append(bin_center[-1] + np.mean(np.diff(bin_center)) / 2)
    bin_means = binned_statistic(x=x, values=values, bins=bin_edges, statistic="mean")[0]
    return bin_means


def get_PS129_CTDs_and_LADCPs():
    # get location of CTD files from the LADCP data files
    LADCP_DIRECTORY = f"/media/sf_VM_Folder/figures/PS129_Plots/ladcp_profiles/"

    ladcp_paths = sorted(helper.IO.get_filepaths_from_directory(LADCP_DIRECTORY, inclusive=(".mat",)))

    # all CTD cast can be identified by a number
    ladcp_cast_numbers = [path.split('/')[-1][7:-4] for path in ladcp_paths]

    # create as many CTDCast objects as casts itself
    list_of_LADCP_casts = [CTDCast() for _ in ladcp_cast_numbers]

    name_to_number_dict = load_stations("/media/sf_VM_Folder/figures/PS129_Plots/conversion.txt")
    number_to_name_dict = {int(v): k for k, v in name_to_number_dict.items()}

    # for every cast object set a location and a name
    for i, cast in enumerate(list_of_LADCP_casts):
        cast_number = ladcp_cast_numbers[i]
        cast_name = number_to_name_dict[int(cast_number)]
        cast.name = cast_name

        # iterate over LADCP data
    for path, cast in zip(ladcp_paths, list_of_LADCP_casts):
        # print(sio.whosmat(path))
        data_struct = sio.loadmat(path)
        data = data_struct["dr"][0][0]
        lat = np.round(np.squeeze(data["lat"]), 3)
        lon = np.round(np.squeeze(data["lon"]), 3)
        cast["u"] = np.squeeze(data["u"]).astype("double")
        cast["v"] = np.squeeze(data["v"]).astype("double")
        cast["depth"] = np.squeeze(data["z"]).astype("double")
        spacing = np.mean(np.diff(cast["depth"]))
        cast["uz"] = np.gradient(cast["u"], spacing)
        cast["vz"] = np.gradient(cast["v"], spacing)

        cast.location = Location(lat=lat, lon=lon)
        cast.date = datetime.datetime(*map(int, np.squeeze(data["date"])))  # convert list of values to datetime object

    # create as many CTDCast objects as casts itself
    list_of_CTD_casts = [CTDCast() for _ in ladcp_cast_numbers]

    # iterate over CTD data
    for LADCP_cast, CTD_cast in zip(list_of_LADCP_casts, list_of_CTD_casts):
        CTD_cast.name = LADCP_cast.name
        CTD_cast.location = LADCP_cast.location
        CTD_cast.date = LADCP_cast.date

        # load CTD data to that LADCP Cast name
        try:
            path = f"/media/sf_VM_Folder/figures/PS129_Plots/ctd_profiles/dps129_{LADCP_cast.name}.cnv"

            # SP = Practical Salinity [PSU]
            # Temperature [ITS-90, °C]

            # add the cast data to the data dictionary
            pressure, in_situ_temperature, practical_salinity = np.genfromtxt(path, skip_header=323, usecols=(0, 1, 5),
                                                                              unpack=True)
            CTD_depth = -1 * gsw.z_from_p(p=pressure, lat=LADCP_cast.location.lat)

            LADCP_cast["t"] = bin_to_10m_resolution(
                x=CTD_depth,
                values=in_situ_temperature,
                bin_center=LADCP_cast["depth"].to_list(),
            )

            LADCP_cast["SP"] = bin_to_10m_resolution(
                x=CTD_depth,
                values=practical_salinity,
                bin_center=LADCP_cast["depth"].to_list(),
            )

            CTD_cast["depth"] = CTD_depth
            CTD_cast["t"] = in_situ_temperature
            CTD_cast["SP"] = practical_salinity

        except ValueError:
            print(f"Not able to load profile PS129_{LADCP_cast.name} at {LADCP_cast.location}")
            continue

    return [list_of_LADCP_casts, list_of_CTD_casts]


def profiles_per_expedition(df):
    # Group by 'expedition' and count the number of unique 'event' occurrences
    event_counts = df.groupby('Expedition')['Event'].nunique().reset_index()

    # Rename the columns for better readability
    event_counts.columns = ['Expedition', 'number_of_events']
    event_counts = event_counts.sort_values(by='Expedition').reset_index(drop=True)
    return event_counts


if __name__ == '__main__':
    try:
        CTDs = load_Joinville_transect_CTDs()
    except FileNotFoundError:
        save_Joinville_transect_CTDs_to_csv()
        CTDs = load_Joinville_transect_CTDs()

    print(profiles_per_expedition(CTDs))