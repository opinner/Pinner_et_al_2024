import src.read_CTDs
import numpy as np
import pandas as pd

def get_transect_CTDs():
    data = src.read_CTDs.get_PS129_CTD_data()

    def load_stations(path):
        transect_names, transect_numbers = np.loadtxt(path, dtype=str, delimiter="\t", unpack=True)
        translate = dict(zip(transect_names, transect_numbers))
        # print(translate)
        return translate

    path = "/media/sf_VM_Folder/figures/PS129_Plots/Weddell_Sea_Transect.txt"
    translate = load_stations(path)
    transect_names = translate.keys()
    transect_numbers = translate.values()

    print(data['Event'])

    # Filter the DataFrame
    transect_df = data[data['Event'].str.replace('PS129_', '').isin(transect_names)]
    return transect_df


transect_df = get_transect_CTDs()
assert not transect_df.empty
# load LADCP data
list_of_LADCP_casts, list_of_CTD_casts = src.read_CTDs.get_PS129_CTDs_and_LADCPs()
# only use CTDS from the Weddell Sea transect
# list_of_LADCP_casts[:] = [cast for cast in list_of_LADCP_casts if
#                           cast.location.lon in transect_df['Longitude'].to_numpy()]
# list_of_CTD_casts[:] = [cast for cast in list_of_CTD_casts if cast.location.lon in transect_df['Longitude'].to_numpy()]

# check order
for LADCP_cast, CTD_cast in zip(list_of_LADCP_casts, list_of_CTD_casts):
    if LADCP_cast.name != CTD_cast.name:
        print(LADCP_cast.name, CTD_cast.name)
        print("Wrong order")
        raise AssertionError
    print(type(LADCP_cast), type(CTD_cast))
