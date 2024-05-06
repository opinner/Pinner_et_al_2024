import numpy as np
import matplotlib.pyplot as plt
import warnings
#warnings.filterwarnings("ignore")  # suppress some warnings about future code changes
import pandas as pd

import src.helper as helper
import src.spectra

def get_integration_intervals_for_tidal_peaks(freq, P, tidal_periods):
    integration_indices_intervals = []

    # get interval of width 2P for each tidal frequency
    for tidal_period in tidal_periods:
        tidal_freq_in_cpd = 24 / tidal_period
        tidal_index = np.argmin(np.abs(freq - (tidal_freq_in_cpd)))
        width = (tidal_index - P, tidal_index + P + 1)
        integration_indices_intervals.append(width)

    # because the peaks all have the same width 2P, I can check the ascending order
    # by checking if the start or the end points are in ascending order. I test here both to be extra sure
    def is_strictly_increasing(L):
        return all(x < y for x, y in zip(L, L[1:]))

    # slicing a list of tuples is complicated
    interval_start_points = list(zip(*integration_indices_intervals))[0]  
    assert is_strictly_increasing(interval_start_points)
    # slicing a list of tuples is complicated
    interval_end_points = list(zip(*integration_indices_intervals))[1]  
    assert is_strictly_increasing(interval_end_points)

    # check intervals for overlaps until none are found
    while True:
        Overlap = False
        # Iteration is over a copy to not change the object iterated over
        for i, interval in enumerate(integration_indices_intervals.copy()):
            if i + 1 == len(integration_indices_intervals):
                break
            next_interval = integration_indices_intervals[i + 1]

            start1, end1 = interval
            start2, end2 = next_interval

            if end1 >= start2 and end2 >= start1:
                Overlap = True
                integration_indices_intervals[i] = (start1, end2)
                integration_indices_intervals[i + 1] = (start1, end2)

        # remove duplicated intervals
        # by creating a dictionary with the tuples as keys and None values and then converting back to list
        integration_indices_intervals = list(
            dict.fromkeys(integration_indices_intervals)
        )

        if not Overlap:
            break
    return integration_indices_intervals

        
import re 
def sorted_nicely( l ): 
    """ 
    Sort the given iterable alphanumerically in the way that humans expect.
    https://stackoverflow.com/questions/2669059/how-to-sort-alpha-numeric-set-in-python
    """ 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)
    
            
# load all 7 moorings as dataframes
list_of_moorings = helper.IO.load_pickle(name="../../data/mooring/list_of_moorings.pkl")

# load Stratification information
N_table = pd.read_pickle("../../data/CTD/N_values.pkl")


data = np.load("../../data/max_depth_dict.npz", allow_pickle = True)
max_depth_dict = data["max_depth_dict"].item()

# drop one mooring measurement way above the sea floor
list_of_moorings[2] = list_of_moorings[2].drop(labels = str(750), axis = "columns")

"""
print("\n\n\n== Summary ==")
for mooring in list_of_moorings:
    print(mooring)
    print(f"{mooring.location = }, {mooring.time_delta = }")
print("----------------------------------------------------------\n")
"""

#load cats_data
cats_df = pd.read_pickle("../../data/CATS/cats_data.pickle")
#load barotropic energies from the CATS2008 model
cats_mean_HKE = cats_df.loc[:, cats_df.columns!='time'].abs().pow(2).mean()
cats_lons = [_[1] for _ in cats_mean_HKE.index.to_numpy()]

semidiurnal_tidal_constits = helper.Constants.get_tidal_frequencies_in_hours(tide_type = "semidiurnal")
print(semidiurnal_tidal_constits)
# sort the constituents in descending order after their tidal periods
# This is important for finding overlapping tidal peaks
semidiurnal_tidal_constits = dict(sorted(semidiurnal_tidal_constits.items(), key=lambda x: x[1], reverse=True))


    
upper_limit_barotropic_energies = [] #data-based estimation as the highest possibly barotropic energy => the lowest measured energy per water column
available_energies = [] # the energy available for local dissipation, baroclinic energy of higher modes
continuum_energies = [] # energy in the continuum without the energy at tidal frequencies
extension_errors = []
latitudes = []
longitudes = []
depths = []

fig, ax = plt.subplots(1)

TIME_BANDWIDTH_PRODUCT = 10
UPPER_INTEGRATION_BOUND_IN_CPD = 10    
for nr, mooring in enumerate(list_of_moorings):

    print(f"\nMooring {nr} at {mooring.location}, CATS data from {cats_lons[nr]}")
    fign, axn = plt.subplots(1)
    
    # assert mooring and model coordinates agree
    assert np.abs(mooring.location.lon - cats_lons[nr]) < 0.02 
    
    coriolis_frequency_in_rads = helper.Constants.get_coriolis_frequency(
        mooring.location.lat, unit="rad/s", absolute=True
        )
    coriolis_frequency_in_cpd = helper.Constants.get_coriolis_frequency(
        mooring.location.lat, unit="cpd", absolute=True
        )
        
        
        
        
        
    # Barotropic semidiurnal tidal energy, predicted from CATS, needed later 
    #select the corresponding column in the CATS dataframe , +1 is needed to skip the time column
    cats_uv = cats_df.iloc[:,nr+1].to_numpy()
    #dt = 1/24 because the model is hourly
    cats_freq, cats_velocity_spectrum = src.spectra.total_multitaper(
        cats_uv , dt= 1 / 24, P=TIME_BANDWIDTH_PRODUCT)
    # The integral over the whole integral yield the variance of the velocity
    # The energy of a signal of mean 0 is then half the variance
    # Therefore we divide by 2 to have the correct physical interpretation of the spectrum      
    cats_HKE_spectrum = cats_velocity_spectrum / 2    
    # == barotropic tidal energy ==
    cats_semidiurnal_barotropic_energy_between_f_and_N = src.spectra.integrate_psd_interval(cats_freq, cats_HKE_spectrum, a = 1.5, b = 2.5)

    axn.semilogy(cats_freq, cats_HKE_spectrum, label = "CATS")

    #--------------------------------------------------------------------------------------------------
    # Calculate the minimal energy at semidiurnal tidal frequencies in the water column 
    # That corresponds to the maximum barotropic energy (assuming we are in a node of the standing wave of the baroclinic tides)
    # iterate over all time series/columns in the mooring dataframe
    print("Columns: ",sorted_nicely(mooring.columns))
    
    horizontal_kinetic_energies_at_tidal_frequencies = []
    measurement_depths = []
    
    for measurement_depth in sorted_nicely(mooring.columns):
        if measurement_depth == "time":
            continue
        
        # Calculate resolved horizontal kinetic energy spectrum between f and N (again) 
        complex_velocity = mooring[measurement_depth]
        complex_velocity_array = helper.Data.cut_trailing_nans(
            complex_velocity.to_numpy()
        )
        freq, velocity_spectrum = src.spectra.total_multitaper(
            complex_velocity_array, dt=1 / 12, P = TIME_BANDWIDTH_PRODUCT
        )
        
        # The integral over the whole integral yield the variance of the velocity
        # The energy of a signal of mean 0 is then half the variance
        # Therefore we divide by 2 to have the correct physical interpretation of the spectrum
        resolved_HKE_spectrum = velocity_spectrum / 2
        #kinetic_psd
        assert not np.any(np.isnan(resolved_HKE_spectrum))

        axn.semilogy(freq, resolved_HKE_spectrum, label = f"{measurement_depth}")
            
        horizontal_kinetic_energy_per_peak = [] #no physical meaning, as it does not differentiate between barotropic and baroclinic tides
        # calculate integration intervals for the tidal peaks
        integration_indices_intervals = get_integration_intervals_for_tidal_peaks(
            freq, P=TIME_BANDWIDTH_PRODUCT, tidal_periods = semidiurnal_tidal_constits.values()
        )

        # iterate over the tidal peaks and calculate energy per peak
        for interval_tuple in integration_indices_intervals:
            a = freq[interval_tuple[0]]
            b = freq[interval_tuple[1]]
            peak_integral = src.spectra.integrate_psd_interval(
                freq,
                resolved_HKE_spectrum,
                a=a,
                b=b
            )
            """
            background_height = min(
                resolved_HKE_spectrum[interval_tuple[0]],
                resolved_HKE_spectrum[interval_tuple[1]]
            )
            background_integral = background_height * (b - a)
            horizontal_kinetic_energy_per_peak.append(peak_integral - background_integral)
            """
            horizontal_kinetic_energy_per_peak.append(peak_integral)
            
        horizontal_kinetic_energy_at_tidal_frequencies = np.sum(horizontal_kinetic_energy_per_peak)   
                
        horizontal_kinetic_energies_at_tidal_frequencies.append(horizontal_kinetic_energy_at_tidal_frequencies)   
        measurement_depths.append(measurement_depth)

    #print(horizontal_kinetic_energies_at_tidal_frequencies)
    maximum_barotropic_energy = np.min(horizontal_kinetic_energies_at_tidal_frequencies)
    instrument_index = np.argmin(horizontal_kinetic_energies_at_tidal_frequencies)
    measurement_depths[instrument_index]
    print(f"maximum barotropic energy found at {measurement_depths[instrument_index]}m of {measurement_depths}")

    axn.set_title(f"mooringt {nr}, at {mooring.location.lon}")
    axn.legend()
    
    #plt.plot(mooring.location.lon, maximum_barotropic_energy, "D")
    ax.plot([mooring.location.lon] * len(horizontal_kinetic_energies_at_tidal_frequencies), horizontal_kinetic_energies_at_tidal_frequencies, "o", color = "gray")
    ax.plot(mooring.location.lon, cats_semidiurnal_barotropic_energy_between_f_and_N, "x")    

                
print("Done")   
plt.show()
             
"""               
# save results as nz file
np.savez(
    "./results_data/",
    barotropic=barotropic_energies,
    baroclinic=baroclinic_energies,
    IW=IW_energies,
    #error = extension_errors,
    lats=latitudes,
    lons=longitudes,
    rounded_depths=depths,
)
"""
