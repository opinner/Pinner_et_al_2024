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


def kinetic_to_total_energy(f, N, omega):
    conversion = (
        2*(N**2 - f**2)
        / (N**2 - omega**2)
        * (omega**2)
        / (omega**2 + f**2)
    )
    return conversion
        
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
N_table = pd.read_pickle("./method_data/N_values.pkl")


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



SEMIDIURNAL_TIDAL_CONSTITUENTS = helper.Constants.get_tidal_frequencies_in_hours(tide_type = "semidiurnal")
# sort the constituents in descending order after their tidal periods
# This is important for finding overlapping tidal peaks
SEMIDIURNAL_TIDAL_CONSTITUENTS = dict(sorted(SEMIDIURNAL_TIDAL_CONSTITUENTS.items(), key=lambda x: x[1], reverse=True))
print(f"{SEMIDIURNAL_TIDAL_CONSTITUENTS = }")

    
barotropic_energies = [] #data-based estimation as the highest possibly barotropic energy => the lowest measured energy per water column
cats_barotropic = [] # barotropic tide prediction of the CATS model
tidal_energies = [] # barotropic + baroclinic semidurnal tidal energy
continuum_energies = [] # energy in the continuum without the energy at tidal frequencies
available_energies = [] # the energy available for local dissipation, baroclinic energy of higher modes

#extension_errors = []
latitudes = []
longitudes = []
depths = []
mabs = []

TIME_BANDWIDTH_PRODUCT = 10
UPPER_INTEGRATION_BOUND_IN_CPD = 10    
for nr, mooring in enumerate(list_of_moorings):

    # pick the shallowest velocity time series
    print(f"\nMooring {nr} at {mooring.location}, CATS data from {cats_lons[nr]}")
    
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
    # The integral over the whole integral yields the variance of the velocity
    # The energy of a signal of mean 0 is then half the variance
    # Therefore we divide by 2 to have the correct physical interpretation of the spectrum      
    cats_HKE_spectrum = cats_velocity_spectrum / 2    
    # == barotropic tidal energy ==
    cats_semidiurnal_barotropic_energy_between_f_and_N = src.spectra.integrate_psd_interval(cats_freq, cats_HKE_spectrum, a = 1.5, b = 2.5)
    
    #--------------------------------------------------------------------------------------------------
    #Calculate the barotropic tide per mooring
    # 1. Calculate the minimal energy at semidiurnal tidal frequencies in the water column 
    #    That corresponds to the measured maximum barotropic energy (assuming we are in a node of the standing wave of the baroclinic tides)
    # 2. Compare with barotropic tide prediction of the CATS model
    

    # iterate over all time series/columns in the mooring dataframe
    print("Columns: ",sorted_nicely(mooring.columns))
    
    horizontal_kinetic_energies_at_tidal_frequencies = []
    barotropic_estimation_depths = []
    
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
        
        horizontal_kinetic_energy_per_peak = [] #no physical meaning, as it does not differentiate between barotropic and baroclinic tides
        # calculate integration intervals for the tidal peaks
        integration_indices_intervals = get_integration_intervals_for_tidal_peaks(
            freq, P=TIME_BANDWIDTH_PRODUCT, tidal_periods = SEMIDIURNAL_TIDAL_CONSTITUENTS.values()
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
            horizontal_kinetic_energy_per_peak.append(peak_integral)
            
        horizontal_kinetic_energy_at_tidal_frequencies = np.sum(horizontal_kinetic_energy_per_peak)        
        horizontal_kinetic_energies_at_tidal_frequencies.append(horizontal_kinetic_energy_at_tidal_frequencies)   
        barotropic_estimation_depths.append(measurement_depth)


    tidal_energies.extend(horizontal_kinetic_energies_at_tidal_frequencies)
    # select the depth, where the lowest amount of energy was measured
    measured_maximum_barotropic_energy = np.min(horizontal_kinetic_energies_at_tidal_frequencies)
    measured_maximum_barotropic_instrument_index = np.argmin(horizontal_kinetic_energies_at_tidal_frequencies)

            
    # if CATS predicts more barotropic energy then full kinetic energy we measured 
    if cats_semidiurnal_barotropic_energy_between_f_and_N > measured_maximum_barotropic_energy:
        # take the measured energy as the barotropic tidal estimation
        semidiurnal_barotropic_kinetic_energy = measured_maximum_barotropic_energy  
        print(f"barotropic energy is taken as energy at {barotropic_estimation_depths[measured_maximum_barotropic_instrument_index]}m of {barotropic_estimation_depths}")
    # else take the CATS prediction
    else:
        semidiurnal_barotropic_kinetic_energy = cats_semidiurnal_barotropic_energy_between_f_and_N         
        print(f"barotropic energy is taken from CATS model")       


 
    
    
    #--------------------------------------------------------------------------------------------------
    #Calculate the wave energy per mooring, which is available for local dissipation
    
    # iterate over all time series/columns in the mooring dataframe
    #print("Columns: ",sorted_nicely(mooring.columns))
    for measurement_depth in sorted_nicely(mooring.columns):
        if measurement_depth == "time":
            continue

        #--------------------------------------------------------------------------------------------------
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



        # get instrument depth in units of meter above the sea floor
        mab_of_measurement = int(max_depth_dict[mooring.location.lon]) - int(measurement_depth) 

        """"
        # if measurement is too far up into the water column, log and skip it
        if mab_of_measurement > 500:
            print(f"Too far away from the sea floor: {mooring.location.lon}, {mab_of_measurement} mab.")
            #continue
        """
        
        if mab_of_measurement < 0:
            if mab_of_measurement > -2:
                print(f"\tInstrument depth was corrected from {mab_of_measurement} to 0 mab.")
                mab_of_measurement = 0
            else:
                raise AssertionError

        # get N value at the geographic locations and depths of the velocity measurement
        # column in N_table is fixed by the mooring latitude, row by the measurement_depth                
        column_name = f"({mooring.location.lat:.2f},{mooring.location.lon:.2f})"    
        avrg_N_in_rads = N_table.loc[
            N_table["mab"] == mab_of_measurement, column_name
        ].item()

        #print(f"{mooring.location.lon},{mab_of_measurement},{avrg_N_in_rads=}")

        avrg_N_in_cpd = avrg_N_in_rads /(2*np.pi) *86400

        fN_freq, resolved_HKE_spectrum_between_f_and_N = helper.Constants.cut_to_f_N_interval(
            freq,
            resolved_HKE_spectrum, 
            f = coriolis_frequency_in_cpd,
            N = avrg_N_in_cpd
        )
        
        

        
        
        
        
        #--------------------------------------------------------------------------------------------------
        # convert resolved kinetic to resolved total (kinetic+potential) energy between f and N

                
        # Correction factor in cpd
        kinetic_to_total_energy_factor = kinetic_to_total_energy(
            f = coriolis_frequency_in_cpd,
            N = avrg_N_in_cpd,
            omega = fN_freq
        )

        # correction factor cannot be larger than 1, as HKE < HKE + APE
        try:
            assert np.all(kinetic_to_total_energy_factor > 0.999)
        except AssertionError:
            plt.plot(fN_freq,kinetic_to_total_energy_factor)
            plt.show()

        # Compare to Correction factor in rad/s
        difference = kinetic_to_total_energy_factor - kinetic_to_total_energy(
            f=coriolis_frequency_in_rads, N=avrg_N_in_rads, omega = fN_freq * (2*np.pi) / 86400
        )
        assert np.all(np.abs(difference) < 1e-10)
    
        # PSD and the frequency is in units of cpd
        resolved_total_energy_spectrum_between_f_and_N = kinetic_to_total_energy_factor * resolved_HKE_spectrum_between_f_and_N





        #--------------------------------------------------------------------------------------------------
        # calculate resolved total energy between f and N only in the continuum (!!!) by subtracting the energy in the peaks
        
        # integrate spectrum from f to N to get a single number
        resolved_total_energy_between_f_and_N = src.spectra.integrate_psd_interval(
            fN_freq, resolved_total_energy_spectrum_between_f_and_N, a=coriolis_frequency_in_cpd, b= avrg_N_in_cpd
        )
        
        total_energy_per_peak = [] #no physical meaning, as it did not differentiate between barotropic and baroclinic tides
        # calculate integration intervals for the tidal peaks
        integration_indices_intervals = get_integration_intervals_for_tidal_peaks(
            fN_freq, P=TIME_BANDWIDTH_PRODUCT, tidal_periods = SEMIDIURNAL_TIDAL_CONSTITUENTS.values()
        )

        # iterate over the tidal peaks and calculate energy per peak
        for interval_tuple in integration_indices_intervals:
            a = fN_freq[interval_tuple[0]]
            b = fN_freq[interval_tuple[1]]
            peak_integral = src.spectra.integrate_psd_interval(
                fN_freq,
                resolved_total_energy_spectrum_between_f_and_N,
                a=a,
                b=b
            )
            # compare values at the edges
            background_height = min(
                resolved_total_energy_spectrum_between_f_and_N[interval_tuple[0]],
                resolved_total_energy_spectrum_between_f_and_N[interval_tuple[1]]
            )
            background_integral = background_height * (b - a)
            total_energy_per_peak.append(peak_integral - background_integral)


        # all spectral energy - peak energy
        resolved_continuums_total_energy = resolved_total_energy_between_f_and_N - np.sum(total_energy_per_peak)
        






        #--------------------------------------------------------------------------------------------------
        # Calculate total wave continuum energy by extending the resolved continuum up to N

        # fit constant slope to the unaltered HKE spectrum 
        last_freq = fN_freq[-1] #highest resolved frequency
        extent_freq = np.linspace(last_freq, avrg_N_in_cpd, 20)

        start_freq = 3.2 #cpd
        start_index = np.argmin(np.abs(fN_freq-start_freq))         
        x = fN_freq[start_index:]
        y = resolved_HKE_spectrum_between_f_and_N[start_index:]
        from scipy.optimize import curve_fit
        def func_powerlaw(x, m, c):
            return x**m * c
        popt, pcov  = curve_fit(func_powerlaw, x, y, p0 = np.asarray([-2,1e4]))
        slope = tuple(popt)[0]
        slope_error = np.sqrt(np.diag(pcov))[0]

        #fit prefactor to resolved total energy spectrum, given the determined slope
        # x (frequency range) stays the same
        y = resolved_total_energy_spectrum_between_f_and_N[start_index:]
        def func_powerlaw_total_psd(x, c):
            return x**slope * c
        slope_height, _pcov  = curve_fit(func_powerlaw_total_psd, x, y, p0 = np.asarray([1e4]))

        # add energy of the spectral extension up to N to the energy calculated from the data
        extension_total_energy = src.spectra.integrate_psd_interval(extent_freq, func_powerlaw_total_psd(extent_freq,slope_height))
        #print(f"{(extension_energy/spectral_energy):.1%}")

        continuum_total_energy = resolved_continuums_total_energy + extension_total_energy 
        print(f"{column_name}, {measurement_depth}m: spectral slope = {slope:.2f}, total extension energy is responsible for {(extension_total_energy/continuum_total_energy):.1%} of the total continuum energy")        
        
        
        
        
        
        


        #--------------------------------------------------------------------------------------------------
        # Calculate available total energy at tidal frequencies, for which the baroclinic tide at this depth has to be determined first
        
        # start again at the HKE
        horizontal_kinetic_energy_per_peak = [] #no physical meaning, as it does not differentiate between barotropic and baroclinic tides
        # calculate integration intervals for the tidal peaks
        integration_indices_intervals = get_integration_intervals_for_tidal_peaks(
            fN_freq, P=TIME_BANDWIDTH_PRODUCT, tidal_periods = SEMIDIURNAL_TIDAL_CONSTITUENTS.values()
        )

        # iterate over the tidal peaks and calculate energy at tidal frequencies (peak + background)
        for interval_tuple in integration_indices_intervals:
            a = fN_freq[interval_tuple[0]]
            b = fN_freq[interval_tuple[1]]
            peak_integral = src.spectra.integrate_psd_interval(
                fN_freq,
                resolved_HKE_spectrum_between_f_and_N,
                a=a,
                b=b
            )
            horizontal_kinetic_energy_per_peak.append(peak_integral)

        # sum over energies at tidal frequencies    
        horizontal_kinetic_energy_at_tidal_frequencies = np.sum(horizontal_kinetic_energy_per_peak)
  
        # subtract barotropic tide to get baroclinic tide
        # the barotropic tide for each mooring was calculated above
        semidiurnal_baroclinic_kinetic_energy = horizontal_kinetic_energy_at_tidal_frequencies - semidiurnal_barotropic_kinetic_energy
        
        # baroclinic energy cannot be negative and is 0 in the worst case
        assert semidiurnal_baroclinic_kinetic_energy >= 0
              
        # As this energy is contained in internal waves, it also must be converted to include the potential energy
        # use conversion factor value at 2 cpd
        semidiurnal_index = np.argmin(np.abs(fN_freq - 2))
        # assert the index is neither the start or end of the frequency array
        assert semidiurnal_index != 0 and semidiurnal_index != len(fN_freq) - 1
        baroclinic_conversion_factor = kinetic_to_total_energy_factor[semidiurnal_index]
  
        semidiurnal_baroclinic_total_energy = baroclinic_conversion_factor * semidiurnal_baroclinic_kinetic_energy
 
        # only higher modes (usually n>4) contain energy available for local dissipation
        # here at first assumed to be 30% of all baroclinic energy
        # Citation: Vic et al, 2019
        available_semidiurnal_baroclinic_energy = 0.3 * semidiurnal_baroclinic_total_energy
  
        


        #--------------------------------------------------------------------------------------------------
        # Combine all energies
        
        available_energy = continuum_total_energy + available_semidiurnal_baroclinic_energy  
        print(f"{column_name}, {measurement_depth}m: available semidiurnal baroclinic energy is responsible for {(available_semidiurnal_baroclinic_energy/available_energy):.1%} of the total available energy")           
        
        # save results
        continuum_energies.append(continuum_total_energy)
        cats_barotropic.append(cats_semidiurnal_barotropic_energy_between_f_and_N)
        barotropic_energies.append(semidiurnal_barotropic_kinetic_energy)
        #baroclinic_energies.append(available_semidiurnal_baroclinic_energy)
        available_energies.append(available_energy)
        latitudes.append(mooring.location.lat)
        longitudes.append(mooring.location.lon)
        depths.append(measurement_depth)
        mabs.append(mab_of_measurement)
    

                
print("Done")                
   
              
# save results as nz file
np.savez(
    "./method_data/results_available_energy",
    continuum = continuum_energies,
    barotropic = barotropic_energies,
    #baroclinic = available_semidiurnal_baroclinic_energy,
    available = available_energies,
    cats = cats_barotropic,
    tidal_energies = tidal_energies, # barotropic + baroclinic semidurnal tidal kinetic energy
    lat = latitudes,
    lon = longitudes,
    depth = depths,
    mab = mabs,      
)
