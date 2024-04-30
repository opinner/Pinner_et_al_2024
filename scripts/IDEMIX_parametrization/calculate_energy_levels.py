import numpy as np
import helper_functions as helper
import matplotlib.pyplot as plt
import warnings
import pandas as pd

warnings.filterwarnings("ignore")  # suppress some warnings about future code changes


def get_integration_intervals_for_tidal_peaks(freq, P, tidal_periods):
    integration_indices_intervals = []

    # get interval for each tidal frequency
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
    Sort the given iterable in the way that humans expect.
    https://stackoverflow.com/questions/2669059/how-to-sort-alpha-numeric-set-in-python
    """ 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)
    
            
# load all 7 moorings as dataframes
list_of_moorings = helper.IO.load_pickle(name="data/list_of_moorings")

# load Stratification information
N_table = pd.read_pickle("/home/ole/Desktop/CTD/N_values.pkl")


data = np.load("data/max_depth_dict.npz", allow_pickle = True)
max_depth_dict = data["max_depth_dict"].item()


list_of_moorings[2] = list_of_moorings[2].drop(labels = str(750), axis = "columns")

print("\n\n\n== Summary ==")
for mooring in list_of_moorings:
    print(mooring)
    print(f"{mooring.location = }, {mooring.time_delta = }")
    print("\n")


semidiurnal_tidal_constits = helper.Constants.get_tidal_frequencies_in_hours(tide_type = "semidiurnal")
print(semidiurnal_tidal_constits)
# sort the constituents in descending order after their tidal periods
# This is important for finding overlapping tidal peaks
semidiurnal_tidal_constits = dict(sorted(semidiurnal_tidal_constits.items(), key=lambda x: x[1], reverse=True))

TIME_BANDWIDTH_PRODUCT = 10

barotropic_energies = []
baroclinic_energies = []
IW_energies = []
extension_errors = []
latitudes = []
longitudes = []
depths = []

for nr, mooring in enumerate(list_of_moorings):

    UPPER_INTEGRATION_BOUND_IN_CPD = 10

    # pick the shallowest velocity time series
    shallowest_depth = mooring.oze.get_shallowest_time_series().name
    print(f"\nMooring {nr}, shallowest measurement at {shallowest_depth}m depth")
    
    
    # == barotropic tidal energy ==
    barotropic_energy_between_f_and_N = 0  # need for the IW energy later

    # iterate over all time series/columns in the mooring dataframe
    print(sorted_nicely(mooring.columns))
    for measurement_depth in sorted_nicely(mooring.columns):
        if measurement_depth == "time":
            continue

        # == total tidal energy ==
        tidal_energy_between_f_and_N = 0

        #print(f"\tmeasurement at {measurement_depth}m")
        complex_velocity = mooring[measurement_depth]
        complex_velocity_array = helper.Data.cut_trailing_nans(
            complex_velocity.to_numpy()
        )
        freq, velocity_spectrum = helper.Spectrum.total_multitaper(
            complex_velocity_array, dt=1 / 12, P=TIME_BANDWIDTH_PRODUCT
        )
        
        # The integral over the whole integral yield the variance of the velocity
        # The energy of a signal of mean 0 is then half the variance
        # Therefore we divide by 2 to have the corecct physical interpretation
        # of the spectrum
        kinetic_psd = velocity_spectrum / 2
        
        assert not np.any(np.isnan(kinetic_psd))

        # convert kinetic to total (kinetic+potential) energy
        # get N value at the geographic locations and depths of the velocity measurement
        # column in N_table is fixed by the mooring latitude, row by the measurement_depth

        # get instrument depth in units of meter above the sea floor
        mab_of_measurement = int(max_depth_dict[mooring.location.lon]) - int(measurement_depth) 
        #print("Test:",mab_of_measurement)

        """"
        # if measurement is too far up into the water column, log and skip it
        if mab_of_measurement > 500:
            print(f"Too far away from the sea floor: {mooring.location.lon}, {mab_of_measurement} mab.")
            #continue
        """
        
        if mab_of_measurement < 0 and mab_of_measurement > -2:
            print(f"\tInstrument depth was corrected from {mab_of_measurement} to 0 mab.")
            mab_of_measurement = 0
    
        column_name = f"({mooring.location.lat:.2f},{mooring.location.lon:.2f})"    
        avrg_N_in_rads = N_table.loc[
            N_table["mab"] == mab_of_measurement, column_name
        ].item()

        #print(f"{mooring.location.lon},{mab_of_measurement},{avrg_N_in_rads=}")

        coriolis_frequency_in_rads = helper.Constants.get_coriolis_frequency(
            mooring.location.lat, unit="rad/s", absolute=True
        )
        coriolis_frequency_in_cpd = helper.Constants.get_coriolis_frequency(
            mooring.location.lat, unit="cpd", absolute=True
        )
        avrg_N_in_cpd = avrg_N_in_rads /(2*np.pi) *86400

        
        fN_freq, fN_kinetic_psd = helper.Constants.cut_to_f_N_interval(
            freq,
            kinetic_psd, 
            f = coriolis_frequency_in_cpd,
            N = avrg_N_in_cpd
        )
        
        
        # fit constant slope to the unaltered spectrum 
        last_freq = fN_freq[-1] #highest resolved frequency
        extent_freq = np.linspace(last_freq, avrg_N_in_cpd, 20)

        start_freq = 3.2 #cpd
        start_index = np.argmin(np.abs(fN_freq-start_freq))         
        x = fN_freq[start_index:]
        y = fN_kinetic_psd[start_index:]
        from scipy.optimize import curve_fit
        def func_powerlaw(x, m, c):
            return x**m * c
        popt, pcov  = curve_fit(func_powerlaw, x, y, p0 = np.asarray([-2,1e4]))
        slope = tuple(popt)[0]
        slope_error = np.sqrt(np.diag(pcov))[0]
        
        # Correction factor in cpd
        kinetic_to_total_energy_factor = kinetic_to_total_energy(
            f = coriolis_frequency_in_cpd,
            N = avrg_N_in_cpd,
            omega = fN_freq
        )

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
        total_psd = kinetic_to_total_energy_factor * fN_kinetic_psd

        #fit prefactor to total psd, given the determined slope
        # x (frequency range) stays the same
        y = total_psd[start_index:]
        def func_powerlaw_total_psd(x, c):
            return x**slope * c
        slope_height, _pcov  = curve_fit(func_powerlaw_total_psd, x, y, p0 = np.asarray([1e4]))


        tidal_energy_per_peak = []
        barotropic_energy_per_peak = []
        # calculate integration intervals for the tidal peaks
        integration_indices_intervals = get_integration_intervals_for_tidal_peaks(
            fN_freq, P=TIME_BANDWIDTH_PRODUCT, tidal_periods = semidiurnal_tidal_constits.values()
        )

        # iterate over the tidal peaks
        for interval_tuple in integration_indices_intervals:
            a = fN_freq[interval_tuple[0]]
            b = fN_freq[interval_tuple[1]]
            peak_integral = helper.Spectrum.integrate_psd_interval(
                fN_freq, total_psd, a=a, b=b
            )
            background_height = min(
                total_psd[interval_tuple[0]], total_psd[interval_tuple[1]]
            )
            background_integral = background_height * (b - a)
            tidal_energy_per_peak.append(peak_integral - background_integral)

            if measurement_depth == shallowest_depth:
                    barotropic_energy_per_peak.append(peak_integral - background_integral)

            tidal_energy_per_peak.append(peak_integral - background_integral)
            
            # if the peak falls between f and N, is has to be later subtracted
            # from the total integral from f to N to obtain the IW energy
            if fN_freq[interval_tuple[0]] > coriolis_frequency_in_cpd:
                tidal_energy_between_f_and_N += peak_integral - background_integral


        # == Internal wave energy ==
        # integrate spectrum from f to N
        spectral_energy = helper.Spectrum.integrate_psd_interval(
            fN_freq, total_psd, a=coriolis_frequency_in_cpd, b= avrg_N_in_cpd
        )
        
        # add energy of the spectral extension up to N to the energy calculated from the data
        extension_energy = helper.Spectrum.integrate_psd_interval(extent_freq, func_powerlaw_total_psd(extent_freq,slope_height))
        #print(f"{(extension_energy/spectral_energy):.1%}")
        print(f"{column_name}, {measurement_depth}m: spectral slope = {slope:.2f}, {(extension_energy/spectral_energy):.1%} added energy")

        spectral_energy += extension_energy
        
        # unelegant implementation 
        #slope_error = np.sqrt(np.diag(pcov))[0]
        #def func_powerlaw_total_psd(x, c):
        #    return x**(slope-slope_error) * c       
        #upper_extension_energy  
       
        
        # subtract the energy from the peaks
        IW_energy = spectral_energy - tidal_energy_between_f_and_N
        # Sanity check: Internal wave energy should always be greater 0
        assert IW_energy > 0

        # go from the energies per peak to the total energy per time series
        barotropic_energy = np.sum(barotropic_energy_per_peak)
        baroclinic_energy = np.sum(tidal_energy_per_peak) - barotropic_energy
        if baroclinic_energy == 0:
            baroclinic_energy = np.nan
        if barotropic_energy == 0:
            barotropic_energy = np.nan

        # every energy level is for a individual time series
        # add levels to lists to they can be saved later together as file
        barotropic_energies.append(barotropic_energy)
        baroclinic_energies.append(baroclinic_energy)
        IW_energies.append(IW_energy)
        #extension_errors.append(extension_error)
        latitudes.append(mooring.location.lat)
        longitudes.append(mooring.location.lon)
        depths.append(measurement_depth)


# save results as nz file
np.savez(
    "./data/energy_levels",
    barotropic=barotropic_energies,
    baroclinic=baroclinic_energies,
    IW=IW_energies,
    #error = extension_errors,
    lats=latitudes,
    lons=longitudes,
    rounded_depths=depths,
)

