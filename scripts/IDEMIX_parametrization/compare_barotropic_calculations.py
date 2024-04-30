import numpy as np
import src.helper as helper
import src.spectra
import matplotlib.pyplot as plt
import pandas as pd
import warnings
warnings.filterwarnings("ignore")  # suppress some warnings about future code changes
import os

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

    interval_start_points = list(zip(*integration_indices_intervals))[0]  # slicing a list of tuples is complicated
    assert is_strictly_increasing(interval_start_points)
    interval_end_points = list(zip(*integration_indices_intervals))[1]  # slicing a list of tuples is complicated
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


# load all 7 moorings as dataframes
list_of_moorings = helper.IO.load_pickle(name="../../data/mooring/list_of_moorings.pkl")


data = np.load("../../data/transect_bathymetry.npz")
bat_depth = data["depth"]
bat_lon = data["lon"]


print("\n\n\n== Summary ==")
for mooring in list_of_moorings:
    print(mooring)
    print(f"{mooring.location = }, {mooring.time_delta = }")
    print("\n")


#TIDAL_CONSTITUENTS = helper.Constants.get_tidal_frequencies_in_hours(Padman = True)
TIDAL_CONSTITUENTS = {"M2": 12.42, "S2": 12.00, "N2":12.66, "K2": 11.96, "K1": 23.93, "O1": 25.82, "P1":24.07, "Q1":26.87}
print(TIDAL_CONSTITUENTS)
# sort the constituents in descending order after their tidal periods
# This is important for finding overlapping tidal peaks
TIDAL_CONSTITUENTS = dict(sorted(TIDAL_CONSTITUENTS.items(), key=lambda x: x[1], reverse=True))
TIME_BANDWIDTH_PRODUCT = 10 #"smoothing" width used in the multitaper 


cats_df = pd.read_pickle("../../data/CATS/cats_data.pickle")

integral_results = []
fit_results = []
longitudes = []
for nr, mooring in enumerate(list_of_moorings):
    coriolis_frequency_in_cpd = helper.Constants.get_coriolis_frequency(
        mooring.location.lat, unit="cpd", absolute=True
    )
    bouyancy_frequency_in_cpd = 10


    fig, ax = plt.subplots(1)

    axins = ax.inset_axes(bounds = [0.02,0.05,0.35,0.8])
    
    
    # == barotropic tidal energy ==

    # pick the shallowest velocity time series
    shallowest_complex_velocity = mooring.oze.get_shallowest_time_series()
    shallowest_depth = shallowest_complex_velocity.name
    print(f"Mooring {nr}, shallowest measurement at {shallowest_depth}m depth")
    shallowest_complex_velocity_array = helper.Data.cut_trailing_nans(
        shallowest_complex_velocity.to_numpy()
    )
    shallowest_freq, shallowest_total_psd = src.spectra.total_multitaper(
        shallowest_complex_velocity_array, dt=1 / 12, P=TIME_BANDWIDTH_PRODUCT
    )
    assert not np.any(np.isnan(shallowest_total_psd))


    # calculate integration intervals for the tidal peaks
    barotropic_energy_per_peak = []
    integration_indices_intervals = get_integration_intervals_for_tidal_peaks(
        shallowest_freq, P=TIME_BANDWIDTH_PRODUCT, tidal_periods=TIDAL_CONSTITUENTS.values()
    )
    # iterate over each peak, integrate it and subtract the background energy, which is here simply calculated as box
    for interval_tuple in integration_indices_intervals:
        a = shallowest_freq[interval_tuple[0]]
        b = shallowest_freq[interval_tuple[1]]
        peak_integral = src.spectra.integrate_psd_interval(
            shallowest_freq, shallowest_total_psd, a=a, b=b
        )
        background_height = min(
            shallowest_total_psd[interval_tuple[0]],
            shallowest_total_psd[interval_tuple[1]],
        )
        background_integral = background_height * (b - a)
        barotropic_energy_per_peak.append(peak_integral - background_integral)
        
        axins.plot([a,b],2*[background_height],"k", zorder = 10)
        #ax.plot([a,b],2*[background_height],"k", zorder = 10)

    integrated_barotropic_energy = np.sum(barotropic_energy_per_peak)

    """
    errcalc = "wboot"
    out  = ttide.t_tide(xin = shallowest_complex_velocity, dt = 2, out_style = None, lat = mooring.location.lat, errcalc = errcalc)
    fitted_barotropic = out["xout"].flatten()
    
    fitted__barotropic_energy = np.mean(1/2 * np.abs(fitted_barotropic)**2)
            
    fitted_freq, fitted_psd = src.spectrum.total_multitaper(
        fitted_barotropic , dt=1 / 12, P=TIME_BANDWIDTH_PRODUCT)
        
    integral_results.append(integrated_barotropic_energy)
    fit_results.append(fitted__barotropic_energy)
    """
    longitudes.append(mooring.location.lon)
        

    axins.vlines(
        helper.Constants.get_coriolis_frequency(
            mooring.location.lat, unit="cpd", absolute=True
        ),
        1e-6,
        1e-1,
        color="k",
        alpha=0.5,
        linestyle="-",
        linewidth=0.75,
    )
    axins.text(1.78, 0.4 * 1e-3, "f", color="k", alpha=0.5)
    axins.vlines(
        24 / 12.4, 1e-6, 1e-1, color="k", alpha=0.5, linestyle="-", linewidth=0.75
    )
    axins.text(1.94, 0.35 * 1e-3, "M2", color="k", alpha=0.5)
    # Turn ticklabels of insets off
    axins.tick_params(
        axis="both",
        which="both",
        labelleft=False,
        labelbottom=False,
        left=False,
        bottom=False,
    )



    #select the corresponding column in the CATS dataframe    
    cats_uv = cats_df.iloc[:,nr+1].to_numpy()
    #dt = 1/24 because the model is hourly
    cats_freq, cats_psd = src.spectra.total_multitaper(
        cats_uv , dt=1 / 24, P=TIME_BANDWIDTH_PRODUCT)
        
        
    ax.loglog(shallowest_freq,shallowest_total_psd, c = "tab:red", label = "integral")
    axins.loglog(shallowest_freq,shallowest_total_psd, c = "tab:red")
    ylims = ax.get_ylim()
    #ax.loglog(fitted_freq,fitted_psd, label = "ttide", c = "tab:green", alpha = 0.5)
    ax.loglog(cats_freq,cats_psd, label = "CATS", c = "tab:blue", alpha = 0.5)
    ax.set_ylim(ylims)
    #axins.loglog(fitted_freq,fitted_psd, c = "tab:green", label = "ttide", alpha = 0.5)
    axins.loglog(cats_freq,cats_psd, c = "tab:blue", label = "CATS", alpha = 0.5)
                
    helper.Plot.add_tidal_and_f_ticks(ax, coriolis_frequency_in_cpd)
    axins.set_xbound(1.73, 2.1)
    axins.set_ybound(1e-4, 0.08)
    ax.legend()
    ax.set_title(f"Mooring {nr} at {mooring.location.pretty_print()}, {shallowest_depth}m")
    ax.set_xlabel("Frequency (cycles/day)")
    ax.set_ylabel("PSD (m$^2$/s$^2$ days)")
    path = os.path.realpath(__file__)           
    helper.Plot.path_as_footnote(fig = fig, path = path)
    fig.tight_layout()
    #plt.savefig(f"./figures/comp_Mooring{nr}", dpi=300)

    

#load barotropic energies from the CATS2008 model
barotropic_cats = cats_df.loc[:, cats_df.columns!='time'].abs().pow(2).mean()
cats_lons = [_[1] for _ in barotropic_cats.index.to_numpy()]

    
f,ax = plt.subplots(2,sharex = True)
ax[0].plot(longitudes, integral_results,"-", c = "tab:red", label = "integral")
ax[0].plot(longitudes, integral_results,"o", c = "tab:red")         
#ax[0].plot(longitudes, fit_results,"-", c = "tab:green", label = "ttide")
#ax[0].plot(longitudes, fit_results,"o", c = "tab:green") 
ax[0].plot(cats_lons, barotropic_cats.values,"-", c = "tab:blue", label = "CATS")
ax[0].plot(cats_lons, barotropic_cats.values,"o", c = "tab:blue") 

for axis in ax:
    a = axis.twinx()
    a.invert_yaxis()
    xlim = axis.get_xlim()
    a.fill_between(bat_lon,np.ones(len(bat_lon))*4900,bat_depth, color = "lightgrey")
    axis.set_xlim(xlim)
    a.set_ylim(4900,10)
    axis.set_zorder(a.get_zorder()+1)
    axis.set_frame_on(False)

ax[1].semilogy(longitudes, integral_results,"-", c = "tab:red", label = "integral")
ax[1].semilogy(longitudes, integral_results,"o", c = "tab:red")         
#ax[1].semilogy(longitudes, fit_results,"-", c = "tab:green", label = "ttide")
#ax[1].semilogy(longitudes, fit_results,"o", c = "tab:green")  
ax[1].semilogy(cats_lons, barotropic_cats.values,"-", c = "tab:blue", label = "CATS")
ax[1].semilogy(cats_lons, barotropic_cats.values,"o", c = "tab:blue") 

ax[0].legend()
ax[1].legend()

ax[0].set_ylabel("Energy / (m²/s²)")
ax[1].set_ylabel("log(Energy) / (m²/s²)")            
ax[1].set_xlabel("Longitude")
    
ax[0].set_title("Comparison of barotropic signal calculations")


path = os.path.realpath(__file__)           
helper.Plot.path_as_footnote(fig = f, path = path)
f.tight_layout()

f.savefig("./barotropic_comparison", dpi = 300)
plt.show()

