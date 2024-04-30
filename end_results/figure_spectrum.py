#!/usr/bin/env python
# coding: utf-8

import numpy as np
import scipy.io as sio
import datetime
import src.helper_functions as helper
import src.spectra as spectra
from src.mooring import Mooring
from src.location import Location
from src.ctd_cast import CTDCast

import matplotlib
import matplotlib.pyplot as plt
plt.rcParams.update({
    "figure.facecolor":  "white",  
    "savefig.facecolor": "white",  
})

import warnings
import pandas as pd
import gsw
#warnings.filterwarnings("ignore")  # suppress some warnings about future code changes



ONE_COLUMN_WIDTH = 8.3
TWO_COLUMN_WIDTH = 12
GOLDEN_RATIO = 1.61
cm = 1/2.54  # centimeters in inches


# load all 7 moorings as dataframes
list_of_moorings = helper.IO.load_pickle(name="../data/mooring/list_of_moorings.pkl")
# load Stratification information
N_table = pd.read_pickle("../data/CTD/N_values.pkl")

data = np.load("../data/mooring/max_depth_dict.npz", allow_pickle = True)
max_depth_dict = data["max_depth_dict"].item()

# Get CTD Data
columns = ['Event', 'Date/Time', 'Latitude', 'Longitude',
       'Depth water [m]', 'Press [dbar]', 'Temp [°C]', 'Sal', 'Expedition']
CTDs = pd.DataFrame(columns=columns)

#definition of box around the transect 
m = (63.21 - 64.22)/(53.8 - 47)
b = 63.21 - m * 53.8
shift = 0.12

def get_PS129_CTD_data():
    #get location of CTD files from the LADCP data files
    LADCP_DIRECTORY = f"/media/sf_VM_Folder/PS129_Plots/ladcp_profiles/"

    ladcp_paths = sorted(helper.IO.get_filepaths_from_directory(LADCP_DIRECTORY, inclusive = (".mat",)))

    #all CTD cast can be identified by a number
    ladcp_cast_numbers = [path.split('/')[-1][7:-4] for path in ladcp_paths]

    #create list of all CTD cast locations
    ctd_locations = []
    ctd_timestamps = []
    for path in ladcp_paths:
        data_struct = sio.loadmat(path) 
        data = data_struct["dr"][0][0]
        lat = np.round(np.squeeze(data["lat"]),3)
        lon = np.round(np.squeeze(data["lon"]),3)
        ctd_locations.append(Location(lat = lat, lon = lon))
        time_stamp = datetime.datetime(*map(int, np.squeeze(data["date"]))) #convert list of values to datetime object
        ctd_timestamps.append(time_stamp)
        
    # Set up conversion from ID number to cast name    
    def load_stations(path):
        transect_names, transect_numbers = np.loadtxt(path, dtype = (str), delimiter = "\t", unpack = True)
        translate = dict(zip(transect_names, transect_numbers))
        #print(translate)
        return translate
    name_to_number_dict = load_stations("/media/sf_VM_Folder/PS129_Plots/conversion.txt")
    number_to_name_dict = {int(v): k for k, v in name_to_number_dict.items()}    

    #create as many CTDCast objects as casts itself
    list_of_PS129_casts= [CTDCast() for _ in ladcp_cast_numbers]

    # for every cast object set a location and a name
    for i,cast in enumerate(list_of_PS129_casts):
        cast_number = ladcp_cast_numbers[i]
        cast_name = number_to_name_dict[int(cast_number)]
        cast.name = cast_name
        cast.location = ctd_locations[i]
        cast.date = ctd_timestamps[i]

    column_names = ['Event', "Latitude", "Longitude", "Press [dbar]", "Sal", "Temp [°C]", "Absolute Salinity", "Conservative Temperature", "Date/Time", "Depth water [m]", "Expedition"]
    data_dict = {name:[] for name in column_names}

    for cast in list_of_PS129_casts:
        #load actual data to that Cast name
        try:   
            path = f"/media/sf_VM_Folder/PS129_Plots/ctd_profiles/dps129_{cast.name}.cnv"
        
            #SP = Practical Salinity [PSU]
            #Temperature [ITS-90, °C]

            #add the cast data to the data dictionary
            pressure,in_situ_temperature,practical_salinity = np.genfromtxt(path, skip_header = 322, usecols = (0,1,5), unpack = True)
            data_dict["Event"].extend([f"PS129_{cast.name}"] * len(pressure))
            data_dict["Latitude"].extend([cast.location.lat] * len(pressure))
            data_dict["Longitude"].extend([cast.location.lon] * len(pressure))
            data_dict["Press [dbar]"].extend(pressure)
            data_dict["Temp [°C]"].extend(in_situ_temperature)
            data_dict["Sal"].extend(practical_salinity)
            data_dict["Date/Time"].extend([cast.date] * len(pressure))
            data_dict["Expedition"].extend(["PS129"] * len(pressure))
            
            #add new attributes
            SA = gsw.SA_from_SP(SP = practical_salinity, p = pressure, lon = cast.location.lon, lat = cast.location.lat)
            data_dict["Absolute Salinity"].extend(SA)
            data_dict["Conservative Temperature"].extend(gsw.CT_from_t(SA = SA, t = in_situ_temperature, p = pressure))
            data_dict["Depth water [m]"].extend(np.abs(gsw.z_from_p(p = pressure, lat = cast.location.lat)))
        
        except ValueError as e:
            print("Error at ",cast.name, cast.location)
            print(e)
            continue
        #print(cast.name)     
    #print(data_dict)
    #for (k,v) in data_dict.items():
    #    print(k,len(v))
    return pd.DataFrame(data = data_dict)

PS129_CTDs = get_PS129_CTD_data()
CTDs = CTDs.merge(PS129_CTDs, on = columns, how = "outer")
#Drop all rows that are not even close to the moorings
CTDs.drop(CTDs[CTDs.Latitude < -64.5].index, inplace = True)
CTDs.drop(CTDs[CTDs.Latitude > -63].index, inplace = True)
CTDs.drop(CTDs[CTDs.Longitude < -54].index, inplace = True)
CTDs.drop(CTDs[CTDs.Longitude > -47].index, inplace = True)

CTDs.drop(CTDs[
    m*CTDs.Longitude-b+shift < CTDs.Latitude
].index, inplace = True)

CTDs.drop(CTDs[
    m*CTDs.Longitude-b-shift > CTDs.Latitude 
].index, inplace = True)

CTDs.reset_index(inplace = True, drop=True)

def find_line_number(filename, target_string):
    with open(filename, 'r') as file:
        for line_number, line in enumerate(file):
            if line.startswith(target_string):
                return line_number
    return None

data_paths = helper.IO.get_filepaths_from_directory(directory = "/media/sf_VM_Folder/data/CTD", inclusive = ".tab", exclusive = ())
#for p in data_paths: print(p)

for i, path in enumerate(data_paths):
    target_string = "Event\tDate/Time\tLatitude\tLongitude\tElevation [m]"
    skiprows = find_line_number(path, target_string)
    if skiprows == None:
        target_string = "Event	Type	Date/Time	Longitude"
        skiprows = find_line_number(path, target_string)

    data = pd.read_csv(path, delimiter = "\t", skiprows= skiprows)
    
    data = data[data.columns.intersection(columns)]
    
    #Drop all rows that are not even close to the moorings
    data.drop(data[data.Latitude < -64.5].index, inplace = True)
    data.drop(data[data.Latitude > -63].index, inplace = True)
    data.drop(data[data.Longitude < -54].index, inplace = True)
    data.drop(data[data.Longitude > -47].index, inplace = True)

    data.drop(data[
        m*data.Longitude-b+shift < data.Latitude
    ].index, inplace = True)
    
    data.drop(data[
        m*data.Longitude-b-shift > data.Latitude 
    ].index, inplace = True)
            
    data.reset_index(inplace = True, drop=True)
    
    data['Date/Time'] =  pd.to_datetime(data['Date/Time'])#, format='%d%b%Y:%H:%M:%S.%f')
    data['Event'] = data['Event'].astype('category')
    
    try:
        if data['Event'].iloc[0][4] != "/": 
            current_expedition = data['Event'].iloc[0][0:5]
        else:
            current_expedition = data['Event'].iloc[0][0:4]
        print(current_expedition, path)
        data['Expedition'] = current_expedition
        CTDs = CTDs.merge(data, on = columns, how = "outer")  
        assert not data.empty
    
    except IndexError as e:
        print(f"Error loading {path}")
        assert data.empty
        continue

        
CTDs['Event'] = CTDs['Event'].astype('category')    
CTDs['Expedition'] = CTDs['Expedition'].astype('category')  

CTDs_grouped = CTDs.groupby("Event")
events = CTDs_grouped.groups.keys()


# # Show Overview

# In[12]:


CTDs_grouped = CTDs.groupby("Event")
events = CTDs_grouped.groups.keys()



# Select Example Mooring
mooring = list_of_moorings[1]
print(f"{mooring.location = }, {mooring.time_delta = }")

measurement_depth  =  "1513"

import matplotlib.dates as mdates

cv = mooring[str(measurement_depth)].to_numpy()
time = mooring["time"].to_numpy()

TIME_BANDWIDTH_PRODUCT = 10

complex_velocity = mooring[measurement_depth]
complex_velocity_array = helper.Data.cut_trailing_nans(
    complex_velocity.to_numpy()
)
freq, velocity_spectrum = spectra.total_multitaper(
    complex_velocity_array, dt=1 / 12, P=TIME_BANDWIDTH_PRODUCT
)
assert not np.any(np.isnan(velocity_spectrum))


def integrate_psd_interval(freq,psd,a = 0,b = 99):
    """
    Integration between a und b using the trapezoidal integration method
    """
    lower = np.argmin(np.abs(freq-a)).astype(int)
    upper = np.argmin(np.abs(freq-b)).astype(int)
    return np.trapz(y = psd[lower:upper], x = freq[lower:upper]) 

print(f"integrated spectrum =\t{integrate_psd_interval(freq,velocity_spectrum):.6e}")
print(f"var(uv) =\t\t{np.var(complex_velocity_array):.3e}")


# To get the interpretation of horizontal kinetic energy per frequency range, we have to divide the spectrum by 2
kinetic_psd = velocity_spectrum/2


u = np.real(complex_velocity_array)
v = np.imag(complex_velocity_array)

# mean removal
u -= np.mean(u)
v -= np.mean(v)
print(f"integrated spectrum =\t{integrate_psd_interval(freq,kinetic_psd):.3e}")
print(f"1/2 |u+iv|^2 =\t{0.5*np.mean(np.abs(complex_velocity_array)**2):.3e}")


# ### Adding tidal frequencies
TIDAL_CONSTITUENTS = {"M2": 12.42, "S2": 12.00, "N2":12.66, "K2": 11.96, "K1": 23.93, "O1": 25.82, "P1":24.07, "Q1":26.87}
print(TIDAL_CONSTITUENTS)
# sort the constituents in descending order after their tidal periods
# This is important for finding overlapping tidal peaks
TIDAL_CONSTITUENTS = dict(sorted(TIDAL_CONSTITUENTS.items(), key=lambda x: x[1], reverse=True))


# ### Adding Coriolis frequency
coriolis_frequency_in_cpd = helper.Constants.get_coriolis_frequency(
    mooring.location.lat, unit="cpd", absolute=True
)
print(f"{coriolis_frequency_in_cpd = :.1f} cpd")
coriolis_frequency_in_rads = helper.Constants.get_coriolis_frequency(
    mooring.location.lat, unit="rad/s", absolute=True
)
print(f"{coriolis_frequency_in_rads = :.1e} rad/s")


# ### Adding Buoyancy frequency N
# get instrument depth in units of meter above the sea floor
mab_of_measurement = int(max_depth_dict[mooring.location.lon]) - int(measurement_depth) 

if mab_of_measurement < 0 and mab_of_measurement > -2:
    print(f"Instrument depth was corrected from {mab_of_measurement} to 0 mab.")
    mab_of_measurement = 0

# get a typical N value, derived from CTD profiles
column_name = f"({mooring.location.lat:.2f},{mooring.location.lon:.2f})"    
avrg_N_in_rads = N_table.loc[
    N_table["mab"] == mab_of_measurement, column_name
].item()

avrg_N_in_cpd = avrg_N_in_rads /(2*np.pi) *86400
print(f"{avrg_N_in_cpd = :.1f} cpd")


# ## Energy Conversion $E_\text{kin} \rightarrow E_\text{tot}$

# From Friederike Pollmanns thesis, section 5.2 we know

# \begin{equation} 
# \mathcal{E} = 2 \cdot \frac{N² - f²}{N² - \omega²}\frac{\omega²}{\omega² + f²} \mathcal{U}
# \end{equation}

# In[38]:

# Coriolis frequency in rad/s
f = helper.Constants.get_coriolis_frequency(
    mooring.location.lat, unit="rad/s", absolute=True
)
print(f"{f=:.1e}") #1.31e-4

# buoyancy frequency in rad/s
N = avrg_N_in_rads #5.06e-4
print(f"{N=:.1e}")


def kinetic_to_total_energy(f, N, omega):
    conversion = (
        2*(N**2 - f**2)
        / (N**2 - omega**2)
        * (omega**2)
        / (omega**2 + f**2)
    )
    return conversion


fN_freq, fN_kinetic_psd = helper.Constants.cut_to_f_N_interval(
    freq,
    kinetic_psd, 
    f = coriolis_frequency_in_cpd,
    N = avrg_N_in_cpd
)
factor = kinetic_to_total_energy(
    f = coriolis_frequency_in_cpd,
    N = avrg_N_in_cpd,
    omega = fN_freq
)
assert np.all(factor>=0.999)
total_psd = factor*fN_kinetic_psd


# ### Extension up to $N$

# We extend the spectra up to $N$ with $A\omega^{s}$, with a constant spectral slope $s<0$. For comparison, the GM model assumes a spectral slope of $-2$.

# Because the conversion factor is not constant close to $N$, we fit the slope to the unaltered horizontal kinetic energy spectrum and "glue" it then to total energy spectrum, which from theory should have the same spectral slope. 

# In[44]:


last_freq = fN_freq[-1]
last_value = np.mean(total_psd[-50:])
extent_freq = np.linspace(last_freq, avrg_N_in_cpd, 20)

def extension(omega, A, s = -2):
    return A*omega**s

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

def GM_func(x, c):
    return x**-2 * c
gm_popt, gm_pcov  = curve_fit(GM_func, x, y, p0 = 1e4)


# In[48]:


integral_f_wmax = integrate_psd_interval(fN_freq, fN_kinetic_psd)
print(f"Measured E_kin:\tIntegral f - w_max = {integral_f_wmax :.1e}")
integral_wmax_N = integrate_psd_interval(extent_freq, func_powerlaw(extent_freq, *popt))
print(
    f"fitted slope:\tIntegral w_max - N = {integral_wmax_N:.1e}, Ratio to E_kin= {(integral_wmax_N)/(integral_f_wmax):.1%}"
)
gm_integral_wmax_N = integrate_psd_interval(extent_freq, GM_func(extent_freq, gm_popt))
print(
    f"GM slope:\tIntegral w_max - N = {gm_integral_wmax_N:.1e}, Ratio to E_kin = {(gm_integral_wmax_N)/(integral_f_wmax):.1%}"
)


# In[50]:


#fit slope height for the total kinetic spectrum, given the calculated slope of the kinetic energy spectrum
x = fN_freq[start_index:]
y = total_psd[start_index:]
def func_powerlaw_total_psd(x, c):
    return x**slope * c
slope_height, _pcov  = curve_fit(func_powerlaw_total_psd, x, y, p0 = np.asarray([1e4]))


# In[51]:


integral_f_wmax = integrate_psd_interval(fN_freq, total_psd)
print(f"Measured E_kin:\tIntegral f - w_max = {integral_f_wmax :.1e}")
integral_wmax_N = integrate_psd_interval(extent_freq, func_powerlaw_total_psd(extent_freq, slope_height))
print(
    f"fitted slope:\tIntegral w_max - N = {integral_wmax_N:.1e}, Ratio to total E= {(integral_wmax_N)/(integral_f_wmax):.1%}"
)


# Put everything together in a complicated figure

#fig, ax = plt.subplots(nrows = 2, figsize=(TWO_COLUMN_WIDTH*cm, 0.8*TWO_COLUMN_WIDTH*cm), height_ratios = (1,2))
fig, ax = plt.subplots(2, figsize = (8,6), height_ratios = (1,2))


cv = mooring[str(1513)].to_numpy()
time = mooring["time"].to_numpy()

#print(np.shape(cv))
time = time[750:1150]
cv = cv[750:1150]

def get_days_between(datePast, dateFuture):
    difference = dateFuture - datePast
    return difference.astype('timedelta64[s]').astype(np.int32) / datetime.timedelta(days=1).total_seconds()

days = [get_days_between(time[0], t) for t in time]

ax[0].plot(days, np.real(cv),"k", label = 'West-East')
ax[0].plot(days, np.imag(cv),"gray", label = 'North-South')
ax[0].legend()

ax[0].set_xlabel(f"Days since {pd.to_datetime(time[0]).strftime('%B %d, %Y')}")
#axis[1].set_ylabel("$v$ velocity / (m/s)")

#plt.xlabel("common X")
ax[0].set_ylabel('Velocity / (m/s)')

fill_label_flag = False

#ax.loglog(24*3600*omega/(2*np.pi), 2*np.pi*(K_omg+P_omg)/3600/24, lw = 3, color = "tab:blue", alpha = 0.5, label = "E = KE + PE")
ax[1].loglog(freq[-4000:],kinetic_psd[-4000:], c = "gray", label = "horizontal\nKE spectrum")
ax[1].loglog(fN_freq,total_psd, c = "k", label = "total\nenergy spectrum")

ax[1].plot(extent_freq, func_powerlaw_total_psd(extent_freq,slope_height), c = "k", ls =':', lw = 4, label = f"extension with\nconstant slope {slope:1.1f}")

# Create an inset figure
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axins = ax[1].inset_axes(bounds = [0.02,0.05,0.35,0.8])
axins.vlines(
    coriolis_frequency_in_cpd,
    1e-6,
    1e1,
    color="r",
    alpha=0.5,
    linestyle="-",
    linewidth=2,
)
axins.text(1.76, 0.7 * 1e-3, "f", color="r", alpha=0.8, size = "x-large")
axins.vlines(
    24 / 12.4, 1e-6, 1e-1, color="tab:red", alpha=0.6, linestyle="--", linewidth=2
)
axins.text(1.95, 0.15 * 1e-3, "M2", color="tab:red", alpha=0.8, size = "x-large")

# Turn ticklabels of insets off
axins.tick_params(
    axis="both",
    which="both",
    labelleft=False,
    labelbottom=False,
    left=False,
    bottom=False,
)

axins.loglog(freq,kinetic_psd, c = "gray", label = "horizontal\nKE spectrum")
axins.loglog(fN_freq,total_psd, c = "k", label = "total\nenergy spectrum")


axins.set_xbound(1.7, 2.3)
axins.set_ybound(1e-4, 0.03)
axins.set_title(r"Zoom onto f/M2", loc = "left")
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ylims = ax[1].get_ylim()
ax[1].set_ylim(ylims)

ax2 = helper.Plot.add_tidal_and_f_ticks(ax[1], coriolis_frequency_in_cpd, upper_bound = 1e1)
ax[1].vlines(
    avrg_N_in_cpd, 1e-6, 1e-1, color="r", alpha=0.6, linestyle="-", linewidth=2
)
ax[1].text(12, 1e-5, "N", color="r", alpha=1, size = "x-large")

ax[1].vlines(
    coriolis_frequency_in_cpd, 1e-6, 1e-1, color="r", alpha=0.6, linestyle="-", linewidth=2
)
ax[1].text(1.5, 1e-5, "f", color="r", alpha=1, size = "x-large")


ax[1].legend(loc = "upper right")
fig.suptitle(f"Measurement at {mooring.location.pretty_print()}, Depth of {measurement_depth}$\,$m, {mab_of_measurement}$\,$m above ground")

ax[1].set_xlabel("Frequency (cycles/day)")
ax[1].set_ylabel("PSD (m$^2$/s$^2$ days)")
ax[1].set_xlim(1e-1,20);
#ax.set_xlim(5e-1,17)


""""""
#connection lines
#y0 = ax.get_ylim()[0]
#ax.plot([coriolis_frequency_in_cpd, 0.039], [y0+0.5e-6,1.65e-5], "k--", lw = 0.5) #f
#ax.plot([24 / 12.4, 0.081], [y0+1.5e-6,1.65e-5],"k--", lw =0.5) #M2       



fig.tight_layout()

plt.show()


