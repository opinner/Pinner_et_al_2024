#!/usr/bin/env python
# coding: utf-8

# import scipy.io as sio
import datetime

# import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import src.helper as helper
import src.spectra as spectra

# from src.mooring import Mooring
# from src.location import Location
# from src.ctd_cast import CTDCast

print(plt.rcParams["font.size"])

plt.style.use('./paper.mplstyle')
legend_font_size = 7

import pandas as pd
#warnings.filterwarnings("ignore")  # suppress some warnings about future code changes

import locale

locale.setlocale(locale.LC_ALL, 'en_US.utf8')

ONE_COLUMN_WIDTH = 8.3
TWO_COLUMN_WIDTH = 12
GOLDEN_RATIO = 1.61
cm = 1 / 2.54  # centimeters in inches

# load all 7 moorings as dataframes
list_of_moorings = helper.IO.load_pickle(name="../data/mooring/list_of_moorings.pkl")
# load Stratification information
N_table = pd.read_pickle("../scripts/IDEMIX_parameterization/method_results/N_values.pkl")

data = np.load("../data/mooring/max_depth_dict.npz", allow_pickle=True)
max_depth_dict = data["max_depth_dict"].item()

# Select Example Mooring
mooring = list_of_moorings[1]
measurement_depth = "1513"
print(f"{mooring.location = }, {mooring.time_delta = }")

cv = mooring[str(measurement_depth)].to_numpy()
time = mooring["time"].to_numpy()

TIME_HALF_BANDWIDTH_PRODUCT = 10

complex_velocity = mooring[measurement_depth]
complex_velocity_array = helper.Data.cut_trailing_nans(
    complex_velocity.to_numpy()
)

#print(np.shape(complex_velocity_array))

freq, velocity_spectrum = spectra.total_multitaper(
    complex_velocity_array, dt=1 / 12, P=TIME_HALF_BANDWIDTH_PRODUCT
)
assert not np.any(np.isnan(velocity_spectrum))


def integrate_psd_interval(freq, psd, a=0, b=99):
    """
    Integration between a und b using the trapezoidal integration method
    """
    lower = np.argmin(np.abs(freq - a)).astype(int)
    upper = np.argmin(np.abs(freq - b)).astype(int)
    return np.trapz(y=psd[lower:upper], x=freq[lower:upper])


print(f"integrated spectrum =\t{integrate_psd_interval(freq, velocity_spectrum):.6e}")
print(f"var(uv) =\t\t{np.var(complex_velocity_array):.3e}")

# To get the interpretation of horizontal kinetic energy per frequency range, we have to divide the spectrum by 2
kinetic_psd = velocity_spectrum / 2

u = np.real(complex_velocity_array)
v = np.imag(complex_velocity_array)

# mean removal
u -= np.mean(u)
v -= np.mean(v)
print(f"integrated spectrum =\t{integrate_psd_interval(freq, kinetic_psd):.3e}")
print(f"1/2 |u+iv|^2 =\t{0.5 * np.mean(np.abs(complex_velocity_array) ** 2):.3e}")

# ### Adding tidal frequencies
TIDAL_CONSTITUENTS = {"M2": 12.42, "S2": 12.00, "N2": 12.66, "K2": 11.96, "K1": 23.93, "O1": 25.82, "P1": 24.07,
                      "Q1": 26.87}
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

avrg_N_in_cpd = avrg_N_in_rads / (2 * np.pi) * 86400
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
print(f"{f=:.1e}")  #1.31e-4

# buoyancy frequency in rad/s
N = avrg_N_in_rads  #5.06e-4
print(f"{N=:.1e}")


def kinetic_to_total_energy(f, N, omega):
    conversion = (
            2 * (N ** 2 - f ** 2)
            / (N ** 2 - omega ** 2)
            * (omega ** 2)
            / (omega ** 2 + f ** 2)
    )
    return conversion


fN_freq, fN_kinetic_psd = helper.Constants.cut_to_f_N_interval(
    freq,
    kinetic_psd,
    f=coriolis_frequency_in_cpd,
    N=avrg_N_in_cpd
)
factor = kinetic_to_total_energy(
    f=coriolis_frequency_in_cpd,
    N=avrg_N_in_cpd,
    omega=fN_freq
)
assert np.all(factor >= 0.999)
total_psd = factor * fN_kinetic_psd

# ### Extension up to $N$

# We extend the spectra up to $N$ with $A\omega^{s}$, with a constant spectral slope $s<0$. For comparison, the GM model assumes a spectral slope of $-2$.

# Because the conversion factor is not constant close to $N$, we fit the slope to the unaltered horizontal kinetic energy spectrum and "glue" it then to total energy spectrum, which from theory should have the same spectral slope. 

# In[44]:


last_freq = fN_freq[-1]
last_value = np.mean(total_psd[-50:])
extent_freq = np.linspace(last_freq, avrg_N_in_cpd, 20)


def extension(omega, A, s=-2):
    return A * omega ** s


start_freq = 3.2  #cpd
start_index = np.argmin(np.abs(fN_freq - start_freq))

x = fN_freq[start_index:]
y = fN_kinetic_psd[start_index:]
from scipy.optimize import curve_fit


def func_powerlaw(x, m, c):
    return x ** m * c


popt, pcov = curve_fit(func_powerlaw, x, y, p0=np.asarray([-2, 1e4]))
slope = tuple(popt)[0]
slope_error = np.sqrt(np.diag(pcov))[0]


def GM_func(x, c):
    return x ** -2 * c


gm_popt, gm_pcov = curve_fit(GM_func, x, y, p0=1e4)

# In[48]:


integral_f_wmax = integrate_psd_interval(fN_freq, fN_kinetic_psd)
print(f"Measured E_kin:\tIntegral f - w_max = {integral_f_wmax :.1e}")
integral_wmax_N = integrate_psd_interval(extent_freq, func_powerlaw(extent_freq, *popt))
print(
    f"fitted slope:\tIntegral w_max - N = {integral_wmax_N:.1e}, Ratio to E_kin= {(integral_wmax_N) / (integral_f_wmax):.1%}"
)
gm_integral_wmax_N = integrate_psd_interval(extent_freq, GM_func(extent_freq, gm_popt))
print(
    f"GM slope:\tIntegral w_max - N = {gm_integral_wmax_N:.1e}, Ratio to E_kin = {(gm_integral_wmax_N) / (integral_f_wmax):.1%}"
)

# In[50]:


#fit slope height for the total kinetic spectrum, given the calculated slope of the kinetic energy spectrum
x = fN_freq[start_index:]
y = total_psd[start_index:]


def func_powerlaw_total_psd(x, c):
    return x ** slope * c


slope_height, _pcov = curve_fit(func_powerlaw_total_psd, x, y, p0=np.asarray([1e4]))

# In[51]:


integral_f_wmax = integrate_psd_interval(fN_freq, total_psd)
print(f"Measured E_kin:\tIntegral f - w_max = {integral_f_wmax :.1e}")
integral_wmax_N = integrate_psd_interval(extent_freq, func_powerlaw_total_psd(extent_freq, slope_height))
print(
    f"fitted slope:\tIntegral w_max - N = {integral_wmax_N:.1e}, Ratio to total E= {(integral_wmax_N) / (integral_f_wmax):.1%}"
)

fig, ax = plt.subplots(nrows=2, figsize=(TWO_COLUMN_WIDTH * cm, 0.8 * TWO_COLUMN_WIDTH * cm), height_ratios=(1, 2))

cv = mooring[str(1513)].to_numpy()
time = mooring["time"].to_numpy()

#print(np.shape(cv))
time = time[750:1150]
cv = cv[750:1150]


def get_days_between(datePast, dateFuture):
    difference = dateFuture - datePast
    return difference.astype('timedelta64[s]').astype(np.int32) / datetime.timedelta(days=1).total_seconds()


days = [get_days_between(time[0], t) for t in time]

ax[0].grid()
ax[0].plot(days, np.real(cv), "k", label='West-East')
ax[0].plot(days, np.imag(cv), "tab:purple", label='North-South')
ax[0].legend(fontsize=legend_font_size)

ax[0].set_xlabel(f"Days since {pd.to_datetime(time[0]).strftime('%B %d, %Y')}")
#axis[1].set_ylabel("$v$ velocity / (m/s)")

#plt.xlabel("common X")
ax[0].set_ylabel(r'Velocity (m$\,$s$^{-1}$)')

fill_label_flag = False

#ax.loglog(24*3600*omega/(2*np.pi), 2*np.pi*(K_omg+P_omg)/3600/24, lw = 3, color = "tab:blue", alpha = 0.5, label = "E = KE + PE")
ax[1].loglog(freq[-4000:], kinetic_psd[-4000:], c="tab:blue", label=r"$\mathcal{U}$")
ax[1].loglog(fN_freq, total_psd, c="k", label=r"$\mathcal{E}$")

ax[1].plot(extent_freq, func_powerlaw_total_psd(extent_freq, slope_height), c="k", ls=':', lw=4,
           label=f"extension")  # with\nconstant slope {slope:1.1f}")

# Create an inset figure
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axins = ax[1].inset_axes(bounds=[0.02, 0.05, 0.35, 0.8])
axins.vlines(
    coriolis_frequency_in_cpd,
    1e-6,
    1e1,
    color="tab:red",
    alpha=1,
    linestyle="-",
    linewidth=2,
)
axins.text(1.76, 0.7 * 1e-3, "f", color="tab:red", alpha=1, size="x-large")
axins.vlines(
    24 / 12.4, 1e-6, 1e-1, color="tab:red", alpha=1, linestyle="--", linewidth=2
)
axins.text(1.95, 0.15 * 1e-3, "M2", color="tab:red", alpha=1, size="x-large")

# Turn ticklabels of insets off
axins.tick_params(
    axis="both",
    which="both",
    labelleft=False,
    labelbottom=False,
    left=False,
    bottom=False,
)

axins.loglog(freq, kinetic_psd, c="tab:blue", label=r"$\mathcal{U}$")
axins.loglog(fN_freq, total_psd, c="k", label=r"$\mathcal{E}$")

axins.set_xbound(1.7, 2.3)
axins.set_ybound(1e-4, 0.03)
axins.set_title(r"Zoom onto f/M2", loc="left", fontsize=9, pad=4.0)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ylims = ax[1].get_ylim()
ax[1].set_ylim(ylims)

ax2 = helper.Plot.add_tidal_and_f_ticks(ax[1], coriolis_frequency_in_cpd, upper_bound=1e1)
ax[1].vlines(
    avrg_N_in_cpd, 1e-6, 1e-1, color="tab:red", alpha=1, linestyle="-", linewidth=2
)
ax[1].text(12, 1e-5, "N", color="tab:red", alpha=1, size="x-large")

ax[1].vlines(
    coriolis_frequency_in_cpd, 1e-6, 1e-1, color="tab:red", alpha=1, linestyle="-", linewidth=2
)
ax[1].text(1.5, 1e-5, "f", color="tab:red", alpha=1, size="x-large")

ax[1].legend(loc="upper right", fontsize=legend_font_size)
# fig.suptitle(f"Measurement at {mooring.location.pretty_print()}, Depth of {measurement_depth}$\,$m, {mab_of_measurement}$\,$m above ground")

ax[1].set_xlabel("Frequency (cycles per day)")
ax[1].set_ylabel(r"Energy density (m$^2\,$s$^{-2}$ days)")
ax[1].set_xlim(1e-1, 20);
#ax.set_xlim(5e-1,17)

#connection lines
#y0 = ax.get_ylim()[0]
#ax.plot([coriolis_frequency_in_cpd, 0.039], [y0+0.5e-6,1.65e-5], "k--", lw = 0.5) #f
#ax.plot([24 / 12.4, 0.081], [y0+1.5e-6,1.65e-5],"k--", lw =0.5) #M2

fig.tight_layout()

#fig.savefig("./TimeSeries_plus_Spectrum.pdf")
plt.show()
