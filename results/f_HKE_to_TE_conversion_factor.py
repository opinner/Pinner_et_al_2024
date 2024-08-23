import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import src.gm_library as gm
import src.helper as helper
import src.spectra as spectra
import src.plots as plots

plt.style.use('./thesis.mplstyle')


def calculate_GM_spectrum(f, N, N0, b):
    """
    Here happens the actual calculation of the GM81 spectrum.
    GM equations are copied from https://github.com/joernc/GM81
    under the GNU General Public License.
    """
    # This module implements the empirical spectrum of internal waves developed by
    # Garrett and Munk, in the incarnation presented in Munk's chapter in Evolution
    # of Physical Oceanography, which can be downloaded here:
    # http://ocw.mit.edu/resources/res-12-000-evolution-of-physical-oceanography-spring-2007/part-2/wunsch_chapter9.pdf
    # The variable names follow Munk's notation.

    # frequency
    omg = np.logspace(np.log10(1.01 * f), np.log10(N), 401)  # Where does the 401 come from?

    # horizontal wavenumber
    k = 2 * np.pi * np.logspace(-6, -2, 401)

    # mode number
    j = np.arange(1, 100)

    # reshape to allow multiplication into 2D array
    Omg = np.reshape(omg, (omg.size, 1))
    K = np.reshape(k, (k.size, 1))
    J = np.reshape(j, (1, j.size))

    # frequency spectra (KE and PE)
    K_omg_j = gm.K_omg_j(Omg, J, f, N, N0, b)
    P_omg_j = gm.P_omg_j(Omg, J, f, N, N0, b)

    # wavenumber spectra (KE and PE)
    K_k_j = gm.K_k_j(K, J, f, N, N0, b)
    P_k_j = gm.P_k_j(K, J, f, N, N0, b)

    # sum over modes (j refers to modes)
    K_omg = np.sum(K_omg_j, axis=1)
    P_omg = np.sum(P_omg_j, axis=1)
    K_k = np.sum(K_k_j, axis=1)
    P_k = np.sum(P_k_j, axis=1)

    return omg, K_omg, P_omg


def kinetic_to_total_energy(f, N, omega):
    conversion = (
            2 * (N ** 2 - f ** 2)
            / (N ** 2 - omega ** 2)
            * (omega ** 2)
            / (omega ** 2 + f ** 2)
    )
    return conversion


def main():
    ###############################
    # Observational Data
    ###############################

    # load all 7 moorings as dataframes
    list_of_moorings = helper.IO.load_pickle(name="../data/mooring/list_of_moorings.pkl")
    # load Stratification information
    N_table = pd.read_pickle("../scripts/IDEMIX_parametrization/method_results/N_values.pkl")

    data = np.load("../data/mooring/max_depth_dict.npz", allow_pickle=True)
    max_depth_dict = data["max_depth_dict"].item()

    # Select Example Mooring
    mooring = list_of_moorings[1]
    measurement_depth = "1513"
    print(f"{mooring.location = }, {mooring.time_delta = }")

    # cv = mooring[str(measurement_depth)].to_numpy()
    # time = mooring["time"].to_numpy()

    coriolis_frequency_in_cpd = helper.Constants.get_coriolis_frequency(
        mooring.location.lat, unit="cpd", absolute=True
    )
    print(f"{coriolis_frequency_in_cpd = :.1f} cpd")
    coriolis_frequency_in_rads = helper.Constants.get_coriolis_frequency(
        mooring.location.lat, unit="rad/s", absolute=True
    )
    print(f"{coriolis_frequency_in_rads = :.1e} rad/s")

    ## Adding Buoyancy frequency N

    # get instrument depth in units of meter above the sea floor
    mab_of_measurement = int(max_depth_dict[mooring.location.lon]) - int(measurement_depth)

    if 0 > mab_of_measurement > -2:
        print(f"Instrument depth was corrected from {mab_of_measurement} to 0 mab.")
        mab_of_measurement = 0

    # get a typical N value, derived from CTD profiles
    column_name = f"({mooring.location.lat:.2f},{mooring.location.lon:.2f})"
    avrg_N_in_rads = N_table.loc[
        N_table["mab"] == mab_of_measurement, column_name
    ].item()

    avrg_N_in_cpd = avrg_N_in_rads / (2 * np.pi) * 86400
    print(f"{avrg_N_in_cpd = :.1f} cpd")

    TIME_BANDWIDTH_PRODUCT = 10

    complex_velocity = mooring[measurement_depth]
    complex_velocity_array = helper.Data.cut_trailing_nans(
        complex_velocity.to_numpy()
    )
    fN_freq, velocity_spectrum = spectra.total_multitaper(
        complex_velocity_array, dt=1 / 12, P=TIME_BANDWIDTH_PRODUCT
    )
    assert not np.any(np.isnan(velocity_spectrum))

    ###############################
    # Garrett Munk model
    ###############################
    # Coriolis frequency in rad/s
    f = helper.Constants.get_coriolis_frequency(
        mooring.location.lat, unit="rad/s", absolute=True
    )
    print(f"{f=:.1e}")  #1.31e-4
    # buoyancy frequency in rad/s
    N = avrg_N_in_rads  #5.06e-4
    print(f"{N=:.1e}")
    # surface-extrapolated buoyancy frequency
    N0 = N  #5.06e-4
    # e-folding scale of N(z)
    b = 1.7e3

    omega, K_omg, P_omg = calculate_GM_spectrum(f, N, N0, b)

    ###############################
    # Figure
    ###############################
    fig, ax = plt.subplots(1, figsize=plots.set_size("thesis", fraction=0.9), layout="tight")
    factor = kinetic_to_total_energy(f, N, omega)
    E_omg = factor * K_omg
    ax.plot(24 * 3600 * omega[:-40] / (2 * np.pi), factor[:-40], color="tab:blue", lw=2, label="from wave theory")
    ax.plot(24 * 3600 * omega / (2 * np.pi), (K_omg + P_omg) / K_omg, zorder=5, lw=2, color="k",
            label="from Garrett-Munk model")
    ax.axvline(6, color="k", alpha=1, linestyle="--", linewidth=2)
    ax.text(5.8, 3.2, "data\nresolution", color="k", ha="right")
    ax.axvline(avrg_N_in_cpd, color="tab:red", alpha=0.6, linestyle="--", linewidth=2)
    ax.text(15.7, 2.5, "N", color="tab:red", alpha=1, size="large")
    ax.axvline(coriolis_frequency_in_cpd, color="r", alpha=0.6, linestyle="--", linewidth=2)
    ax.text(1.3, 2.5, "f", color="tab:red", alpha=1, size="large")

    ax.legend(loc="upper left", framealpha=0.9)
    ax.set_xlabel("Frequency (cycles per day)")
    ax.set_ylabel("HKE to TE factor")
    fig.savefig(f"./HKE_TE_factor.pdf")

if __name__ == "__main__":
    main()
    plt.show()
