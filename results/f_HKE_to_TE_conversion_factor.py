import numpy as np
import src.gm_library as gm
import src.helper_functions as helper
import src.spectra as spectra

def calculate_GM_spectrum(f,N,N0,b):
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
    omg = np.logspace(np.log10(1.01*f), np.log10(N), 401) #Where does the 401 come from?

    # horizontal wavenumber
    k = 2*np.pi*np.logspace(-6, -2, 401)

    # mode number
    j = np.arange(1, 100)

    # reshape to allow multiplication into 2D array
    Omg = np.reshape(omg, (omg.size,1))
    K = np.reshape(k, (k.size,1))
    J = np.reshape(j, (1,j.size))

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



# Coriolis frequency in rad/s
f = helper.Constants.get_coriolis_frequency(
    mooring.location.lat, unit="rad/s", absolute=True
)
print(f"{f=:.1e}") #1.31e-4

# buoyancy frequency in rad/s
N = avrg_N_in_rads #5.06e-4
print(f"{N=:.1e}")

# surface-extrapolated buoyancy frequency
N0 = N #5.06e-4

# e-folding scale of N(z)
b = 1.7e3


# Correction factor in cpd
kinetic_to_total_energy_factor = kinetic_to_total_energy(
    f=coriolis_frequency_in_cpd,
    N=avrg_N_in_cpd,
    omega=fN_freq
)


# correction factor cannot be larger than 1, as HKE < HKE + APE
assert np.all(kinetic_to_total_energy_factor > 0.999)

# Compare to Correction factor in rad/s
difference = kinetic_to_total_energy_factor - kinetic_to_total_energy(
    f=coriolis_frequency_in_rads, N=avrg_N_in_rads, omega=fN_freq * (2 * np.pi) / 86400
)
assert np.all(np.abs(difference) < 1e-10)

# PSD and the frequency is in units of cpd
resolved_total_energy_spectrum_between_f_and_N = kinetic_to_total_energy_factor * resolved_HKE_spectrum_between_f_and_N
