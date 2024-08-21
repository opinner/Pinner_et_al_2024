


def kinetic_to_total_energy(f, N, omega):
    conversion = (
            2 * (N ** 2 - f ** 2)
            / (N ** 2 - omega ** 2)
            * (omega ** 2)
            / (omega ** 2 + f ** 2)
    )
    return conversion

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
