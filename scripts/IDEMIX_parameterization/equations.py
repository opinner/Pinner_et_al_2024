import numpy as np


def get_dissipation_rate(coriolis_frequency, buoyancy_frequency, energy_level):
    """
    eq. 18 from Olbers & Eden, 2013
    """
    MIXING_EFFICIENCY = 0.2
    MU_0 = 1 / 3  # value recommended by Pollmann et al., 2017
    m_star = 0.01
    effective_coriolis_frequency = (
            np.abs(coriolis_frequency) *
            np.arccosh(buoyancy_frequency / np.abs(coriolis_frequency))
    )

    dissipation_rate = (
            1 / (1 + MIXING_EFFICIENCY)
            * MU_0
            * effective_coriolis_frequency
            * m_star ** 2
            * energy_level ** 2
            / buoyancy_frequency ** 2
    )

    return dissipation_rate


def get_multiplicative_error_of_dissipation_rate(coriolis_frequency, buoyancy_frequency, error_buoyancy_frequency, energy_level,
                                  error_energy_level):
    """
    calculates the error to the dissipation rate order of magnitude
    dissipation rate is calculated with eq. 18 from Olbers & Eden, 2013
    for derivation of the equations see /docs/energy_level_error_propagation.ipynb
    """

    first_N_summand = (
        1 / (
            np.sqrt(buoyancy_frequency**2 - coriolis_frequency**2)
            * np.log(10)
            * np.arccosh(buoyancy_frequency / np.abs(coriolis_frequency))
        )
    )

    second_N_summand = (
        2
        / (buoyancy_frequency * np.log(10))
    )

    error_due_to_N = (
            error_buoyancy_frequency
            * (first_N_summand - second_N_summand)
    )

    error_due_to_E = (
            error_energy_level
            * 2
            / (energy_level * np.log(10))
    )

    error_to_magnitude = np.sqrt(error_due_to_N ** 2 + error_due_to_E ** 2)
    # return error of the observable instead of error to the magnitude
    return 10**error_to_magnitude


def get_additive_error_of_dissipation_rate(coriolis_frequency, buoyancy_frequency, error_buoyancy_frequency, energy_level,
                                  error_energy_level):
    """
    eq. 18 from Olbers & Eden, 2013
    """
    MIXING_EFFICIENCY = 0.2
    MU_0 = 1 / 3  # value recommeded by Pollmann et al 2017
    m_star = 0.01

    arccosh_N_derivative = (
            1 / (
            np.sqrt(buoyancy_frequency / np.abs(coriolis_frequency) - 1)
            * np.sqrt(buoyancy_frequency / np.abs(coriolis_frequency) + 1)
    )
            * m_star ** 2
            * energy_level ** 2
            / buoyancy_frequency ** 2
    )

    exponent_N_derivative = (
            np.arccosh(buoyancy_frequency / np.abs(coriolis_frequency))
            * -2 * m_star ** 2
            * energy_level ** 2
            / buoyancy_frequency ** 3
    )

    error_due_to_N = (
            1 / (1 + MIXING_EFFICIENCY)
            * MU_0 * np.abs(coriolis_frequency)
            * (arccosh_N_derivative + exponent_N_derivative)
            * error_buoyancy_frequency
    )

    effective_coriolis_frequency = (
            np.abs(coriolis_frequency) *
            np.arccosh(buoyancy_frequency / np.abs(coriolis_frequency))
    )

    error_due_to_E = (
            1 / (1 + MIXING_EFFICIENCY)
            * MU_0
            * effective_coriolis_frequency
            * 2 * m_star ** 2
            * energy_level
            / (buoyancy_frequency ** 2)
            * error_energy_level
    )

    dissipation_total_error = np.sqrt(error_due_to_N ** 2 + error_due_to_E ** 2)

    return dissipation_total_error
