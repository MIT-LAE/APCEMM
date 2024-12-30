import numpy as np


######################### WATER THERMODYNAMIC FUNCTIONS #########################
def compute_p_sat_liq(T_K):
    # Calculates the liquid saturation pressure "p_sat_liq" (in units of Pa).
    # This equation can be found in Schumann 2012.
    #
    # The inputs are as follows:
    # - T_K is the absolute temperature in Kelvin

    a = -6096.9385
    b = 16.635794
    c = -0.02711193
    d = 1.673952e-5
    e = 2.433502

    return 100 * np.exp(a / T_K + b + c * T_K + d * T_K * T_K + e * np.log(T_K))


def compute_p_sat_ice(T_K):
    # Calculates the ice saturation pressure "p_sat_ice" (in units of Pa)
    # The equation can be found in Schumann 2012.
    #
    # The inputs are as follows
    # - T_K is the absolute temperature in Kelvin

    a = -6024.5282
    b = 24.7219
    c = 0.010613868
    d = -1.3198825e-5
    e = -0.49382577

    return 100 * np.exp(a / T_K + b + c * T_K + d * T_K * T_K + e * np.log(T_K))


def convert_RH_to_RHi(T_K, RH):
    return RH * compute_p_sat_liq(T_K) / compute_p_sat_ice(T_K)


def convert_RHi_to_RH(T_K, RHi):
    return RHi * compute_p_sat_ice(T_K) / compute_p_sat_liq(T_K)


def convert_RH_to_SH(T_K, p_Pa, RH_PC):
    # Implements the equations from https://earthscience.stackexchange.com/a/2361
    # The output specific humidity is in ratio format.
    # The input RH must be in %.

    RH_ratio = RH_PC / 100.0
    R_d = 287  # J kg^-1 K^-1 from https://en.wikipedia.org/wiki/Gas_constant#:~:text=Specific%20gas%20constant,-Rspecific&text=for%20dry%20air%20of%2028.964917%20g%2Fmol.
    R_v = 461.5  # J kg^-1 K^-1 from https://en.wikipedia.org/wiki/Water_vapor

    e_s_Pa = compute_p_sat_liq(T_K)  # Saturation pressure of water vapour
    w_s = (
        (R_d / R_v) * e_s_Pa / (p_Pa - e_s_Pa)
    )  # Mass mixing ratio of water vapor to dry air at equilibrium (dimensionless)

    w = RH_ratio * w_s  # Mass mixing ratio of water vapor to dry air

    return w / (w + 1)


def convert_RHi_to_SH(T_K, p_Pa, RHi_PC):
    # Implements the equations from https://earthscience.stackexchange.com/a/2361
    # The output specific humidity is in ratio format.
    # The input RHi must be in %.

    RH_PC = convert_RHi_to_RH(T_K, RHi_PC)

    return convert_RH_to_SH(T_K, p_Pa, RH_PC)


def convert_SH_to_RH(T_K, p_Pa, SH):
    # Implements the equations from https://earthscience.stackexchange.com/a/2361
    # The specific humidity must be in ratio format. The output RH is in %

    R_d = 287  # J kg^-1 K^-1 from https://en.wikipedia.org/wiki/Gas_constant#:~:text=Specific%20gas%20constant,-Rspecific&text=for%20dry%20air%20of%2028.964917%20g%2Fmol.
    R_v = 461.5  # J kg^-1 K^-1 from https://en.wikipedia.org/wiki/Water_vapor

    e_s_Pa = compute_p_sat_liq(T_K)  # Saturation pressure of water vapour
    w_s = (
        (R_d / R_v) * e_s_Pa / (p_Pa - e_s_Pa)
    )  # Mass mixing ratio of water vapor to dry air at equilibrium (dimensionless)

    w = SH / (1 - SH)

    return w / w_s * 100


def convert_SH_to_RHi(T_K, p_Pa, SH):
    # Implements the equations from https://earthscience.stackexchange.com/a/2361
    # The specific humidity must be in ratio format. The output RHi is in %

    RH_PC = convert_SH_to_RH(T_K, p_Pa, SH)

    return convert_RH_to_RHi(T_K, RH_PC)
