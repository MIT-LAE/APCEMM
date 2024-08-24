import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from src.ISA import *
from src.moist import *
from src.thermo import *

def met_clean_slate(desired_altitudes, desired_RHis = None, RHi_background = 0., shear = 0.002, T_offset_K = 0):
    # Creates a clean slate xr.Dataset object with the weather parameters required by APCEMM at
    # each altitude. The vertical resolution is determined by the "desired_altitudes" input.
    #
    # The met created by this function is time-invariant and is only defined for 24 hours.
    #
    # This function uses the ISA. However, the temperature offset ("T_offset_K") does not affect the
    # ISA calculations (which assume zero temperature offset), except for the temperature and density.
    #
    # The following weather parameters can be changed:
    #   - RHi: Must be in %. Set by the use of a vertical profile (set in "desired_RHis") or by the
    #          use of a constant RHi (set in "RHi_background"; only used when "desired_RHis" is "None").
    #   - Shear: In units of 1/s. A uniform shear field is assumed. (Set by the "shear" input.)
    #   - Temperature offset: In units of K. A linear shift in the ISA temperature, can be either
    #         positive or negative.
    #
    # Returns a xr.Dataset object.

    desired_altitudes = desired_altitudes / 1e3 # m to km

    # Time vector, just in case it is needed
    time_ip = np.concatenate((np.linspace(15, 23, 9), np.linspace(0, 14, 15))) # hours

    # Initialise the met variables
    len_alt = len(desired_altitudes)
    pressure_ip = np.zeros(len_alt) # hPa
    temperature_ip = np.zeros(len_alt) # Degrees Kelvin
    relative_humidity_ip = np.zeros(len_alt) # %
    relative_humidity_ice_ip = np.zeros(len_alt) # %
    shear_ip = np.zeros(len_alt) # s^-1

    # Initialise the dataset with bad weather data
    met_data_temp = xr.Dataset(
        data_vars=dict(
            pressure=(["altitude"], pressure_ip, {'units':'hPa'}),
            temperature=(["altitude"], temperature_ip, {'units':'K'}),
            relative_humidity=(["altitude"], relative_humidity_ip, {'units':'Pct'}),
            relative_humidity_ice=(["altitude"], relative_humidity_ice_ip, {'units':'Pct'}),
            shear=(["altitude"], shear_ip, {'units':'s^-1'}),
        ),
        coords=dict(
            altitude=(["altitude"], desired_altitudes, {'units':'km'}),
            time=(["time"], time_ip, {'units':'h'}),         
        ),
    )

    if desired_RHis is None:
        met_data_ISA = met_from_ISA(
            met = met_data_temp,
            RHi_background = RHi_background,
            T_offset_K = T_offset_K,
            shear = shear
            )
    else:
        RHi_background = np.min(desired_RHis)
        met_data_ISA = met_from_ISA(
            met = met_data_temp,
            RHi_background = RHi_background,
            T_offset_K = T_offset_K,
            shear = shear
            )
        RHi = met_data_ISA.relative_humidity_ice
        RH = met_data_ISA.relative_humidity
        Ts_K = met_data_ISA.temperature.to_numpy()
        desired_RHs = convert_RHi_to_RH(Ts_K, desired_RHis)
        
        # Enforce the desired RHis to each altitude
        for i in range(len(desired_altitudes)):
            current_alt = desired_altitudes[i]

            mask = RHi.altitude != current_alt
            RHi = RHi.where(mask, desired_RHis[i]) 
            RH = RH.where(mask, desired_RHs[i])

        met_data_ISA["relative_humidity_ice"] = RHi
        met_data_ISA["relative_humidity"] = RH

    return met_data_ISA

def truncate_met(met, truncation_alt_km = 20.):
    # Truncates the APCEMM met data to remain beneath the tropopause.

    altitudes = met.altitude
    truncation_idx = len(altitudes)

    # Find the index at which the altitude surpasses the truncation threshold
    for i in range(len(altitudes)):
        if (altitudes[i] > truncation_alt_km):
            truncation_idx = i
            break

    return met.isel(altitude = slice(0, truncation_idx))

def met_from_ISA(met, RHi_background = 0., shear = 0.002, T_offset_K = 0):
    # Clears the existing temperature values of an input xr.Dataset ("met"), and
    # enforces ISA conditions at each point.
    #
    # The met created by this function is time-invariant and is only defined for 24 hours.
    #
    # The temperature offset ("T_offset_K") does not affect the ISA calculations (which 
    # assume zero temperature offset), except for the temperature itself.
    #
    # The following weather parameters can be changed:
    #   - Background RHi: Must be in %. Set by the use of "RHi_background".
    #   - Shear: In units of 1/s. A uniform shear field is assumed. (Set by the "shear" input.)
    #   - Temperature offset: In units of K. A shift in the ISA temperature, can be either
    #         positive or negative. Does not affect the other ISA parameters.
    #
    # Returns a xr.Dataset object.

    altitudes = met.altitude

    T = met["temperature"]
    RH = met["relative_humidity"]
    p = met["pressure"]

    # Enforce the ISA conditions at each altitude
    for i in range(len(altitudes)):
        h_current = altitudes[i] * 1000 # km to m
        T_current = compute_T_ISA(h_current, T_offset_K)
        p_current = compute_p_ISA(h_current) / 100. # Pa to hPa

        mask = met.altitude != altitudes[i]
        T = T.where(mask, T_current)
        RH = RH.where(mask, convert_RHi_to_RH(T_current, RHi_background))
        p = p.where(mask, p_current)

    met["temperature"] = T
    met["relative_humidity"] = RH
    met["pressure"] = p

    # Set the RHi of the met data to the RHi default value
    RHi_met = met["relative_humidity_ice"]
    RHi_met = RHi_met.where(False, RHi_background)
    met["relative_humidity_ice"] = RHi_met

    # Set the shear of the met data to the default value
    shear_met = met["shear"]
    shear_met = shear_met.where(False, shear)
    met["shear"] = shear_met

    return met

def create_and_save_met(alts_m, RHis_PC, T_offset_K = 0, shear = 0.002, prefix = "00"):
    met = met_clean_slate(
        desired_altitudes = alts_m,
        desired_RHis = RHis_PC,
        T_offset_K = T_offset_K,
        shear = shear
        )
    
    Path("outputs/").mkdir(parents=True, exist_ok=True)
    met.to_netcdf('outputs/' + prefix + '-met.nc', engine="netcdf4")
    return met

def make_idealised_met(centre_altitude_m, moist_depth_m, RHi_background_PC = 0., RHi_peak_PC = 150, 
                grad_RHi_PC_per_m = None, alt_resolution_m = 100, upper_alt_m = 11001, shear = 2e-3,
                T_offset_K = 0, filename_prefix = "default"):

    # Creates an idealised met file with a trapezium moist layer.
    #
    # The altitude is sampled at intervals with spacing of alt_resolution_m, up until 
    # upper_alt_m (not inclusive of this)
    #
    # The depth of the moist region indicates where the relative humidity is above background
    # levels.
    #
    # The generated trapezium is defined from the sides inwards, meaning that if it is too
    # narrow, the peak RHi will not be reached anywhere in the moist region.
    #
    # The met created by this function is time-invariant and is only defined for 24 hours.
    
    altitudes, RHis = trapezium_moist_layer(centre_altitude_m, moist_depth_m, RHi_background_PC, 
                                            RHi_peak_PC, grad_RHi_PC_per_m, alt_resolution_m, 
                                            upper_alt_m)
    
    met = create_and_save_met(alts_m=altitudes,RHis_PC=RHis,T_offset_K=T_offset_K,shear=shear,prefix=filename_prefix)

    return met
