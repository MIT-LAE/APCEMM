import xarray as xr
import numpy as np

from pathlib import Path
from src.ISA import compute_T_ISA, compute_p_ISA
from src.moist import trapezium_moist_layer


def met_clean_slate(
    desired_altitudes_m,
    desired_RHis_PC=None,
    RHi_background_PC=0.0,
    shear_over_s=0.002,
    T_offset_K=0,
    w_Pa_per_s=0,
):
    # Creates a clean slate xr.Dataset object with the weather parameters required by APCEMM at
    # each altitude. The vertical resolution is determined by the "desired_altitudes_m" input.
    #
    # The met created by this function is time-invariant and is only defined for 24 hours.
    #
    # This function uses the ISA. However, the temperature offset ("T_offset_K") does not affect the
    # ISA calculations (which assume zero temperature offset), except for the temperature and
    # density.
    #
    # The following weather parameters can be changed:
    #   - RHi: Must be in %. Set by the use of a vertical profile (set in "desired_RHis_PC") or by
    #          the
    #          use of a constant RHi (set in "RHi_background_PC"; only used when "desired_RHis_PC"
    #          is "None").
    #   - Shear: In units of 1/s. A uniform shear field is assumed. (Set by the "shear_over_s"
    #            input.)
    #   - Temperature offset: In units of K. A linear shift in the ISA temperature, can be either
    #         positive or negative.
    #   - Vertical velocity of the atmospheric air: Must be in Pa per s set in "w_Pa_per_s".
    #
    # Returns a xr.Dataset object.

    desired_altitudes_km = desired_altitudes_m / 1e3  # m to km

    # Time vector, just in case it is needed
    time_ip = np.concatenate((np.linspace(15, 23, 9), np.linspace(0, 14, 15)))  # hours

    # Initialise the dataset with bad weather data
    met_data_temp = xr.Dataset(
        data_vars=dict(
            shear=(["altitude"], np.zeros_like(desired_altitudes_m), {"units": "m/s/m"}),
            stretch=(["time"], np.ones_like(time_ip), {"units": "-"}),
            pressure=(["altitude"], np.zeros_like(desired_altitudes_m), {"units": "hPa"}),
            temperature=(["altitude"], np.zeros_like(desired_altitudes_m), {"units": "K"}),
            w=(["altitude"], np.zeros_like(desired_altitudes_m), {"units": "Pa s-1"}),
            relative_humidity_ice=(
                ["altitude"],
                np.zeros_like(desired_altitudes_m),
                {"units": "pct"},
            ),
        ),
        coords=dict(
            altitude=(["altitude"], desired_altitudes_km, {"units": "km"}),
            time=(["time"], time_ip, {"units": "h"}),
            reference_time=(np.datetime64("2019-01-01T00:00")),
        ),
    )

    if desired_RHis_PC is None:
        met_data_ISA = met_from_ISA(
            met=met_data_temp,
            RHi_background_PC=RHi_background_PC,
            T_offset_K=T_offset_K,
            shear_over_s=shear_over_s,
            w_Pa_per_s=w_Pa_per_s,
        )
    else:
        RHi_background_PC = np.min(desired_RHis_PC)
        met_data_ISA = met_from_ISA(
            met=met_data_temp,
            RHi_background_PC=RHi_background_PC,
            T_offset_K=T_offset_K,
            shear_over_s=shear_over_s,
            w_Pa_per_s=w_Pa_per_s,
        )

        # Enforce the desired RHis to each altitude
        RHi = met_data_ISA.relative_humidity_ice

        for i in range(len(desired_altitudes_km)):
            current_alt = desired_altitudes_km[i]
            mask = RHi.altitude != current_alt
            RHi = RHi.where(mask, desired_RHis_PC[i])

        met_data_ISA["relative_humidity_ice"] = RHi

    return met_data_ISA


def truncate_met(met, truncation_alt_km=20.0):
    # Truncates the APCEMM met data to remain beneath the tropopause.

    altitudes = met.altitude  # km
    truncation_idx = len(altitudes)

    # Find the index at which the altitude surpasses the truncation threshold
    for i in range(len(altitudes)):
        if altitudes[i] > truncation_alt_km:
            truncation_idx = i
            break

    return met.isel(altitude=slice(0, truncation_idx))


def met_from_ISA(met, RHi_background_PC=0.0, shear_over_s=0.002, T_offset_K=0, w_Pa_per_s=0):
    # Clears the existing temperature values of an input xr.Dataset ("met"), and
    # enforces ISA conditions at each point.
    #
    # The met created by this function is time-invariant and is only defined for 24 hours.
    #
    # The temperature offset ("T_offset_K") does not affect the ISA calculations (which
    # assume zero temperature offset), except for the temperature itself.
    #
    # The following weather parameters can be changed:
    #   - Background RHi: Must be in %. Set by the use of "RHi_background_PC".
    #   - Shear: In units of 1/s. A uniform shear field is assumed. (Set by the "shear_over_s"
    #            input.)
    #   - Temperature offset: In units of K. A shift in the ISA temperature, can be either
    #         positive or negative. Does not affect the other ISA parameters.
    #
    # Returns a xr.Dataset object.

    altitudes = met.altitude

    T = met["temperature"]
    p = met["pressure"]

    # Enforce the ISA conditions at each altitude
    for i in range(len(altitudes)):
        h_current = altitudes[i] * 1000  # km to m
        T_current = compute_T_ISA(h_current, T_offset_K)
        p_current = compute_p_ISA(h_current) / 100.0  # Pa to hPa

        mask = met.altitude != altitudes[i]
        T = T.where(mask, T_current)
        p = p.where(mask, p_current)

    met["temperature"] = T
    met["pressure"] = p

    # Set the RHi of the met data to the RHi default value
    RHi_met = met["relative_humidity_ice"]
    RHi_met = RHi_met.where(False, RHi_background_PC)
    met["relative_humidity_ice"] = RHi_met

    # Set the shear of the met data to the default value
    shear_met = met["shear"]
    shear_met = shear_met.where(False, shear_over_s)
    met["shear"] = shear_met

    # Set the vertical velocity of the met data to the RHi default value
    w_met = met["w"]
    w_met = w_met.where(False, w_Pa_per_s)
    met["w"] = w_met

    return met


def create_and_save_met(
    alts_m, RHis_PC, T_offset_K=0, shear_over_s=0.002, w_Pa_per_s=0, prefix="00"
):
    met = met_clean_slate(
        desired_altitudes_m=alts_m,
        desired_RHis_PC=RHis_PC,
        T_offset_K=T_offset_K,
        shear_over_s=shear_over_s,
        w_Pa_per_s=w_Pa_per_s,
    )

    Path("outputs/").mkdir(parents=True, exist_ok=True)
    met.to_netcdf("outputs/" + prefix + "-met.nc", engine="netcdf4")
    return met


def make_idealised_met(
    centre_altitude_m,
    moist_depth_m,
    RHi_background_PC=0.0,
    RHi_peak_PC=150,
    grad_RHi_PC_per_m=None,
    alt_resolution_m=100,
    upper_alt_m=11001,
    shear_over_s=2e-3,
    T_offset_K=0,
    w_Pa_per_s=0,
    filename_prefix="default",
):
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

    altitudes, RHis = trapezium_moist_layer(
        centre_altitude_m=centre_altitude_m,
        moist_depth_m=moist_depth_m,
        RHi_background_PC=RHi_background_PC,
        RHi_peak_PC=RHi_peak_PC,
        grad_RHi_PC_per_m=grad_RHi_PC_per_m,
        alt_resolution_m=alt_resolution_m,
        upper_alt_m=upper_alt_m,
    )

    met = create_and_save_met(
        alts_m=altitudes,
        RHis_PC=RHis,
        T_offset_K=T_offset_K,
        shear_over_s=shear_over_s,
        w_Pa_per_s=w_Pa_per_s,
        prefix=filename_prefix,
    )

    return met
