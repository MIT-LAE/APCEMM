import numpy as np

######
##    Moist layer shape Creation
######


def trapezium_moist_layer(
    centre_altitude_m,
    moist_depth_m,
    RHi_background_PC=0.0,
    RHi_peak_PC=150,
    grad_RHi_PC_per_m=None,
    alt_resolution_m=100,
    upper_alt_m=11001,
):
    # Creates a RHi (%) and altitude (m) vector that represents a trapezium in an alt vs RHi
    # graph to simulate a moist region.
    #
    # The altitude is sampled at intervals with spacing of alt_resolution_m, up until
    # upper_alt_m (not inclusive of this)
    #
    # The depth of the moist region indicates where the relative humidity is above background
    # levels.
    #
    # The generated trapezium is defined from the sides inwards, meaning that if it is too
    # narrow, the peak RHi will not be reached anywhere in the moist region.

    altitudes_m = np.arange(start=0, stop=upper_alt_m, step=alt_resolution_m)

    # Add the missing key points to the altutude array
    if not np.any(altitudes_m == centre_altitude_m):
        altitudes_m = np.append(altitudes_m, centre_altitude_m)

    thresh_upper_m = centre_altitude_m + moist_depth_m / 2.0
    thresh_lower_m = centre_altitude_m - moist_depth_m / 2.0

    if not np.any(altitudes_m == thresh_upper_m):
        altitudes_m = np.append(altitudes_m, thresh_upper_m)

    if not np.any(altitudes_m == thresh_lower_m):
        altitudes_m = np.append(altitudes_m, thresh_lower_m)

    # Add points adjacent to the threshold in order to reduce interpolation distortion of the trapezium
    if grad_RHi_PC_per_m is None:
        if not np.any(altitudes_m == thresh_upper_m + 1e-6):
            altitudes_m = np.append(altitudes_m, thresh_upper_m + 1e-6)

        if not np.any(altitudes_m == thresh_lower_m - 1e-6):
            altitudes_m = np.append(altitudes_m, thresh_lower_m - 1e-6)

    # Create the inner trapezium thresholds if necessary
    if grad_RHi_PC_per_m is not None:
        thresh_2_m = thresh_lower_m + (RHi_peak_PC - RHi_background_PC) / grad_RHi_PC_per_m
        thresh_3_m = thresh_upper_m - (RHi_peak_PC - RHi_background_PC) / grad_RHi_PC_per_m

        if thresh_3_m < thresh_2_m:
            thresh_2_m = centre_altitude_m
            thresh_3_m = centre_altitude_m
        else:
            if not np.any(altitudes_m == thresh_2_m):
                altitudes_m = np.append(altitudes_m, thresh_2_m)
            if not np.any(altitudes_m == thresh_3_m):
                altitudes_m = np.append(altitudes_m, thresh_3_m)

    # Sort the altitude array and create the RHi array
    altitudes_m = np.sort(altitudes_m)
    RHis_PC = np.zeros(altitudes_m.shape)

    # Enforce the moisture step
    mask = (altitudes_m >= thresh_lower_m) & (altitudes_m <= thresh_upper_m)
    RHis_PC = np.where(mask, RHi_peak_PC, RHi_background_PC)

    # Add the moisture gradient from the sides inwards.
    if grad_RHi_PC_per_m is not None:
        for i in range(len(altitudes_m)):
            current_alt_m = altitudes_m[i]

            if (current_alt_m < thresh_lower_m) or (current_alt_m > thresh_upper_m):
                continue

            if current_alt_m <= thresh_2_m:
                RHis_PC[i] = (
                    RHi_background_PC + (current_alt_m - thresh_lower_m) * grad_RHi_PC_per_m
                )

            if current_alt_m >= thresh_3_m:
                RHis_PC[i] = (
                    RHi_background_PC + (thresh_upper_m - current_alt_m) * grad_RHi_PC_per_m
                )

    return altitudes_m, RHis_PC
