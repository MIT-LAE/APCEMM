from src.met import make_idealised_met
import matplotlib.pyplot as plt

######################### MAIN FUNCTION (EXAMPLE) #########################

if __name__ == "__main__":
    met = make_idealised_met(
        centre_altitude_m = 10000, # m
        moist_depth_m = 1000, # m
        RHi_background_PC = 40, # m
        RHi_peak_PC = 140, # %
        grad_RHi_PC_per_m = None, # % RHi / m; None indicates a rectangular moist region
        alt_resolution_m = 100, # m
        upper_alt_m = 11001, # m
        shear_over_s = 2e-3, # 1/s
        T_offset_K = 0, # K
        filename_prefix = "example"
    )

    RHi_PC = met.relative_humidity_ice.to_numpy()
    alt_m = met.altitude.to_numpy() * 1e3 # km to m

    plt.plot(RHi_PC, alt_m, color = 'b', lw = 1.5)
    plt.ylabel("Altitude, m")
    plt.xlabel('RHi, %')
    plt.yticks([9000, 9500, 10000, 10500, 11000])
    plt.ylim(9000,11000)
    plt.xlim(0,150)
    plt.show()
