import os
import chaospy
import os.path
import pickle
import time
import numpy as np
import netCDF4 as nc
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable



"""
**********************************
WRITING APCEMM VARIABLES FUNCTIONS
**********************************
"""
def set_temp_K(lines : list, T : float) -> list:
    """Sets thetemperature (T / K) in the lines from input.yaml (lines)"""
    newlines = lines.copy()

    line =  "    Temperature [K] (double): " + str(T) + "\n"
    newlines[38] = line

    return newlines

def set_RH_percent(lines : list, RH : float) -> list:
    """Sets the relative humidity (RH / %) in the lines from input.yaml (lines)"""
    newlines = lines.copy()

    line =  "    R.Hum. wrt water [%] (double): " + str(RH) + "\n"
    newlines[39] = line

    return newlines

def set_p_hPa(lines : list, p : float) -> list:
    """Sets the pressure (p / hPa) in the lines from input.yaml (lines)"""
    newlines = lines.copy()

    line =  "    Pressure [hPa] (double): " + str(p) + "\n"
    newlines[40] = line

    return newlines

def set_coords_deg(lines : list, lon : float, lat : float) -> list:
    """Sets the longitude (lon / deg) and latitude (lat / deg) 
    in the lines from input.yaml (lines)"""
    newlines = lines.copy()

    line =  "    LON [deg] (double): " + str(lon) + "\n"
    newlines[46] = line

    line =  "    LAT [deg] (double): " + str(lat) + "\n"
    newlines[47] = line

    return newlines

def set_day(lines : list, day : int) -> list:
    """Sets the day (1-365) in the lines from input.yaml (lines)"""
    newlines = lines.copy()

    line =  "    Emission day [1-365] (int): " + str(day) + "\n"
    newlines[48] = line

    return newlines

def set_time_hrs_UTC(lines : list, hr : float) -> list:
    """Sets the time (24hr format UTC) in the lines from input.yaml (lines)"""
    newlines = lines.copy()

    line =  "    Emission time [hr] (double) : " + str(hr) + "\n"
    newlines[49] = line

    return newlines

def set_EI_soot_gPerkg(lines : list, EI_soot : float) -> list:
    """Sets the soot Emissions Index (EI_soot / g of soot per kg of fuel) 
    in the lines from input.yaml (lines)"""
    newlines = lines.copy()

    line =  "    Soot [g/kg_fuel] (double): " + str(EI_soot) + "\n"
    newlines[63] = line

    return newlines

def set_fuel_flow_kgPers(lines : list, fuel_flow : float) -> list:
    """Sets the fuel flow (fuel_flow / kg per s) in the lines from input.yaml (lines)"""
    newlines = lines.copy()

    line =  "  Total fuel flow [kg/s] (double) : " + str(fuel_flow) + "\n"
    newlines[65] = line

    return newlines

def set_aircraft_mass_kg(lines : list, aircraft_mass : float) -> list:
    """Sets the aircraft mass (aircraft_mass / kg) in the lines from input.yaml (lines)"""
    newlines = lines.copy()

    line =  "  Aircraft mass [kg] (double): " + str(aircraft_mass) + "\n"
    newlines[66] = line

    return newlines

def set_flight_speed_mPers(lines : list, flight_speed : float) -> list:
    """Sets the flight speed (flight_speed / m/s) in the lines from input.yaml (lines)"""
    newlines = lines.copy()

    line =  "  Flight speed [m/s] (double): " + str(flight_speed) + "\n"
    newlines[67] = line

    return newlines

def set_core_exit_temp_K(lines : list, T_core_exit : float) -> list:
    """Sets the core exit temperature (T_core_exit / K) in the lines from input.yaml (lines)"""
    newlines = lines.copy()

    line =  "  Core exit temp. [K] (double): " + str(T_core_exit) + "\n"
    newlines[70] = line

    return newlines

def default_APCEMM_vars():
    location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

    # Read the input file
    ip_file = open(os.path.join(location,'original.yaml'), 'r')
    op_lines = ip_file.readlines()
    ip_file.close()

    op_file = open(os.path.join(location,'input.yaml'), 'w')
    op_file.writelines(op_lines)
    op_file.close()

def write_APCEMM_vars(temp_K = 217, RH_percent = 63.94, p_hPa = 250.0, lat_deg = 20.2, 
               lon_deg = 20.2, day = 20, time_hrs_UTC = 20.0, EI_soot_gPerkg = 0.008,
               fuel_flow_kgPers = 2.8, aircraft_mass_kg = 3.10e+05, 
               flight_speed_mPers = 250.0, core_exit_temp_K = 560.0):
    
    location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

    # Read the input file
    ip_file = open(os.path.join(location,'input.yaml'), 'r')
    op_lines = ip_file.readlines()
    ip_file.close()

    # Write the variables to input.yaml
    op_lines = set_temp_K(op_lines, temp_K)
    op_lines = set_RH_percent(op_lines, RH_percent)
    op_lines = set_p_hPa(op_lines, p_hPa)
    op_lines = set_coords_deg(op_lines, lon_deg, lat_deg)
    op_lines = set_day(op_lines, day)
    op_lines = set_time_hrs_UTC(op_lines, time_hrs_UTC)
    op_lines = set_EI_soot_gPerkg(op_lines, EI_soot_gPerkg)
    op_lines = set_fuel_flow_kgPers(op_lines, fuel_flow_kgPers)
    op_lines = set_aircraft_mass_kg(op_lines, aircraft_mass_kg)
    op_lines = set_flight_speed_mPers(op_lines, flight_speed_mPers)
    op_lines = set_core_exit_temp_K(op_lines, core_exit_temp_K)

    op_file = open(os.path.join(location,'input.yaml'), 'w')
    op_file.writelines(op_lines)
    op_file.close()

def write_APCEMM_nipc_vars(nipc_vars):
    
    location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

    # Read the input file
    ip_file = open(os.path.join(location,'input.yaml'), 'r')
    op_lines = ip_file.readlines()
    ip_file.close()

    for var in nipc_vars:
        if var.name == "temp_K":
            op_lines = set_temp_K(op_lines, var.data)
            continue
        if var.name == "RH_percent":
            op_lines = set_RH_percent(op_lines, var.data)
            continue
        if var.name == "EI_soot_gPerkg":
            op_lines = set_EI_soot_gPerkg(op_lines, var.data)
            continue
        if var.name == "fuel_flow_kgPers":
            op_lines = set_fuel_flow_kgPers(op_lines, var.data)
            continue
        if var.name == "aircraft_mass_kg":
            op_lines = set_aircraft_mass_kg(op_lines, var.data)
            continue
        if var.name == "flight_speed_mPers":
            op_lines = set_flight_speed_mPers(op_lines, var.data)
            continue
        if var.name == "core_exit_temp_K":
            op_lines = set_core_exit_temp_K(op_lines, var.data)
            continue
        if var.name == "time_hrs_UTC":
            op_lines = set_time_hrs_UTC(op_lines, var.data)
            continue
        if var.name == "p_hPa":
            op_lines = set_p_hPa(op_lines, var.data)
            continue

    op_file = open(os.path.join(location,'input.yaml'), 'w')
    op_file.writelines(op_lines)
    op_file.close()



"""
**********************************
READING APCEMM OUTPUTS FUNCTIONS
**********************************
"""
# This code is taken from the APCEMM example files
def read_nc_file(filename):
    ip_filepath_base = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    ip_filepath_base = ip_filepath_base + "/APCEMM_out"

    ip_filepath = os.path.join(ip_filepath_base,filename)

    ds = nc.Dataset(ip_filepath)

    print(ds.variables)

    print(ds.variables['intOD'][:])

    return ds
    
def read_apcemm_data(directory):
    t_mins = []
    # optical_depth_vert_int = []
    ice_particles = []
    # ds_t = []
    # optical_depth_vert = []
    # optical_depth_horiz = []

    for file in sorted(os.listdir(directory)):
        if(file.startswith('ts_aerosol') and file.endswith('.nc')):
            file_path = os.path.join(directory,file)
            ds = xr.open_dataset(file_path, engine = "netcdf4", decode_times = False)
            tokens = file_path.split('.')
            mins = int(tokens[-2][-2:])
            hrs = int(tokens[-2][-4:-2])
            t_mins.append(hrs*60 + mins)
            #print(ds.variables['Number Ice Particles'])
            ice_particles.append(ds.variables['Number Ice Particles'][:].values[0])
            # optical_depth_vert_int.append(ds.variables['intOD'][:].values[0])
            # optical_depth_horiz.append(ds["Horizontal optical depth"])
            # optical_depth_vert.append(ds["Vertical optical depth"])
            # ds_t.append(ds)

    return t_mins, ice_particles
    # return apce_data_struct(t_mins, ds_t, optical_depth_vert, optical_depth_horiz)

def removeLow(arr, cutoff = 1e-3):
    func = lambda x: (x > cutoff) * x
    vfunc = np.vectorize(func)
    return vfunc(arr)



"""
**********************************
NIPC FUNCTIONS
**********************************
"""
class NIPC_var:
    def __init__(self, name, data):
        self.name = name
        self.data = data

# The model NIPC is being applied on
def eval_model(sample):
    if sample.ndim < 1:
        sample = np.array([sample])

    var_rh = NIPC_var("RH_percent", sample[0])
    nipc_vars = [var_rh]

    # Default the variables
    default_APCEMM_vars()

    # Write the specific variables one by one
    write_APCEMM_nipc_vars(nipc_vars)

    # Run APCEMM
    os.system('./../../Code.v05-00/APCEMM input.yaml')

    # Read the output
    directory = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    directory = directory + "/APCEMM_out"
    t_mins, output = read_apcemm_data(directory)

    # Return the integrated optical depth
    return output

# See Equation A.5 from FYR
def transform(samples, distribution_input, distribution_germ):
    return distribution_input.inv(distribution_germ.fwd(samples))



"""
**********************************
MAIN FUNCTION
**********************************
"""
if __name__ == "__main__" :
    timing = False

    # Chaospy code from https://chaospy.readthedocs.io/en/master/user_guide/advanced_topics/generalized_polynomial_chaos.html
    # Using point collocation
    directory = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    directory = directory + "/APCEMM_out"

    var_rh = NIPC_var("RH_percent", 150)
    default_APCEMM_vars() # Default the variables
    write_APCEMM_nipc_vars([var_rh]) # Write the specific variables one by one
    os.system('./../../Code.v05-00/APCEMM input.yaml') # Run APCEMM
    times, optical_depth_int = read_apcemm_data(directory)

    # Input variable distribution
    dist_input = chaospy.Uniform(95,105)
    # Multivar: dist_input = chaospy.MvNormal(mu=[10, 1], sigma=[[1.0, 0.09], [0.09, 0.1]])
    
    # Germ distribution
    dist_germ = chaospy.Uniform(-1, 1)
    # Multivar: dist_germ = chaospy.J(chaospy.Normal(0, 1), chaospy.Normal(0, 1))

    # Sample the germ
    samples_r = dist_germ.sample(1, rule="sobol", seed=0)

    # Match the germ samples to the input variable space
    samples_q = transform(samples_r, dist_input, dist_germ)

    # Create orthogonal polynomial expansion according to the distribution (7th order)
    expansion = chaospy.generate_expansion(7, dist_germ)

    if timing:
        start = time.time()
    
    # Evaluate the deterministic samples of the output variable
    evaluations = np.array([eval_model(sample) for sample in samples_q.T])

    if timing:
        end = time.time()
        print("\n\n" + str(end - start))

    if np.any(evaluations < 0):
        print("Negative optical depth value found!")

    # Save the evaluations
    DF = pd.DataFrame(evaluations)
    DF.to_csv("APCEMM-NIPC-evaluations.csv")

    # Save the time vector
    DF = pd.DataFrame(times)
    DF.to_csv("APCEMM-NIPC-times.csv")

    # Pickle the PCE of the germ
    file = open("APCEMM-NIPC-expansion.chp", 'wb')
    pickle.dump(expansion, file)
    file.close()

    # Pickle the germ samples
    file = open("APCEMM-NIPC-samples_r.chp", 'wb')
    pickle.dump(samples_r, file)
    file.close()
    