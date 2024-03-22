import os
import chaospy
import shutil
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

def write_APCEMM_NIPC_vars(NIPC_vars):
    
    location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

    # Read the input file
    ip_file = open(os.path.join(location,'input.yaml'), 'r')
    op_lines = ip_file.readlines()
    ip_file.close()

    for var in NIPC_vars:
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
    
def read_APCEMM_data(directory, output_id):
    """ 
    Supported output_id values:
        - "Horizontal optical depth"
        - "Vertical optical depth"
        - "Number Ice Particles" (#/m)
        - "Ice Mass" (Ice mass of contrail section per unit length (kg/m))
        - "intOD" (Vertical optical depth integrated over the grid)
    """
    t_mins = []
    output = []

    for file in sorted(os.listdir(directory)):
        if(file.startswith('ts_aerosol') and file.endswith('.nc')):
            file_path = os.path.join(directory,file)
            ds = xr.open_dataset(file_path, engine = "netcdf4", decode_times = False)
            tokens = file_path.split('.')
            mins = int(tokens[-2][-2:])
            hrs = int(tokens[-2][-4:-2])
            t_mins.append(hrs*60 + mins)

            if (output_id == "Horizontal optical depth") | (output_id == "Vertical optical depth"):
                output.append(ds[output_id])
            else:
                output.append(ds.variables[output_id][:].values[0])

    if len(t_mins) == 0:
        t_mins.append(0)

    while len(t_mins) < 37:
        t_mins.append(t_mins[-1] + 10)

    while len(output) < 37:
        output.append(0)

    return t_mins, output

def reset_APCEMM_outputs(directory):
    for file in sorted(os.listdir(directory)):
        if(file.startswith('ts_aerosol') and file.endswith('.nc')):
            file_path = os.path.join(directory,file)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print('Failed to delete %s. Reason: %s' % (file_path, e))

def removeLow(arr, cutoff = 1e-3):
    func = lambda x: (x > cutoff) * x
    vfunc = np.vectorize(func)
    return vfunc(arr)

"""
**********************************
THERMODYNAMICS FUNCTIONS
**********************************
"""
def compute_p_sat_liq(T_K):
    # Calculates the liquid saturation pressure "p_sat_liq" (in units of Pa).
    # This equation can be found in Schumann 2012.
    # % 
    # % The inputs are as follows
    # %	- T_K is the absolute temperature in Kelvin
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

"""
**********************************
NIPC FUNCTIONS
**********************************
"""
class NIPC_var:
    """
    Supported NIPC_var.names:
       - "temp_K": Air Temperature in K
       - "RH_percent": Relative humidity wrt water in %
       - "EI_soot_gPerkg": Soot emissions index in grams per kg
       - "fuel_flow_kgPers": Fuel flow in kg per seconds
       - "aircraft_mass_kg": Aircraft mass in kg
       - "flight_speed_mPers": TAS at cruise in m/s
       - "core_exit_temp_K": Engine core exit temperature in K
       - "time_hrs_UTC": Time in hours UTC
       - "p_hPa": Pressure at the flight altitude in hPa
    """
        
    def __init__(self, name, data):
        self.name = name
        self.data = data

def eval_APCEMM(NIPC_vars, directory, output_id = "Number Ice Particles"):
    # Supported NIPC_var.names:
    #   - "temp_K"
    #   - "RH_percent"
    #   - "EI_soot_gPerkg"
    #   - "fuel_flow_kgPers"
    #   - "aircraft_mass_kg"
    #   - "flight_speed_mPers"
    #   - "core_exit_temp_K"
    #   - "time_hrs_UTC"
    #   - "p_hPa"
    #
    #
    # Supported output_id values:
    #     - "Horizontal optical depth"
    #     - "Vertical optical depth"
    #     - "Number Ice Particles" (#/m)
    #     - "Ice Mass" (Ice mass of contrail section per unit length (kg/m))
    #     - "intOD" (Vertical optical depth integrated over the grid)

    # Default the variables
    default_APCEMM_vars()

    # Write the specific variables one by one
    write_APCEMM_NIPC_vars(NIPC_vars)

    # Eliminate the output files
    reset_APCEMM_outputs(directory)

    # Run APCEMM
    os.system('./../../Code.v05-00/APCEMM input.yaml')

    # Read the output
    t_mins, output = read_APCEMM_data(directory, output_id=output_id)

    # Return the output
    return t_mins, output

# The model NIPC is being applied on
def eval_model_NIPC(sample, directory, output_id = "Number Ice Particles"):
    if sample.ndim < 1:
        sample = np.array([sample])

    var_RH = NIPC_var("RH_percent", sample[0])
        
    if len(sample) > 1:
        var_T = NIPC_var("temp_K", sample[1])
        NIPC_vars = [var_RH, var_T]
    else:
        NIPC_vars = [var_RH]

    # Read the output
    t_mins, output = eval_APCEMM(NIPC_vars=NIPC_vars, 
                                 directory=directory, output_id=output_id)

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
    # Chaospy code from https://chaospy.readthedocs.io/en/master/user_guide/advanced_topics/generalized_polynomial_chaos.html
    # Using point collocation
    timing = False
    output_id = "Number Ice Particles"
    runs = 100

    directory = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    directory = directory + "/APCEMM_out"

    RHi_default = 150
    T_default = 217
    var_RH = NIPC_var("RH_percent", convert_RHi_to_RH(T_default, RHi_default))
    var_T = NIPC_var("temp_K", T_default)
    times, evaluations = eval_APCEMM([var_RH, var_T], directory = directory, output_id=output_id)

    # Save the evaluations
    DF = pd.DataFrame(evaluations)
    DF.to_csv(directory + "APCEMM-debug-evaluations.csv")

    # Save the time vector
    DF = pd.DataFrame(times)
    DF.to_csv(directory + "APCEMM-debug-times.csv")

