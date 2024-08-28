import numpy as np

######
##    ISA equations and constants from https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Fundamentals_of_Aerospace_Engineering_(Arnedo)/02%3A_Generalities/2.03%3A_Standard_atmosphere
######

# TODO: Include humidity in these calculations...

class ISA:
    # ISA temperature lapse rate magnitude at the troposphere
    alpha = 6.5e-3 # K / m

    # Perfect gas constant for air
    R = 287.053 # J / kg K

    # Acceleration due to gravity
    g = 9.80665 # m / s^2

    # Sea level quantities
    T_0 = 288.15 # K
    p_0 = 101325 # Pa
    rho_0 = 1.225 # kg / m^3

    # Sea level temperature offfset
    T_offset = 0 # K

    # Quantities at the tropopause (11000m AMSL)
    h_11 = 11000 # m AMSL
    p_11 = np.nan # Calculated during the object initialisation
    T_11 = np.nan # Calculated during the object initialisation

    def __init__(self, T_offset_K = 0):
        # Evaluate the quantities at the tropopause (11000m AMSL) upon initialisation
        # according to the International Standard Atmosphere equations from the following link:
        #  
        # https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Fundamentals_of_Aerospace_Engineering_(Arnedo)/02%3A_Generalities/2.03%3A_Standard_atmosphere    

        self.T_offset = T_offset_K
        
        self.p_11 = self.p_0 * np.power(1 - self.alpha / self.T_0 * self.h_11,
                                        self.g / self.R / self.alpha)
        T_11 = self.T_0 - self.alpha * self.h_11 + self.T_offset
        
    def compute_p_from_h(self, h_m):
        # Returns the pressure (in Pa) at a given altitude h (in metres AMSL)
        # according to the International Standard Atmosphere equations from the following link:
        #  
        # https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Fundamentals_of_Aerospace_Engineering_(Arnedo)/02%3A_Generalities/2.03%3A_Standard_atmosphere        
        
        p_troposphere_Pa = self.p_0 * np.power(1 - self.alpha / self.T_0 * h_m,
                                            self.g / self.R / self.alpha)
        p_stratosphere_Pa = self.p_11 * np.exp(-self.g / self.R / self.T_11 * (h_m - self.h_11))

        p_out_Pa = np.where(h_m <= self.h_11, p_troposphere_Pa, p_stratosphere_Pa)

        # Enforce altitudes beneath SL to be equal to SL values.
        return np.where(p_out_Pa > self.p_0, self.p_0, p_out_Pa)

    def compute_T_from_h(self, h_m):
        # Returns the temperature (in degrees K) at a given altitude h (in metres AMSL)
        # according to the International Standard Atmosphere equations from the following link:
        #    
        # https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Fundamentals_of_Aerospace_Engineering_(Arnedo)/02%3A_Generalities/2.03%3A_Standard_atmosphere
        
        T_troposphere_K = self.T_0 - self.alpha * h_m
        T_stratosphere_K = self.T_11

        T_out_K = np.where(h_m <= self.h_11, T_troposphere_K, T_stratosphere_K)

        # Enforce altitudes beneath SL to be equal to SL values. Also apply the T offset.
        return self.T_offset + np.where(T_out_K > self.T_0, self.T_0, T_out_K)

    def compute_rho_from_h(self, h_m):
        # Returns the air density (in kg / m^3) at a given altitude h (in metres AMSL)
        # according to the ideal gas law.
        #
        # NOTE: NO DISTINCTION IS MADE BETWEEN DRY AND HUMID AIR DENSITY HERE!
        
        p_Pa =  self.compute_p_from_h(h_m)
        T_K =  self.compute_T_from_h(h_m)

        return p_Pa / (self.R * T_K)

    def compute_h_from_p(self, p_Pa):
        # Returns the altitude (in metres AMSL) at a given altitude ambient pressure (in Pa)
        # according to the International Standard Atmosphere equations from the following link:
        #    
        # https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Fundamentals_of_Aerospace_Engineering_(Arnedo)/02%3A_Generalities/2.03%3A_Standard_atmosphere
        
        h_troposphere_m = self.T_0 / self.alpha * ( 1 - np.power(p_Pa / self.p_0,
                                                               self.R * self.alpha / self.g) )
        h_stratosphere_m = self.h_11 - self.R * self.T_11 / self.g * np.log(p_Pa / self.p_11)
        
        h_out_m = np.where(p_Pa >= self.p_11, h_troposphere_m, h_stratosphere_m)

        # Enforce altitudes beneath SL to be equal to the SL.
        return np.where(h_out_m < 0, 0, h_out_m)

def compute_p_ISA(h_m):
    # Returns the pressure (in Pa) at a given altitude h (in metres AMSL)
    # according to the International Standard Atmosphere equations from the following link:
    #
    # https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Fundamentals_of_Aerospace_Engineering_(Arnedo)/02%3A_Generalities/2.03%3A_Standard_atmosphere        
    
    ISA_computer = ISA()
    return ISA_computer.compute_p_from_h(h_m)

def compute_T_ISA(h_m, T_offset_K = 0):
    # Returns the temperature (in degrees K) at a given altitude h (in metres AMSL)
    # according to the International Standard Atmosphere equations from the following link:
    #  
    # https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Fundamentals_of_Aerospace_Engineering_(Arnedo)/02%3A_Generalities/2.03%3A_Standard_atmosphere
    
    ISA_computer = ISA(T_offset_K)
    return ISA_computer.compute_T_from_h(h_m)

def compute_rho_ISA(h_m, T_offset_K = 0):
    # Returns the air density (in kg / m^3) at a given altitude h (in metres AMSL)
    # according to the ideal gas law.
    #
    # NOTE: NO DISTINCTION IS MADE BETWEEN DRY AND HUMID AIR DENSITY HERE!

    ISA_computer = ISA(T_offset_K)
    return ISA_computer.compute_rho_from_h(h_m)

def compute_h_from_p_ISA(p_Pa):
    # Returns the altitude (in metres AMSL) at a given altitude ambient pressure (in Pa)
    # according to the International Standard Atmosphere equations from the following link:
    #    
    # https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Fundamentals_of_Aerospace_Engineering_(Arnedo)/02%3A_Generalities/2.03%3A_Standard_atmosphere
    
    ISA_computer = ISA()
    return ISA_computer.compute_h_from_p(p_Pa)
