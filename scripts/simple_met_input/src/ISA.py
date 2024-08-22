import numpy as np

######
##    ISA Equations and constants from https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Fundamentals_of_Aerospace_Engineering_(Arnedo)/02%3A_Generalities/2.03%3A_Standard_atmosphere
######

# TODO: Include humidity in these calculations...

class ISA:
    # ISA Class needed to store the constants without making them global

    # ISA Temperature Lapse Rate Magnitude (Troposphere)
    alpha = 6.5e-3 # K / m

    # Perfect Gas Constant for Air
    R = 287.053 # J / kg K

    # Acceleration Due to Gravity
    g = 9.80665 # m / s^2

    # Sea level quantities
    T_0 = 288.15 # K
    p_0 = 101325 # Pa
    rho_0 = 1.225 # kg / m^3

    # Sea level temperature offfset
    T_offset = 0 # K

    # Quantities at the boundary between the troposphere and the stratosphere (11000m)
    h_11 = 11000 # m
    p_11 = np.nan # Calculated during the object init
    T_11 = np.nan # Calculated during the object init

    def __init__(self, T_offset_K = 0):
        # Evaluate the quantities at the boundary between the troposphere and the stratosphere (11000m) upon initialisation
        self.T_offset = T_offset_K
        
        self.p_11 = self.p_0 * np.power(1 - self.alpha / self.T_0 * self.h_11,
                                        self.g / self.R / self.alpha)
        T_11 = self.T_0 - self.alpha * self.h_11 + self.T_offset
        
    def compute_p_from_h(self, h):
        # Returns the pressure (in Pa) at a given altitude h (in metres)
        # according to the International Standard Atmosphere from the following link:
        #  
        # https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Fundamentals_of_Aerospace_Engineering_(Arnedo)/02%3A_Generalities/2.03%3A_Standard_atmosphere        
        
        p_troposphere = self.p_0 * np.power(1 - self.alpha / self.T_0 * h,
                                            self.g / self.R / self.alpha)
        p_stratosphere = self.p_11 * np.exp(-self.g / self.R / self.T_11 * (h - self.h_11))

        p_out = np.where(h <= self.h_11, p_troposphere, p_stratosphere)

        # Enforce altitudes beneath SL to be equal to SL values.
        return np.where(p_out > self.p_0, self.p_0, p_out)

    def compute_T_from_h(self, h):
        # Returns the temperature (in degrees K) at a given altitude h (in metres)
        # according to the International Standard Atmosphere from the following link:
        #    
        # https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Fundamentals_of_Aerospace_Engineering_(Arnedo)/02%3A_Generalities/2.03%3A_Standard_atmosphere
        
        T_troposphere = self.T_0 - self.alpha * h
        T_stratosphere = self.T_11

        T_out = np.where(h <= self.h_11, T_troposphere, T_stratosphere)

        # Enforce altitudes beneath SL to be equal to SL values. Also apply the T offset.
        return self.T_offset + np.where(T_out > self.T_0, self.T_0, T_out)

    def compute_rho_from_h(self, h):
        # Returns the air density (in kg / m^3) at a given altitude h (in metres)
        # according to the International Standard Atmosphere from the following link:
        #  
        # https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Fundamentals_of_Aerospace_Engineering_(Arnedo)/02%3A_Generalities/2.03%3A_Standard_atmosphere
        
        p =  self.compute_p_from_h(h)
        T =  self.compute_T_from_h(h)

        return p / (self.R * T)

    def compute_h_from_p(self, p):
        # Returns the pressure (in kPa) at a given altitude h (in metres)
        # according to the International Standard Atmosphere from the following link:
        #    
        # https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Fundamentals_of_Aerospace_Engineering_(Arnedo)/02%3A_Generalities/2.03%3A_Standard_atmosphere
        
        h_troposphere = self.T_0 / self.alpha * ( 1 - np.power(p / self.p_0,
                                                               self.R * self.alpha / self.g) )
        h_stratosphere = self.h_11 - self.R * self.T_11 / self.g * np.log(p/self.p_11)
        
        h_out = np.where(p >= self.p_11, h_troposphere, h_stratosphere)

        # Enforce altitudes beneath SL to be equal to the SL.
        return np.where(h_out < 0, 0, h_out)

def compute_p_ISA(h):
    # Returns the pressure (in Pa) at a given altitude h (in metres)
    # according to the International Standard Atmosphere from the following link:
    #
    # https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Fundamentals_of_Aerospace_Engineering_(Arnedo)/02%3A_Generalities/2.03%3A_Standard_atmosphere        
    
    ISA_computer = ISA()
    return ISA_computer.compute_p_from_h(h)

def compute_T_ISA(h, T_offset_K = 0):
    # Returns the temperature (in degrees K) at a given altitude h (in metres)
    # according to the International Standard Atmosphere from the following link:
    #  
    # https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Fundamentals_of_Aerospace_Engineering_(Arnedo)/02%3A_Generalities/2.03%3A_Standard_atmosphere
    
    ISA_computer = ISA(T_offset_K)
    return ISA_computer.compute_T_from_h(h)

def compute_rho_ISA(h, T_offset_K = 0):
    # Returns the air density (in kg / m^3) at a given altitude h (in metres)
    # according to the International Standard Atmosphere from the following link:
    #
    # NOTE: NO DISTINCTION IS MADE BETWEEN DRY AND HUMID AIR DENSITY HERE!
    #
    # https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Fundamentals_of_Aerospace_Engineering_(Arnedo)/02%3A_Generalities/2.03%3A_Standard_atmosphere

    ISA_computer = ISA(T_offset_K)
    return ISA_computer.compute_rho_from_h(h)

def compute_h_from_p_ISA(p):
    # Returns the pressure (in kPa) at a given altitude h (in metres)
    # according to the International Standard Atmosphere from the following link:
    #    
    # https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Fundamentals_of_Aerospace_Engineering_(Arnedo)/02%3A_Generalities/2.03%3A_Standard_atmosphere
    
    ISA_computer = ISA()
    return ISA_computer.compute_h_from_p(p)
