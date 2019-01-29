/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Input Program File                                               */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 11/19/2018                                */
/* File                 : Input.cpp                                 */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Core/Input.hpp"

Input::Input( const RealDouble temperature_K, \
              const RealDouble pressure_Pa,   \
              const RealDouble relHumidity_w, \
              const RealDouble long_deg,      \
              const RealDouble lat_deg,       \
              const UInt dayGMT,              \
              const RealDouble emissionTime,  \
              const Vector_1D emissionInput,  \
              const Vector_1D backgMixRatio,  \
              const std::string fileName,     \
              const std::string fileName_ADJ ):
    temperature_K_ ( temperature_K ),
    pressure_Pa_   ( pressure_Pa) ,
    relHumidity_w_ ( relHumidity_w ),
    longitude_deg_ ( long_deg ),
    latitude_deg_  ( lat_deg ),
    emissionDOY_   ( dayGMT ),
    emissionTime_  ( emissionTime ),
    EI_NOx_        ( emissionInput[0] ),
    EI_CO_         ( emissionInput[1] ),
    EI_HC_         ( emissionInput[2] ),
    EI_SO2_        ( emissionInput[3] ),
    EI_Soot_       ( emissionInput[4] ),
    sootRad_       ( emissionInput[5] ),
    fuelFlow_      ( emissionInput[6] ),
    backgNOx_      ( backgMixRatio[0] ),
    backgHNO3_     ( backgMixRatio[1] ),
    backgO3_       ( backgMixRatio[2] ),
    backgCO_       ( backgMixRatio[3] ),
    backgCH4_      ( backgMixRatio[4] ),
    backgSO2_      ( backgMixRatio[5] ),
    fileName_      ( fileName ),
    fileName_ADJ_  ( fileName_ADJ )
{

    /* Constructor */

    while ( longitude_deg_ > 180 )
        longitude_deg_ -= 360;
    
    while ( longitude_deg_ < -180 )
        longitude_deg_ += 360;

    while ( latitude_deg_ > 90 )
        latitude_deg_ -= 180;
    
    while ( latitude_deg_ < -90 )
        latitude_deg_ += 180;

    while ( emissionDOY_ > 365 )
        emissionDOY_ -= 365;

    while ( emissionTime_ > 24.0 )
        emissionTime_ -= 24.0;

    if ( temperature_K_ >= 300.0 || temperature_K_ <= 160.0 ) {
        std::cout << " In Input::Input:";
        std::cout << " temperature_K takes an odd value: temperature_K = ";
        std::cout << temperature_K_ << " [K]" << std::endl;
        exit(-1);
    }
    
    if ( pressure_Pa_ >= 1.00E+05 || pressure_Pa_ <= 1.00E+03 ) {
        std::cout << " In Input::Input:";
        std::cout << " pressure_Pa takes an odd value: pressure_Pa_ = ";
        std::cout << pressure_Pa_ << " [Pa]" << std::endl;
        exit(-1);
    }

    if ( relHumidity_w_ >= 9.50E+01 || pressure_Pa_ <= 0.00E+00 ) {
        std::cout << " In Input::Input:";
        std::cout << " relHumidity_w takes an unrealisable value: relHumidity_w = ";
        std::cout << relHumidity_w_ << " [%]" << std::endl;
        exit(-1);
    }

    if ( emissionDOY_ <= 0.0 ) {
        std::cout << " In Input::Input:";
        std::cout << " emissionDOY takes an unrealisable value: emissionDOY = ";
        std::cout << emissionDOY_ << " [-]" << std::endl;
        exit(-1);
    }

    if ( EI_NOx_ < 0.0E+00 || EI_NOx_ > 5.0E+01 ) {
        std::cout << " In Input::Input:";
        std::cout << " EI_NOx takes an unrealisable value: EI_NOx = ";
        std::cout << EI_NOx_ << " [g/kg_fuel]" << std::endl;
        exit(-1);
    }
    
    if ( EI_CO_ < 0.0E+00 || EI_CO_ > 3.0E+01 ) {
        std::cout << " In Input::Input:";
        std::cout << " EI_CO takes an unrealisable value: EI_CO = ";
        std::cout << EI_CO_ << " [g/kg_fuel]" << std::endl;
        exit(-1);
    }
    
    if ( EI_HC_ < 0.0E+00 || EI_HC_ > 1.0E+01 ) {
        std::cout << " In Input::Input:";
        std::cout << " EI_HC takes an unrealisable value: EI_HC = ";
        std::cout << EI_HC_ << " [g/kg_fuel]" << std::endl;
        exit(-1);
    }
    
    if ( EI_SO2_ < 0.0E+00 || EI_SO2_ > 1.0E+02 ) {
        std::cout << " In Input::Input:";
        std::cout << " EI_SO2 takes an unrealisable value: EI_SO2 = ";
        std::cout << EI_SO2_ << " [g/kg_fuel]" << std::endl;
        exit(-1);
    }

    if ( EI_Soot_ < 0.0E+00 || EI_Soot_ > 2.0E-01 ) {
        std::cout << " In Input::Input:";
        std::cout << " EI_Soot takes an unrealisable value: EI_Soot = ";
        std::cout << EI_Soot_ << " [g/kg_fuel]" << std::endl;
        exit(-1);
    }
    
    if ( ( ( sootRad_ < 1.0E-10 ) && ( sootRad_ != 0.0E+00 ) ) || sootRad_ > 1.0E-07 ) {
        std::cout << " In Input::Input:";
        std::cout << " sootRad takes an unrealisable value: sootRad = ";
        std::cout << sootRad_ * 1.0E+09 << " [nm]" << std::endl;
        exit(-1);
    }
    
    if ( fuelFlow_ < 0.0E+00 ) {
        std::cout << " In Input::Input:";
        std::cout << " fuelFlow takes an unrealisable value: fuelFlow = ";
        std::cout << fuelFlow_ << " [kg/s]" << std::endl;
        exit(-1);
    }
    
    if ( backgNOx_ < 0.0E+00 || backgNOx_ > 1.0E+09 ) {
        std::cout << " In Input::Input:";
        std::cout << " backgNOx takes an unrealisable value: backgNOx = ";
        std::cout << backgNOx_ << " [ppb]" << std::endl;
        exit(-1);
    }
    
    if ( backgHNO3_ < 0.0E+00 || backgHNO3_ > 1.0E+09 ) {
        std::cout << " In Input::Input:";
        std::cout << " backgHNO3 takes an unrealisable value: backgHNO3 = ";
        std::cout << backgHNO3_ << " [ppb]" << std::endl;
        exit(-1);
    }
    
    if ( backgO3_ < 0.0E+00 || backgO3_ > 1.0E+09 ) {
        std::cout << " In Input::Input:";
        std::cout << " backgO3 takes an unrealisable value: backgO3 = ";
        std::cout << backgO3_ << " [ppb]" << std::endl;
        exit(-1);
    }
    
    if ( backgCO_ < 0.0E+00 || backgCO_ > 1.0E+09 ) {
        std::cout << " In Input::Input:";
        std::cout << " backgCO takes an unrealisable value: backgCO = ";
        std::cout << backgCO_ << " [ppb]" << std::endl;
        exit(-1);
    }
    
    if ( backgCH4_ < 0.0E+00 || backgCH4_ > 1.0E+09 ) {
        std::cout << " In Input::Input:";
        std::cout << " backgCH4 takes an unrealisable value: backgCH4 = ";
        std::cout << backgCH4_ << " [ppb]" << std::endl;
        exit(-1);
    }
    
    if ( backgSO2_ < 0.0E+00 || backgSO2_ > 1.0E+09 ) {
        std::cout << " In Input::Input:";
        std::cout << " backgSO2 takes an unrealisable value: backgCH4 = ";
        std::cout << backgSO2_ << " [ppb]" << std::endl;
        exit(-1);
    }

    if ( emissionDOY_ <= 31 ) {
        emissionMonth_ = 1;
        emissionDay_   = emissionDOY_;
    } else if ( ( emissionDOY_ >  31 ) && ( emissionDOY_ <=  59 ) ) {
        emissionMonth_ = 2;
        emissionDay_   = emissionDOY_ - 31;
    } else if ( ( emissionDOY_ >  59 ) && ( emissionDOY_ <=  90 ) ) {
        emissionMonth_ = 3;
        emissionDay_   = emissionDOY_ - 59;
    } else if ( ( emissionDOY_ >  90 ) && ( emissionDOY_ <= 120 ) ) {
        emissionMonth_ = 4;
        emissionDay_   = emissionDOY_ - 90;
    } else if ( ( emissionDOY_ > 120 ) && ( emissionDOY_ <= 151 ) ) {
        emissionMonth_ = 5;
        emissionDay_   = emissionDOY_ - 120;
    } else if ( ( emissionDOY_ > 151 ) && ( emissionDOY_ <= 181 ) ) {
        emissionMonth_ = 6;
        emissionDay_   = emissionDOY_ - 151;
    } else if ( ( emissionDOY_ > 181 ) && ( emissionDOY_ <= 212 ) ) {
        emissionMonth_ = 7;
        emissionDay_   = emissionDOY_ - 181;
    } else if ( ( emissionDOY_ > 212 ) && ( emissionDOY_ <= 243 ) ) {
        emissionMonth_ = 8;
        emissionDay_   = emissionDOY_ - 212;
    } else if ( ( emissionDOY_ > 243 ) && ( emissionDOY_ <= 273 ) ) {
        emissionMonth_ = 9;
        emissionDay_   = emissionDOY_ - 243;
    } else if ( ( emissionDOY_ > 273 ) && ( emissionDOY_ <= 304 ) ) {
        emissionMonth_ = 10;
        emissionDay_   = emissionDOY_ - 273;
    } else if ( ( emissionDOY_ > 304 ) && ( emissionDOY_ <= 334 ) ) {
        emissionMonth_ = 11;
        emissionDay_   = emissionDOY_ - 304;
    } else if ( ( emissionDOY_ > 334 ) && ( emissionDOY_ <= 365 ) ) {
        emissionMonth_ = 12;
        emissionDay_   = emissionDOY_ - 334;
    } else {
        std::cout << " emissionDOY = " << emissionDOY_ << std::endl;
        std::cout << " Could not figure out what month this is" << std::endl;
        exit(-1);
    }

} /* End of Input::Input */

Input::Input( unsigned int iCase,          \
              const Vector_2D &parameters, \
              const std::string fileName,  \
              const std::string fileName_ADJ ):
    temperature_K_ ( parameters[0][iCase] ),
    relHumidity_w_ ( parameters[1][iCase] ),
    longitude_deg_ ( parameters[2][iCase] ),
    latitude_deg_  ( parameters[3][iCase] ),
    pressure_Pa_   ( parameters[4][iCase] ),
    emissionDOY_   ( parameters[5][iCase] ),
    emissionTime_  ( parameters[6][iCase] ),
    EI_NOx_        ( parameters[7][iCase] ),
    EI_CO_         ( parameters[8][iCase] ),
    EI_HC_         ( parameters[9][iCase] ),
    EI_SO2_        ( parameters[10][iCase] ),
    EI_SO2TOSO4_   ( parameters[11][iCase] ),
    EI_Soot_       ( parameters[12][iCase] ),
    sootRad_       ( parameters[13][iCase] ),
    fuelFlow_      ( parameters[14][iCase] ),
    backgNOx_      ( parameters[15][iCase] ),
    backgHNO3_     ( parameters[16][iCase] ),
    backgO3_       ( parameters[17][iCase] ),
    backgCO_       ( parameters[18][iCase] ),
    backgCH4_      ( parameters[19][iCase] ),
    backgSO2_      ( parameters[20][iCase] ),
    fileName_      ( fileName ),
    fileName_ADJ_  ( fileName_ADJ )
{

    /* Constructor */
 
    while ( longitude_deg_ > 180 )
        longitude_deg_ -= 360;
    
    while ( longitude_deg_ < -180 )
        longitude_deg_ += 360;

    while ( latitude_deg_ > 90 )
        latitude_deg_ -= 180;
    
    while ( latitude_deg_ < -90 )
        latitude_deg_ += 180;

    while ( emissionDOY_ > 365 )
        emissionDOY_ -= 365;

    while ( emissionTime_ > 24.0 )
        emissionTime_ -= 24.0;

    if ( temperature_K_ >= 300.0 || temperature_K_ <= 160.0 ) {
        std::cout << " In Input::Input:";
        std::cout << " temperature_K takes an odd value: temperature_K = ";
        std::cout << temperature_K_ << " [K]" << std::endl;
        exit(-1);
    }
    
    if ( pressure_Pa_ >= 1.00E+05 || pressure_Pa_ <= 1.00E+03 ) {
        std::cout << " In Input::Input:";
        std::cout << " pressure_Pa takes an odd value: pressure_Pa_ = ";
        std::cout << pressure_Pa_ << " [Pa]" << std::endl;
        exit(-1);
    }

    if ( relHumidity_w_ >= 9.50E+01 || pressure_Pa_ <= 0.00E+00 ) {
        std::cout << " In Input::Input:";
        std::cout << " relHumidity_w takes an unrealisable value: relHumidity_w = ";
        std::cout << relHumidity_w_ << " [%]" << std::endl;
        exit(-1);
    }

    if ( emissionDOY_ <= 0.0 ) {
        std::cout << " In Input::Input:";
        std::cout << " emissionDOY takes an unrealisable value: emissionDOY = ";
        std::cout << emissionDOY_ << " [-]" << std::endl;
        exit(-1);
    }

    if ( EI_NOx_ < 0.0E+00 || EI_NOx_ > 5.0E+01 ) {
        std::cout << " In Input::Input:";
        std::cout << " EI_NOx takes an unrealisable value: EI_NOx = ";
        std::cout << EI_NOx_ << " [g/kg_fuel]" << std::endl;
        exit(-1);
    }
    
    if ( EI_CO_ < 0.0E+00 || EI_CO_ > 3.0E+01 ) {
        std::cout << " In Input::Input:";
        std::cout << " EI_CO takes an unrealisable value: EI_CO = ";
        std::cout << EI_CO_ << " [g/kg_fuel]" << std::endl;
        exit(-1);
    }
    
    if ( EI_HC_ < 0.0E+00 || EI_HC_ > 1.0E+01 ) {
        std::cout << " In Input::Input:";
        std::cout << " EI_HC takes an unrealisable value: EI_HC = ";
        std::cout << EI_HC_ << " [g/kg_fuel]" << std::endl;
        exit(-1);
    }
    
    if ( EI_SO2_ < 0.0E+00 || EI_SO2_ > 1.0E+02 ) {
        std::cout << " In Input::Input:";
        std::cout << " EI_SO2 takes an unrealisable value: EI_SO2 = ";
        std::cout << EI_SO2_ << " [g/kg_fuel]" << std::endl;
        exit(-1);
    }

    if ( EI_Soot_ < 0.0E+00 || EI_Soot_ > 2.0E-01 ) {
        std::cout << " In Input::Input:";
        std::cout << " EI_Soot takes an unrealisable value: EI_Soot = ";
        std::cout << EI_Soot_ << " [g/kg_fuel]" << std::endl;
        exit(-1);
    }
    
    if ( ( ( sootRad_ < 1.0E-10 ) && ( sootRad_ != 0.0E+00 ) ) || sootRad_ > 1.0E-07 ) {
        std::cout << " In Input::Input:";
        std::cout << " sootRad takes an unrealisable value: sootRad = ";
        std::cout << sootRad_ * 1.0E+09 << " [nm]" << std::endl;
        exit(-1);
    }
    
    if ( fuelFlow_ < 0.0E+00 ) {
        std::cout << " In Input::Input:";
        std::cout << " fuelFlow takes an unrealisable value: fuelFlow = ";
        std::cout << fuelFlow_ << " [kg/s]" << std::endl;
        exit(-1);
    }
    
    if ( backgNOx_ < 0.0E+00 || backgNOx_ > 1.0E+09 ) {
        std::cout << " In Input::Input:";
        std::cout << " backgNOx takes an unrealisable value: backgNOx = ";
        std::cout << backgNOx_ << " [ppb]" << std::endl;
        exit(-1);
    }
    
    if ( backgHNO3_ < 0.0E+00 || backgHNO3_ > 1.0E+09 ) {
        std::cout << " In Input::Input:";
        std::cout << " backgHNO3 takes an unrealisable value: backgHNO3 = ";
        std::cout << backgHNO3_ << " [ppb]" << std::endl;
        exit(-1);
    }
    
    if ( backgO3_ < 0.0E+00 || backgO3_ > 1.0E+09 ) {
        std::cout << " In Input::Input:";
        std::cout << " backgO3 takes an unrealisable value: backgO3 = ";
        std::cout << backgO3_ << " [ppb]" << std::endl;
        exit(-1);
    }
    
    if ( backgCO_ < 0.0E+00 || backgCO_ > 1.0E+09 ) {
        std::cout << " In Input::Input:";
        std::cout << " backgCO takes an unrealisable value: backgCO = ";
        std::cout << backgCO_ << " [ppb]" << std::endl;
        exit(-1);
    }
    
    if ( backgCH4_ < 0.0E+00 || backgCH4_ > 1.0E+09 ) {
        std::cout << " In Input::Input:";
        std::cout << " backgCH4 takes an unrealisable value: backgCH4 = ";
        std::cout << backgCH4_ << " [ppb]" << std::endl;
        exit(-1);
    }
    
    if ( backgSO2_ < 0.0E+00 || backgSO2_ > 1.0E+09 ) {
        std::cout << " In Input::Input:";
        std::cout << " backgSO2 takes an unrealisable value: backgCH4 = ";
        std::cout << backgSO2_ << " [ppb]" << std::endl;
        exit(-1);
    }

    if ( emissionDOY_ <= 31 ) {
        emissionMonth_ = 1;
        emissionDay_   = emissionDOY_;
    } else if ( ( emissionDOY_ >  31 ) && ( emissionDOY_ <=  59 ) ) {
        emissionMonth_ = 2;
        emissionDay_   = emissionDOY_ - 31;
    } else if ( ( emissionDOY_ >  59 ) && ( emissionDOY_ <=  90 ) ) {
        emissionMonth_ = 3;
        emissionDay_   = emissionDOY_ - 59;
    } else if ( ( emissionDOY_ >  90 ) && ( emissionDOY_ <= 120 ) ) {
        emissionMonth_ = 4;
        emissionDay_   = emissionDOY_ - 90;
    } else if ( ( emissionDOY_ > 120 ) && ( emissionDOY_ <= 151 ) ) {
        emissionMonth_ = 5;
        emissionDay_   = emissionDOY_ - 120;
    } else if ( ( emissionDOY_ > 151 ) && ( emissionDOY_ <= 181 ) ) {
        emissionMonth_ = 6;
        emissionDay_   = emissionDOY_ - 151;
    } else if ( ( emissionDOY_ > 181 ) && ( emissionDOY_ <= 212 ) ) {
        emissionMonth_ = 7;
        emissionDay_   = emissionDOY_ - 181;
    } else if ( ( emissionDOY_ > 212 ) && ( emissionDOY_ <= 243 ) ) {
        emissionMonth_ = 8;
        emissionDay_   = emissionDOY_ - 212;
    } else if ( ( emissionDOY_ > 243 ) && ( emissionDOY_ <= 273 ) ) {
        emissionMonth_ = 9;
        emissionDay_   = emissionDOY_ - 243;
    } else if ( ( emissionDOY_ > 273 ) && ( emissionDOY_ <= 304 ) ) {
        emissionMonth_ = 10;
        emissionDay_   = emissionDOY_ - 273;
    } else if ( ( emissionDOY_ > 304 ) && ( emissionDOY_ <= 334 ) ) {
        emissionMonth_ = 11;
        emissionDay_   = emissionDOY_ - 304;
    } else if ( ( emissionDOY_ > 334 ) && ( emissionDOY_ <= 365 ) ) {
        emissionMonth_ = 12;
        emissionDay_   = emissionDOY_ - 334;
    } else {
        std::cout << " emissionDOY = " << emissionDOY_ << std::endl;
        std::cout << " Could not figure out what month this is" << std::endl;
        exit(-1);
    }


} /* End of Input::Input */

Input::~Input()
{

    /* Destructor */

} /* End of Input::~Input */

/* End of Input.cpp */
