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
    dayGMT_        ( dayGMT ),
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

    while ( dayGMT_ > 365 )
        dayGMT_ -= 365;

    while ( emissionTime_ > 24.0 )
        emissionTime_ -= 24.0;

    if ( temperature_K_ >= 300.0 || temperature_K_ <= 160.0 )
        std::cout << " In Input::Input: temperature_K takes an odd value: temperature_K = " << temperature_K_ << " [K]\n";
    
    if ( pressure_Pa_ >= 8.00E+04 || pressure_Pa_ <= 2.00E+03 )
        std::cout << " In Input::Input: pressure_Pa takes an odd value: pressure_Pa_ = " << pressure_Pa_/100.0 << " [hPa]\n";

    if ( relHumidity_w_ >= 9.50E+01 || pressure_Pa_ <= 0.00E+00 )
        std::cout << " In Input::Input: relHumidity_w takes an unrealisable value: relHumidity_w = " << relHumidity_w_ << " [%]\n";

    if ( dayGMT_ <= 0.0 )
        std::cout << " In Input::Input: dayGMT takes an unrealisable value: dayGMT = " << dayGMT_ << " [-]\n";

    if ( EI_NOx_ < 0.0E+00 || EI_NOx_ > 3.0E+01 )
        std::cout << " In Input::Input: EI_NOx takes an unrealisable value: EI_NOx = " << EI_NOx_ << " [g/kg_fuel]\n";
    
    if ( EI_CO_ < 0.0E+00 || EI_CO_ > 3.0E+01 )
        std::cout << " In Input::Input: EI_CO takes an unrealisable value: EI_CO = " << EI_CO_ << " [g/kg_fuel]\n";
    
    if ( EI_HC_ < 0.0E+00 || EI_HC_ > 1.0E+01 )
        std::cout << " In Input::Input: EI_HC takes an unrealisable value: EI_HC = " << EI_HC_ << " [g/kg_fuel]\n";
    
    if ( EI_SO2_ < 0.0E+00 || EI_SO2_ > 1.0E+02 )
        std::cout << " In Input::Input: EI_SO2 takes an unrealisable value: EI_SO2 = " << EI_SO2_ << " [g/kg_fuel]\n";
    
    if ( EI_Soot_ < 0.0E+00 || EI_Soot_ > 2.0E-01 )
        std::cout << " In Input::Input: EI_Soot takes an unrealisable value: EI_Soot = " << EI_Soot_ << " [g/kg_fuel]\n";
    
    if ( ( ( sootRad_ < 1.0E-10 ) && ( sootRad_ != 0.0E+00 ) ) || sootRad_ > 1.0E-07 )
        std::cout << " In Input::Input: sootRad takes an unrealisable value: sootRad = " << sootRad_ * 1.0E+09 << " [nm]\n";
    
    if ( fuelFlow_ < 0.0E+00 )
        std::cout << " In Input::Input: fuelFlow takes an unrealisable value: fuelFlow = " << fuelFlow_ << " [kg/s/engine]\n";
    
    if ( backgNOx_ < 0.0E+00 || backgNOx_ > 1.0E+09 )
        std::cout << " In Input::Input: backgNOx takes an unrealisable value: backgNOx = " << backgNOx_ << " [ppb]\n";
    
    if ( backgHNO3_ < 0.0E+00 || backgHNO3_ > 1.0E+09 )
        std::cout << " In Input::Input: backgHNO3 takes an unrealisable value: backgHNO3 = " << backgHNO3_ << " [ppb]\n";
    
    if ( backgO3_ < 0.0E+00 || backgO3_ > 1.0E+09 )
        std::cout << " In Input::Input: backgO3 takes an unrealisable value: backgO3 = " << backgO3_ << " [ppb]\n";
    
    if ( backgCO_ < 0.0E+00 || backgCO_ > 1.0E+09 )
        std::cout << " In Input::Input: backgCO takes an unrealisable value: backgCO = " << backgCO_ << " [ppb]\n";
    
    if ( backgCH4_ < 0.0E+00 || backgCH4_ > 1.0E+09 )
        std::cout << " In Input::Input: backgCH4 takes an unrealisable value: backgCH4 = " << backgCH4_ << " [ppb]\n";
    
    if ( backgSO2_ < 0.0E+00 || backgSO2_ > 1.0E+09 )
        std::cout << " In Input::Input: backgSO2 takes an unrealisable value: backgCH4 = " << backgSO2_ << " [ppb]\n";

} /* End of Input::Input */

Input::Input( unsigned int iCase,          \
              const Vector_2D &parameters, \
              const std::string fileName,  \
              const std::string fileName_ADJ ):
    temperature_K_ ( parameters[0][iCase] ),
    pressure_Pa_   ( parameters[1][iCase] ),
    relHumidity_w_ ( parameters[2][iCase] ),
    longitude_deg_ ( parameters[3][iCase] ),
    latitude_deg_  ( parameters[4][iCase] ),
    dayGMT_        ( parameters[5][iCase] ),
    emissionTime_  ( parameters[6][iCase] ),
    EI_NOx_        ( parameters[7][iCase] ),
    EI_CO_         ( parameters[8][iCase] ),
    EI_HC_         ( parameters[9][iCase] ),
    EI_SO2_        ( parameters[10][iCase] ),
    EI_Soot_       ( parameters[11][iCase] ),
    sootRad_       ( parameters[12][iCase] ),
    fuelFlow_      ( parameters[13][iCase] ),
    backgNOx_      ( parameters[14][iCase] ),
    backgHNO3_     ( parameters[15][iCase] ),
    backgO3_       ( parameters[16][iCase] ),
    backgCO_       ( parameters[17][iCase] ),
    backgCH4_      ( parameters[18][iCase] ),
    backgSO2_      ( parameters[19][iCase] ),
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

    while ( dayGMT_ > 365 )
        dayGMT_ -= 365;

    while ( emissionTime_ > 24.0 )
        emissionTime_ -= 24.0;

    if ( temperature_K_ >= 300.0 || temperature_K_ <= 160.0 )
        std::cout << " In Input::Input: temperature_K takes an odd value: temperature_K = " << temperature_K_ << " [K]\n";
    
    if ( pressure_Pa_ >= 8.00E+04 || pressure_Pa_ <= 2.00E+03 )
        std::cout << " In Input::Input: pressure_Pa takes an odd value: pressure_Pa_ = " << pressure_Pa_/100.0 << " [hPa]\n";

    if ( relHumidity_w_ >= 9.50E+01 || pressure_Pa_ <= 0.00E+00 )
        std::cout << " In Input::Input: relHumidity_w takes an unrealisable value: relHumidity_w = " << relHumidity_w_ << " [%]\n";

    if ( dayGMT_ <= 0.0 )
        std::cout << " In Input::Input: dayGMT takes an unrealisable value: dayGMT = " << dayGMT_ << " [-]\n";

    if ( EI_NOx_ < 0.0E+00 || EI_NOx_ > 3.0E+01 )
        std::cout << " In Input::Input: EI_NOx takes an unrealisable value: EI_NOx = " << EI_NOx_ << " [g/kg_fuel]\n";
    
    if ( EI_CO_ < 0.0E+00 || EI_CO_ > 3.0E+01 )
        std::cout << " In Input::Input: EI_CO takes an unrealisable value: EI_CO = " << EI_CO_ << " [g/kg_fuel]\n";
    
    if ( EI_HC_ < 0.0E+00 || EI_HC_ > 1.0E+01 )
        std::cout << " In Input::Input: EI_HC takes an unrealisable value: EI_HC = " << EI_HC_ << " [g/kg_fuel]\n";
    
    if ( EI_SO2_ < 0.0E+00 || EI_SO2_ > 1.0E+02 )
        std::cout << " In Input::Input: EI_SO2 takes an unrealisable value: EI_SO2 = " << EI_SO2_ << " [g/kg_fuel]\n";

    if ( EI_Soot_ < 0.0E+00 || EI_Soot_ > 2.0E-01 )
        std::cout << " In Input::Input: EI_Soot takes an unrealisable value: EI_Soot = " << EI_Soot_ << " [g/kg_fuel]\n";
    
    if ( ( ( sootRad_ < 1.0E-10 ) && ( sootRad_ != 0.0E+00 ) ) || sootRad_ > 1.0E-07 )
        std::cout << " In Input::Input: sootRad takes an unrealisable value: sootRad = " << sootRad_ * 1.0E+09 << " [nm]\n";
    
    if ( fuelFlow_ < 0.0E+00 )
        std::cout << " In Input::Input: fuelFlow takes an unrealisable value: fuelFlow = " << fuelFlow_ << " [kg/s/engine]\n";
    
    if ( backgNOx_ < 0.0E+00 || backgNOx_ > 1.0E+09 )
        std::cout << " In Input::Input: backgNOx takes an unrealisable value: backgNOx = " << backgNOx_ << " [ppb]\n";
    
    if ( backgHNO3_ < 0.0E+00 || backgHNO3_ > 1.0E+09 )
        std::cout << " In Input::Input: backgHNO3 takes an unrealisable value: backgHNO3 = " << backgHNO3_ << " [ppb]\n";
    
    if ( backgO3_ < 0.0E+00 || backgO3_ > 1.0E+09 )
        std::cout << " In Input::Input: backgO3 takes an unrealisable value: backgO3 = " << backgO3_ << " [ppb]\n";
    
    if ( backgCO_ < 0.0E+00 || backgCO_ > 1.0E+09 )
        std::cout << " In Input::Input: backgCO takes an unrealisable value: backgCO = " << backgCO_ << " [ppb]\n";
    
    if ( backgCH4_ < 0.0E+00 || backgCH4_ > 1.0E+09 )
        std::cout << " In Input::Input: backgCH4 takes an unrealisable value: backgCH4 = " << backgCH4_ << " [ppb]\n";
    
    if ( backgSO2_ < 0.0E+00 || backgSO2_ > 1.0E+09 )
        std::cout << " In Input::Input: backgSO2 takes an unrealisable value: backgCH4 = " << backgSO2_ << " [ppb]\n";

} /* End of Input::Input */

Input::~Input()
{

    /* Destructor */

} /* End of Input::~Input */

/* End of Input.cpp */
