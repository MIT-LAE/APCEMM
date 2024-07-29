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

Input::Input( unsigned int iCase,               \
              const Vector_2D &parameters,      \
              const std::string fileName,       \
              const std::string fileName_ADJ,   \
              const std::string fileName_BOX,   \
              const std::string fileName_micro, \
              const std::string author          ):
    Case_          ( iCase                 ),
    simulationTime_( parameters[ 0][iCase] ),
    temperature_K_ ( parameters[ 1][iCase] ),
    relHumidity_w_ ( parameters[ 2][iCase] ),
    horizDiff_     ( parameters[ 3][iCase] ),
    vertiDiff_     ( parameters[ 4][iCase] ),
    shear_         ( parameters[ 5][iCase] ),
    longitude_deg_ ( parameters[ 6][iCase] ),
    latitude_deg_  ( parameters[ 7][iCase] ),
    pressure_Pa_   ( parameters[ 8][iCase] ),
    emissionDOY_   ( parameters[ 9][iCase] ),
    emissionTime_  ( parameters[10][iCase] ),
    EI_NOx_        ( parameters[11][iCase] ),
    EI_CO_         ( parameters[12][iCase] ),
    EI_HC_         ( parameters[13][iCase] ),
    EI_SO2_        ( parameters[14][iCase] ),
    EI_SO2TOSO4_   ( parameters[15][iCase] ),
    EI_Soot_       ( parameters[16][iCase] ),
    sootRad_       ( parameters[17][iCase] ),
    fuelFlow_      ( parameters[18][iCase] ),
    aircraftMass_  ( parameters[19][iCase] ),
    backgNOx_      ( parameters[20][iCase] ),
    backgHNO3_     ( parameters[21][iCase] ),
    backgO3_       ( parameters[22][iCase] ),
    backgCO_       ( parameters[23][iCase] ),
    backgCH4_      ( parameters[24][iCase] ),
    backgSO2_      ( parameters[25][iCase] ),
    flightSpeed_   ( parameters[26][iCase] ),
    numEngines_    ( parameters[27][iCase] ),
    wingspan_      ( parameters[28][iCase] ),
    coreExitTemp_  ( parameters[29][iCase] ),
    bypassArea_    ( parameters[30][iCase] ),
    fileName_      ( fileName ),
    fileName_ADJ_  ( fileName_ADJ ),
    fileName_BOX_  ( fileName_BOX ),
    fileName_micro_ ( fileName_micro ),
    author_        ( author )
{

    /* Constructor */
    adjustLatLong();
    checkInputValidity();
    calculate_emissionMonth();

} /* End of Input::Input */

Input::Input( unsigned int iCase,               \
        const std::vector<std::unordered_map<std::string, double>> &parameters,      \
        const std::string fileName,       \
        const std::string fileName_ADJ,   \
        const std::string fileName_BOX,   \
        const std::string fileName_micro, \
        const std::string author          ):
    Case_          ( iCase                 ),
    simulationTime_( parameters[iCase].at("PLUMEPROCESS")),
    temperature_K_ ( parameters[iCase].at("TEMPERATURE")),
    relHumidity_w_ ( parameters[iCase].at("RHW")),
    horizDiff_     ( parameters[iCase].at("DH")),
    vertiDiff_     ( parameters[iCase].at("DV")),
    shear_         ( parameters[iCase].at("SHEAR")),

    longitude_deg_ ( parameters[iCase].at("LONGITUDE")),
    latitude_deg_  ( parameters[iCase].at("LATITUDE")),
    pressure_Pa_   ( parameters[iCase].at("PRESSURE")),

    emissionDOY_   ( parameters[iCase].at("EDAY")),
    emissionTime_  ( parameters[iCase].at("ETIME")),

    EI_NOx_        ( parameters[iCase].at("EI_NOX")),
    EI_CO_         ( parameters[iCase].at("EI_CO")),
    EI_HC_         ( parameters[iCase].at("EI_UHC")),
    EI_SO2_        ( parameters[iCase].at("EI_SO2")),
    EI_SO2TOSO4_   ( parameters[iCase].at("EI_SO2TOSO4")),
    EI_Soot_       ( parameters[iCase].at("EI_SOOT")),
    sootRad_       ( parameters[iCase].at("EI_SOOTRAD")),

    fuelFlow_      ( parameters[iCase].at("FF")),
    aircraftMass_  ( parameters[iCase].at("AMASS")),

    backgNOx_      ( parameters[iCase].at("BACKG_NOX")),
    backgHNO3_     ( parameters[iCase].at("BACKG_HNO3")),
    backgO3_       ( parameters[iCase].at("BACKG_O3")),
    backgCO_       ( parameters[iCase].at("BACKG_CO")),
    backgCH4_      ( parameters[iCase].at("BACKG_CH4")),
    backgSO2_      ( parameters[iCase].at("BACKG_SO2")),

    flightSpeed_   ( parameters[iCase].at("FSPEED")),
    numEngines_    ( parameters[iCase].at("NUMENG")),
    wingspan_      ( parameters[iCase].at("WINGSPAN")),
    coreExitTemp_  ( parameters[iCase].at("COREEXITTEMP")),
    bypassArea_    ( parameters[iCase].at("BYPASSAREA")),
    fileName_      ( fileName ),
    fileName_ADJ_  ( fileName_ADJ ),
    fileName_BOX_  ( fileName_BOX ),
    fileName_micro_ ( fileName_micro ),
    author_        ( author ),

    nBV_           ( parameters[iCase].at("NBV"))

{

}
Input::~Input()
{

    /* Destructor */

} /* End of Input::~Input */

/* End of Input.cpp */

void Input::checkInputValidity(){
        if ( simulationTime_ <= 0.0 || simulationTime_ >= 48.0 ) {
        std::cout << " In Input::Input:";
        std::cout << " simulationTime takes an odd value: simulationTime = ";
        std::cout << simulationTime_ << " [hrs]" << std::endl;
        exit(-1);
    }

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

    if ( relHumidity_w_ >= 9.50E+01 || relHumidity_w_ <= 0.00E+00 ) {
        std::cout << " In Input::Input:";
        std::cout << " relHumidity_w takes an unrealisable value: relHumidity_w = ";
        std::cout << relHumidity_w_ << " [%]" << std::endl;
        exit(-1);
    }
    
    if ( horizDiff_ >= 4.00E+01 || horizDiff_ < 1.00E+00 ) {
        std::cout << " In Input::Input:";
        std::cout << " horizDiff takes an odd value: horizDiff_ = ";
        std::cout << horizDiff_ << " [m^2/s]" << std::endl;
        if ( horizDiff_ < 0.00E+00 )
            exit(-1);
    }

    if ( vertiDiff_ >= 4.00E-01 || vertiDiff_ < 1.00E-02 ) {
        std::cout << " In Input::Input:";
        std::cout << " vertiDiff takes an odd value: vertiDiff_ = ";
        std::cout << vertiDiff_ << " [m^2/s]" << std::endl;
        if ( vertiDiff_ < 0.00E+00 )
            exit(-1);
    }

    if ( shear_ >= 5.0E-02 || shear_ <= -5.0E-02 ) {
        std::cout << " In Input::Input:";
        std::cout << " shear takes an unrealisable value: shear = ";
        std::cout << shear_ << " [1/s]" << std::endl;
        exit(-1);
    }

    if ( emissionDOY_ < 1  || emissionDOY_ > 365) {
        std::cout << " In Input::Input:";
        std::cout << " emissionDOY takes an unrealisable value: emissionDOY = ";
        std::cout << emissionDOY_ << " [-]" << std::endl;
        exit(-1);
    }

    if ( emissionTime_ < 0 || emissionTime_ > 24) {
        std::cout << " In Input::Input:";
        std::cout << " emissionTime takes an unrealisable value: emissionDOY = ";
        std::cout << emissionTime_ << " [-]" << std::endl;
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

    /* if ( aircraftMass_ < 50.0E+03 ) { */
    if ( aircraftMass_ < 0.0E+00 ) {
        std::cout << " In Input::Input:";
        std::cout << " aircraftMass takes an unrealisable value: aircraftMass = ";
        std::cout << aircraftMass_ << " [kg]" << std::endl;
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

    if ( flightSpeed_ < 0.0E+00 ) {
        std::cout << " In Input::Input:";
        std::cout << " Flight speed takes an unrealisable value: flightSpeed = ";
        std::cout << flightSpeed_ << " [m/s]" << std::endl;
        exit(-1);
    }
    
    if ( !( numEngines_ == 2 || numEngines_ == 3 ||  numEngines_ == 4 ) ) {
        std::cout << " In Input::Input:";
        std::cout << " Number of engines takes an unrealisable value: numEngines = ";
        std::cout << numEngines_ << " []" << std::endl;
        exit(-1);
    }

    if ( wingspan_ < 0.0E+00 ) {
        std::cout << " In Input::Input:";
        std::cout << " Wingspan takes an unrealisable value: wingspan_ = ";
        std::cout << wingspan_ << " [m]" << std::endl;
        exit(-1);
    }

    if ( coreExitTemp_ < 0.0E+00 ) {
        std::cout << " In Input::Input:";
        std::cout << " Exhaust temperature takes an unrealisable value: coreExitTemp = ";
        std::cout << coreExitTemp_ << " [K]" << std::endl;
        exit(-1);
    }

    if ( bypassArea_ < 0.0E+00 ) {
        std::cout << " In Input::Input:";
        std::cout << " Exhaust area takes an unrealisable value: bypassArea = ";
        std::cout << bypassArea_ << " [m^2]" << std::endl;
        exit(-1);
    }
}

void Input::adjustLatLong(){
    while ( longitude_deg_ > 180 )
        longitude_deg_ -= 360;
    
    while ( longitude_deg_ < -180 )
        longitude_deg_ += 360;

    while ( latitude_deg_ > 90 )
        latitude_deg_ -= 180;
    
    while ( latitude_deg_ < -90 )
        latitude_deg_ += 180;
}

void Input::calculate_emissionMonth(){
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
}
