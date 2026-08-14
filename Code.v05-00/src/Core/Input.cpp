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

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "Core/Input.hpp"

void Input::checkInputValidity(){

    /* Collect all error messages before throwing */
    std::vector<std::string> errors;

    /* Value is outside of normal range, report and crash the run */
    auto reject = [&errors]( const char* name, auto value,
                             const char* unit, const char* expected ){
        std::ostringstream os;
        os << name << " = " << value << " [" << unit << "]"
           << ", expected " << expected;
        errors.push_back( os.str() );
    };

    /* Value is unusual, so we warn */
    auto warn = []( const char* name, auto value,
                    const char* unit, const char* typical ){
        std::cout << " In Input::checkInputValidity: " << name << " = " << value
                  << " [" << unit << "] is outside the typical range "
                  << typical << std::endl;
    };

    if ( simulationTime_ <= 0.0 || simulationTime_ >= 48.0 )
        reject( "simulationTime", simulationTime_, "hrs", "0 < simulationTime < 48" );

    if ( pressure_Pa_ >= 1.00E+05 || pressure_Pa_ <= 1.00E+03 )
        reject( "pressure_Pa", pressure_Pa_, "Pa", "1.0E+03 < pressure_Pa < 1.0E+05" );

    if ( horizDiff_ < 0.00E+00 )
        reject( "horizDiff", horizDiff_, "m^2/s", "horizDiff >= 0" );
    else if ( horizDiff_ >= 4.00E+01 || horizDiff_ < 1.00E+00 )
        warn( "horizDiff", horizDiff_, "m^2/s", "1 to 40" );

    if ( vertiDiff_ < 0.00E+00 )
        reject( "vertiDiff", vertiDiff_, "m^2/s", "vertiDiff >= 0" );
    else if ( vertiDiff_ >= 4.00E-01 || vertiDiff_ < 1.00E-02 )
        warn( "vertiDiff", vertiDiff_, "m^2/s", "0.01 to 0.4" );

    if ( emissionDOY_ < 1 || emissionDOY_ > 365 )
        reject( "emissionDOY", emissionDOY_, "-", "1 to 365" );

    if ( emissionTime_ < 0 || emissionTime_ > 24 )
        reject( "emissionTime", emissionTime_, "hrs", "0 to 24" );

    if ( EI_NOx_ < 0.0E+00 || EI_NOx_ > 5.0E+01 )
        reject( "EI_NOx", EI_NOx_, "g/kg_fuel", "0 to 50" );

    if ( EI_CO_ < 0.0E+00 || EI_CO_ > 3.0E+01 )
        reject( "EI_CO", EI_CO_, "g/kg_fuel", "0 to 30" );

    if ( EI_HC_ < 0.0E+00 || EI_HC_ > 1.0E+01 )
        reject( "EI_HC", EI_HC_, "g/kg_fuel", "0 to 10" );

    if ( EI_SO2_ < 0.0E+00 || EI_SO2_ > 1.0E+02 )
        reject( "EI_SO2", EI_SO2_, "g/kg_fuel", "0 to 100" );

    if ( EI_Soot_ < 0.0E+00 || EI_Soot_ > 2.0E-01 )
        reject( "EI_Soot", EI_Soot_, "g/kg_fuel", "0 to 0.2" );

    if ( ( ( sootRad_ < 1.0E-10 ) && ( sootRad_ != 0.0E+00 ) ) || sootRad_ > 1.0E-07 )
        reject( "sootRad", sootRad_ * 1.0E+09, "nm", "0, or 0.1 to 100" );

    if ( fuelFlow_ < 0.0E+00 )
        reject( "fuelFlow", fuelFlow_, "kg/s", "fuelFlow >= 0" );

    if ( aircraftMass_ < 0.0E+00 )
        reject( "aircraftMass", aircraftMass_, "kg", "aircraftMass >= 0" );

    if ( backgNOx_ < 0.0E+00 || backgNOx_ > 1.0E+09 )
        reject( "backgNOx", backgNOx_, "ppb", "0 to 1.0E+09" );

    if ( backgHNO3_ < 0.0E+00 || backgHNO3_ > 1.0E+09 )
        reject( "backgHNO3", backgHNO3_, "ppb", "0 to 1.0E+09" );

    if ( backgO3_ < 0.0E+00 || backgO3_ > 1.0E+09 )
        reject( "backgO3", backgO3_, "ppb", "0 to 1.0E+09" );

    if ( backgCO_ < 0.0E+00 || backgCO_ > 1.0E+09 )
        reject( "backgCO", backgCO_, "ppb", "0 to 1.0E+09" );

    if ( backgCH4_ < 0.0E+00 || backgCH4_ > 1.0E+09 )
        reject( "backgCH4", backgCH4_, "ppb", "0 to 1.0E+09" );

    if ( backgSO2_ < 0.0E+00 || backgSO2_ > 1.0E+09 )
        reject( "backgSO2", backgSO2_, "ppb", "0 to 1.0E+09" );

    if ( flightSpeed_ < 0.0E+00 )
        reject( "flightSpeed", flightSpeed_, "m/s", "flightSpeed >= 0" );

    if ( !( numEngines_ == 2 || numEngines_ == 3 || numEngines_ == 4 ) )
        reject( "numEngines", numEngines_, "-", "2, 3 or 4" );

    if ( wingspan_ < 0.0E+00 )
        reject( "wingspan", wingspan_, "m", "wingspan >= 0" );

    if ( coreExitTemp_ < 0.0E+00 )
        reject( "coreExitTemp", coreExitTemp_, "K", "coreExitTemp >= 0" );

    if ( bypassArea_ < 0.0E+00 )
        reject( "bypassArea", bypassArea_, "m^2", "bypassArea >= 0" );

    if ( !errors.empty() ) {
        std::ostringstream os;
        os << "Invalid value(s) in the PARAMETER MENU:";
        for ( const std::string& error : errors )
            os << "\n    - " << error;
        throw std::invalid_argument( os.str() );
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
