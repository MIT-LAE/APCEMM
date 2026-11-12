/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Emission Program File                                            */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Emission.cpp                              */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <iomanip>
#include "Core/Emission.hpp"

Emission::Emission( const Engine &engine, const Fuel &fuel, const double SO2ToSO4Ratio )
{

    /* Constructor */

    Populate_withEngine( engine );

    Populate_withFuel( fuel );
    
    /* Carbon footprint */
    const double carbContent = 0.95; 
    CO2 -= ( 44.095 / 28.010 * CO + 44.095 / (82.0/5.0) * HC + 44.095 / 12.0 * carbContent * Soot);
    /* Assume that soot is 95% carbon by weight */
    /* Add organic carbon .. TBD */

    /* HC splitting */
    CH4  = 0.268248 * HC;
    C2H6 = 0.013750 * HC;
    PRPE = 0.384804 * HC;
    ALK4 = 0.013497 * HC;
    CH2O = 0.184800 * HC;
    ALD2 = 0.086150 * HC;
    GLYX = 0.017129 * HC;
    MGLY = 0.009899 * HC;

    if ( CO2 <= 0.0 ) {
        std::cout << "CO2 emissions are negative: EI(CO_2) = " << CO2 << "g/kg" << std::endl;
        CO2 = 0.0;
    }
    if ( H2O <= 0.0 ) {
        std::cout << "H2O emissions are negative: EI(H_2O) = " << H2O << "g/kg" << std::endl;
        H2O = 0.0;
    }
    if ( SO2 < 0.0 ) {
        std::cout << "SO2 emissions are negative: EI(SO_2) = " << SO2 << "g/kg" << std::endl;
        SO2 = 0.0;
    }
    if ( NOx <= 0.0 ) {
        std::cout << "NOx emissions are negative: EI(NO_x) = " << NOx << "g/kg" << std::endl;
        NOx = 0.0;
    }
    if ( CO <= 0.0 ) {
        std::cout << "CO emissions are negative: EI(CO) = " << CO << "g/kg" << std::endl;
        CO = 0.0;
    }
    if ( HC <= 0.0 ) {
        std::cout << "HC emissions are negative: EI(HC) = " << HC << "g/kg" << std::endl;
        HC = 0.0;
    }

    if ( Soot < 0.0 ) {
        std::cout << "Soot emissions are negative: EI(Soot) = " << Soot << "g/kg" << std::endl;
        Soot = 0.0;
    }

    SO2toSO4 = SO2ToSO4Ratio;

    engineName = engine.getName();
    fuelChem   = fuel.getChemFormula();

} /* End of Emission::Emission */

void Emission::Populate_withEngine( const Engine &engine )
{

    /* Gaseous species */
    /* Engine characteristics */
    NOx  = engine.getEI_NOx();  /* g(NO2) /kg */
    NO   = engine.getEI_NO();   /* g(NO)  /kg */
    NO2  = engine.getEI_NO2();  /* g(NO2) /kg */
    HNO2 = engine.getEI_HNO2(); /* g(HNO2)/kg */
    CO   = engine.getEI_CO();
    HC   = engine.getEI_HC();

    /* BC */
    Soot = engine.getEI_Soot();
    SootRad = engine.getSootRad();

} /* End of Emission::Populate_withEngine */

void Emission::Populate_withFuel( const Fuel &fuel )
{
    /*** Fuel characteristics */
    CO2 = ( 12.0107 * fuel.getAtomC() ) / ( 12.0107 * fuel.getAtomC() + 1.0079 * fuel.getAtomH() + 10.0067 * fuel.getAtomN() + 32.065 * fuel.getAtomS() ) * ( 44.0095 / ( 1 * 12.0107 ) ); /* [ g/g fuel ] */
    H2O = EI_H2O( fuel ); /* [ g/kg fuel ] */
    /* Convert to g/kg */
    CO2 *= 1000; /* [ g/kg fuel ]*/
    SO2 = fuel.EI_SO2(); /* [ g/kg fuel ] */

} /* End of Emission::Populate_withFuel */

double Emission::EI_H2O( const Fuel &fuel )
{

    const double EI = (  1.0079 * fuel.getAtomH() ) / ( 12.0107 * fuel.getAtomC() + 1.0079 * fuel.getAtomH() + 10.0067 * fuel.getAtomN() + 32.065 * fuel.getAtomS() ) * ( 18.0153 / ( 2 *  1.0079 ) ); /* [ g/g fuel ] */
    return EI * 1000; /* [ g/kg fuel ] */

} /* End of Emission::EI_H2O */

double Emission::getCO2( ) const
{

    return CO2;

} /* End of Emission::getCO2 */

double Emission::getH2O( ) const
{

    return H2O;

} /* End of Emission::getH2O */

double Emission::getNOx( ) const
{

    return NOx;

} /* End of Emission::getNOx */

double Emission::getNO( ) const
{

    return NO;

} /* End of Emission::getNO */

double Emission::getNO2( ) const
{

    return NO2;

} /* End of Emission::getNO2 */

double Emission::getHNO2( ) const
{

    return HNO2;

} /* End of Emission::getHNO2 */

double Emission::getSO2( ) const
{

    return SO2;

} /* End of Emission::getSO2 */

double Emission::getCO( ) const
{

    return CO;

} /* End of Emission::getCO */

double Emission::getHC( ) const
{

    return HC;

} /* End of Emission::getHC */

double Emission::getCH4( ) const
{

    return CH4;

} /* End of Emission::getCH4 */

double Emission::getC2H6( ) const
{

    return C2H6;

} /* End of Emission::getC2H6 */

double Emission::getPRPE( ) const
{

    return PRPE;

} /* End of Emission::getPRPE */

double Emission::getALK4( ) const
{

    return ALK4;

} /* End of Emission::getALK4 */

double Emission::getCH2O( ) const
{

    return CH2O;

} /* End of Emission::getCH2O */

double Emission::getALD2( ) const
{

    return ALD2;

} /* End of Emission::getALD2 */

double Emission::getGLYX( ) const
{

    return GLYX;

} /* End of Emission::getGLYX */

double Emission::getMGLY( ) const
{

    return MGLY;

} /* End of Emission::getMGLY */

double Emission::getSO2toSO4( ) const
{

    return SO2toSO4;
    
} /* End of Emission::getSO2toSO4 */

double Emission::getSoot( ) const
{

    return Soot;

} /* End of Emission::getSoot */

double Emission::getSootRad( ) const
{

    return SootRad;

} /* End of Emission::getSootRad */

std::string Emission::getEngineName( ) const
{

    return engineName;

} /* End of Emission::getEngineName */

std::string Emission::getFuelChem( ) const
{

    return fuelChem;

} /* End of Emission::getFuelChem */


void Emission::Debug( ) const
{
    std::cout << std::endl;
    std::cout << "**** Input Debugger ****" << std::endl;
    std::cout << "Emission indices: " << std::endl;
    std::cout << std::endl;
    std::cout << std::setw(10);
    std::cout << "Species";
    std::cout << std::setw(13);
    std::cout << "Value";
    std::cout << std::setw(14);
    std::cout << "Units " << std::endl;
    std::cout << std::setw(12);
    std::cout << " -> CO2 : ";
    std::cout << std::setw(12);
    std::cout << CO2 << "   [ g/kg fuel ] " << std::endl;
    std::cout << std::setw(12);
    std::cout << " -> H2O : ";
    std::cout << std::setw(12);
    std::cout << H2O << "   [ g/kg fuel ] " << std::endl;
    std::cout << std::setw(12);
    std::cout << " -> NOx : ";
    std::cout << std::setw(12);
    std::cout << NOx << "   [ g/kg fuel ] " << std::endl;
    std::cout << std::setw(14);
    std::cout << " ---> NO  : ";
    std::cout << std::setw(10);
    std::cout << NO << "   [ g/kg fuel ] " << std::endl;
    std::cout << std::setw(14);
    std::cout << " ---> NO2 : ";
    std::cout << std::setw(10);
    std::cout << NO2 << "   [ g/kg fuel ] " << std::endl;
    std::cout << std::setw(14);
    std::cout << " ---> HNO2: ";
    std::cout << std::setw(10);
    std::cout << HNO2 << "   [ g/kg fuel ] " << std::endl;
    std::cout << std::setw(12);
    std::cout << " -> SO2 : ";
    std::cout << std::setw(12);
    std::cout << SO2 << "   [ g/kg fuel ] " << std::endl;
    std::cout << std::setw(12);
    std::cout << " -> CO  : ";
    std::cout << std::setw(12);
    std::cout << CO << "   [ g/kg fuel ] " << std::endl;
    std::cout << std::setw(12);
    std::cout << " -> HC  : ";
    std::cout << std::setw(12);
    std::cout << HC << "   [ g/kg fuel ] " << std::endl;
    std::cout << std::setw(14);
    std::cout << " ---> CH4 : ";
    std::cout << std::setw(10);
    std::cout << CH4 << "   [ g/kg fuel ] " << std::endl;
    std::cout << std::setw(14);
    std::cout << " ---> C2H6: ";
    std::cout << std::setw(10);
    std::cout << C2H6 << "  [ g/kg fuel ] " << std::endl;
    std::cout << std::setw(14);
    std::cout << " ---> PRPE: ";
    std::cout << std::setw(10);
    std::cout << PRPE << "   [ g/kg fuel ] " << std::endl;
    std::cout << std::setw(14);
    std::cout << " ---> ALK4: ";
    std::cout << std::setw(10);
    std::cout << ALK4 << "  [ g/kg fuel ] " << std::endl;
    std::cout << std::setw(14);
    std::cout << " ---> GLYX: ";
    std::cout << std::setw(10);
    std::cout << GLYX << "  [ g/kg fuel ] " << std::endl;
    std::cout << std::setw(14);
    std::cout << " ---> MGLY: ";
    std::cout << std::setw(10);
    std::cout << MGLY << "  [ g/kg fuel ] " << std::endl;
    std::cout << std::setw(12);
    std::cout << " -> Soot: ";
    std::cout << std::setw(12);
    std::cout << Soot << "   [ g/kg fuel ] " << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

} /* End of Emission::Debug */

/* End of Emission.cpp */
