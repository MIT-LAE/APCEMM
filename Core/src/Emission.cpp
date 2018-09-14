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
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Emission.hpp"

Emission::Emission( )
{
    /* Constructor */

} /* End of Emission::Emission */

void Emission::Populate( const Engine &engine, const Fuel &fuel )
{
    /* Constructor */

    Populate_withEngine( engine );

    Populate_withFuel( fuel );
    
    /* Carbon footprint */
    double carbContent = 0.95; 
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
    if ( SO2 <= 0.0 ) {
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

    if ( Soot <= 0.0 ) {
        std::cout << "Soot emissions are negative: EI(Soot) = " << Soot << "g/kg" << std::endl;
        Soot = 0.0;
    }
    if ( SootRad <= 1.0E-09 ) {
        std::cout << "Soot emissions are negative: Rad_Soot = " << SootRad << std::endl;
        SootRad = 1.0E-09;
    }

    engineName = engine.Name;

    fuelChem =  fuel.ChemFormula;

} /* End of Emission::Populate */

Emission::~Emission( )
{
    /* Destructor */

} /* End of Emission::~Emission */

void Emission::Populate_withEngine( const Engine &engine )
{
    /* Gaseous species */
    /*** Engine characteristics */
    NOx  = engine.EI_NOx;
    NO   = engine.EI_NO;
    NO2  = engine.EI_NO2;
    HNO2 = engine.EI_HNO2;
    CO   = engine.EI_CO;
    HC   = engine.EI_HC;

    /* BC */
    Soot = engine.EI_Soot;
    SootRad = engine.SootRad;

} /* End of Emission::Populate_withEngine */

void Emission::Populate_withFuel( const Fuel &fuel )
{
    /*** Fuel characteristics */
    CO2 = ( 12.0107 * fuel.atomC ) / ( 12.0107 * fuel.atomC + 1.0079 * fuel.atomH + 10.0067 * fuel.atomN + 32.065 * fuel.atomS ) * ( 44.0095 / ( 1 * 12.0107 ) ); /* [ g/g fuel ] */
    H2O = (  1.0079 * fuel.atomH ) / ( 12.0107 * fuel.atomC + 1.0079 * fuel.atomH + 10.0067 * fuel.atomN + 32.065 * fuel.atomS ) * ( 18.0153 / ( 2 *  1.0079 ) ); /* [ g/g fuel ] */
    CO2 *= 1000; /* [ g/kg fuel ]*/
    H2O *= 1000; /* [ g/kg fuel ]*/
    SO2 = fuel.FSC / 2.0 / 1000.0;

} /* End of Emission::Populate_withFuel */

Emission Emission::operator+( const Emission &emission_add )
{
    Emission emission;
    emission.CO2 = this->CO2 + emission_add.CO2 ;
    emission.H2O = this->H2O + emission_add.H2O ;
    emission.NOx = this->NOx + emission_add.NOx ;
    emission.NO  = this->NO  + emission_add.NO  ;
    emission.NO2 = this->NO2 + emission_add.NO2 ;
    emission.HNO2= this->HNO2+ emission_add.HNO2;
    emission.SO2 = this->SO2 + emission_add.SO2 ;
    emission.CO  = this->CO  + emission_add.CO  ;
    emission.HC  = this->HC  + emission_add.HC  ;
    emission.CH4 = this->CH4 + emission_add.CH4 ;
    emission.C2H6= this->C2H6+ emission_add.C2H6;
    emission.PRPE= this->PRPE+ emission_add.PRPE;
    emission.ALK4= this->ALK4+ emission_add.ALK4;
    emission.CH2O= this->CH2O+ emission_add.CH2O;
    emission.ALD2= this->ALD2+ emission_add.ALD2;
    emission.GLYX= this->GLYX+ emission_add.GLYX;
    emission.MGLY= this->MGLY+ emission_add.MGLY;
    emission.Soot= this->Soot+ emission_add.Soot;

    if ( this->SootRad != emission_add.SootRad )
        std::cout << " Emission:operator+: SootRad differs for both Emission instances" << std::endl;

    emission.SootRad= this->SootRad;
    return emission;

} /* End of Emission::operator+ */

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

