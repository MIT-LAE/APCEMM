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

Emission::Emission( Engine engine, Fuel fuel )
{
    /* Constructor */

    Populate_withEngine( engine );

    Populate_withFuel( fuel );

    if ( EI_CO2 <= 0.0 ) {
        std::cout << "CO2 emissions are negative: EI(CO_2) = " << EI_CO2 << "g/kg" << std::endl;
        EI_CO2 = 0.0;
    }
    if ( EI_H2O <= 0.0 ) {
        std::cout << "H2O emissions are negative: EI(H_2O) = " << EI_H2O << "g/kg" << std::endl;
        EI_H2O = 0.0;
    }
    if ( EI_SO2 <= 0.0 ) {
        std::cout << "SO2 emissions are negative: EI(SO_2) = " << EI_SO2 << "g/kg" << std::endl;
        EI_SO2 = 0.0;
    }
    if ( EI_NOx <= 0.0 ) {
        std::cout << "NOx emissions are negative: EI(NO_x) = " << EI_NOx << "g/kg" << std::endl;
        EI_NOx = 0.0;
    }
    if ( EI_CO <= 0.0 ) {
        std::cout << "CO emissions are negative: EI(CO) = " << EI_CO << "g/kg" << std::endl;
        EI_CO = 0.0;
    }
    if ( EI_HC <= 0.0 ) {
        std::cout << "HC emissions are negative: EI(HC) = " << EI_HC << "g/kg" << std::endl;
        EI_HC = 0.0;
    }
    if ( EI_OH <= 0.0 ) {
        std::cout << "OH emissions are negative: EI(OH) = " << EI_OH << "g/kg" << std::endl;
        EI_OH = 0.0;
    }

    if ( EI_Soot <= 0.0 ) {
        std::cout << "Soot emissions are negative: EI(Soot) = " << EI_Soot << "g/kg" << std::endl;
        EI_Soot = 0.0;
    }
    if ( SootRad <= 0.0 ) {
        std::cout << "Soot emissions are negative: Rad_Soot = " << SootRad << std::endl;
        SootRad = 1.0E-08;
    }


}

Emission::~Emission( )
{
    /* Destructor */
}

void Emission::Populate_withEngine( Engine engine )
{
    /* Gaseous species */
    /*** Engine characteristics */
    EI_NOx = engine.EI_NOx;
    EI_CO = engine.EI_CO;
    EI_HC = engine.EI_HC;
    EI_OH = engine.EI_OH;

    /* BC */
    EI_Soot = engine.EI_Soot;
    SootRad = engine.SootRad;

}

void Emission::Populate_withFuel( Fuel fuel )
{
    /*** Fuel characteristics */
    EI_CO2 = ( 12.0107 * fuel.atomC ) / ( 12.0107 * fuel.atomC + 1.0079 * fuel.atomH + 10.0067 * fuel.atomN + 32.065 * fuel.atomS ) * ( 44.0095 / 12.0107);
    EI_H2O = (  1.0079 * fuel.atomH ) / ( 12.0107 * fuel.atomC + 1.0079 * fuel.atomH + 10.0067 * fuel.atomN + 32.065 * fuel.atomS ) * ( 18.0153 /  1.0079 );
    EI_SO2 = fuel.FSC / 2.0 / 1000.0;

}

/* End of Emission.cpp */

