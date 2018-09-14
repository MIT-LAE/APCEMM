/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Aircraft Program File                                            */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Aircraft.cpp                              */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Aircraft.hpp"

Aircraft::Aircraft( const char *aircraftName, double temperature_K, double pressure_Pa, double relHumidity_w )
{
    /* Constructor */

    Name = aircraftName;

    if ( Name.compare( "B747" ) >= 0 ) {
        /* Flight characteristics */
        vFlight_ms = 250.0;
        machNumber = vFlight_ms / sqrt( gamma * R_Air * temperature_K );

        /* Engine characteristics */
        const char *engineName = "GEnx-2B67B";
        engine.Build( engineName, temperature_K, pressure_Pa, relHumidity_w, machNumber );
        engNumber = 4;

        /* Dimensions */
        wingSpan = 69.8; /* [m] */

        /* Weight */
        MTOW = 404.810E+03; /* [kg] */
        currMass = 0.9 * MTOW; /* [kg] */

    }

} /* End of Aircraft::Aircraft */

Aircraft::~Aircraft( )
{
    /* Destructor */

} /* End of Aircraft::~Aircraft */

void Aircraft::Debug( ) const
{

    std::streamsize ss = std::cout.precision();

    std::cout.precision(5);
    
    std::cout << std::endl;
    std::cout << "**** Input Debugger ****" << std::endl;
    std::cout << "Aircraft and Engine parameters: " << std::endl;
    std::cout << std::endl;
    std::cout << std::setw(15);
    std::cout << "Variable";
    std::cout << std::setw(19);
    std::cout << "Value";
    std::cout << std::setw(12);
    std::cout << "Units" << std::endl;
    std::cout << "  -Aircraft Name          : " << Name << std::endl;
    std::cout << "  -Aircraft Engine Number : " << engNumber << std::endl;
    std::cout << "  -Aircraft Flight Speed  : " << vFlight_ms << "         [ m/s ]" << std::endl;
    std::cout << "  -Aircraft Mach Number   : " << machNumber << "     [ - ]" << std::endl;
    std::cout << "  -Aircraft WingSpan      : " << wingSpan << "        [ m ]" << std::endl;
    std::cout << "  -Aircraft MTOW          : " << MTOW << "  [ kg ]" << std::endl;
    std::cout << "  -Aircraft Current Mass  : " << currMass << "  [ kg ]" << std::endl;
    std::cout << "      +Engine Name        : " << engine.Name << std::endl;
    std::cout << "      +Engine EI NOx      : " << engine.EI_NOx << "      [ g/kg ]" << std::endl;
    std::cout << "      +Engine EI NO       : " << engine.EI_NO << "      [ g/kg ]" << std::endl;
    std::cout << "      +Engine EI NO2      : " << engine.EI_NO2 << "      [ g/kg ]" << std::endl;
    std::cout << "      +Engine EI HNO2     : " << engine.EI_HNO2 << "     [ g/kg ]" << std::endl;
    std::cout << "      +Engine EI CO       : " << engine.EI_CO << "      [ g/kg ]" << std::endl;
    std::cout << "      +Engine EI HC       : " << engine.EI_HC << "   [ g/kg ]" << std::endl;
    std::cout << "      +Engine EI Soot     : " << engine.EI_Soot << "        [ g/kg ]" << std::endl;
    std::cout << "      +Engine SootRad     : " << engine.SootRad << "       [ m ]" << std::endl;
    std::cout << "      +Engine FuelFlow    : " << engine.fuelflow << "         [ kg/s ]" << std::endl;
    std::cout << "      +Engine theta       : " << engine.theta << "     [ - ]" << std::endl;
    std::cout << "      +Engine delta       : " << engine.delta << "     [ - ]" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout.precision(ss);

} /* End of Aircraft::Debug */

/* End of Aircraft.cpp */
