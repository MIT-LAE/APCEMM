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

Aircraft::Aircraft( )
{

    /* Default Constructor */

} /* End of Aircraft::Aircraft */

Aircraft::Aircraft( const char *aircraftName, double temperature_K, double pressure_Pa, double relHumidity_w )
{
    /* Constructor */

    Name = aircraftName;

    if ( Name.compare( "B747" ) >= 0 ) {
        /* Flight characteristics */
        vFlight_ms = 250.0;
        machNumber = vFlight_ms / sqrt( GAMMA_Air * R_Air * temperature_K );

        /* Engine characteristics */
        const char *engineName = "GEnx-2B67B";
        engine = Engine( engineName, temperature_K, pressure_Pa, relHumidity_w, machNumber );
        engNumber = 4;

        /* Dimensions */
        wingSpan = 69.8; /* [m] */

        /* Weight */
        MTOW = 404.810E+03; /* [kg] */
        currMass = 0.9 * MTOW; /* [kg] */

    }

    vortex = Vortex( temperature_K, pressure_Pa, 1.3E-02, wingSpan, currMass, vFlight_ms  );

} /* End of Aircraft::Aircraft */

Aircraft::Aircraft( const Aircraft &ac )
{

    Name = ac.Name;
    vFlight_ms = ac.vFlight_ms;
    machNumber = ac.machNumber;
    wingSpan = ac.wingSpan;
    MTOW = ac.MTOW;
    currMass = ac.currMass;
    engine = ac.engine;
    engNumber = ac.engNumber;

} /* End of Aircraft::Aircraft */

Aircraft& Aircraft::operator=( const Aircraft &ac )
{

    if ( &ac == this )
        return *this;

    Name = ac.Name;
    vFlight_ms = ac.vFlight_ms;
    machNumber = ac.machNumber;
    wingSpan = ac.wingSpan;
    MTOW = ac.MTOW;
    currMass = ac.currMass;
    engine = ac.engine;
    engNumber = ac.engNumber;
    return *this;

} /* End of Aircraft::operator= */
Aircraft::~Aircraft( )
{
    /* Destructor */

} /* End of Aircraft::~Aircraft */

std::string Aircraft::getName() const
{

    return Name;

} /* End of Aircraft::getName */

double Aircraft::getVFlight() const
{

    return vFlight_ms;

} /* End of Aircraft::getVFlight */

double Aircraft::getMach() const
{

    return machNumber;

} /* End of Aircraft::getMach */

double Aircraft::getWingSpan() const
{

    return wingSpan;

} /* End of Aircraft::getWingSpan */

double Aircraft::getMTOW() const
{

    return MTOW;

} /* End of Aircraft::getMTOW */

double Aircraft::getCurrMass() const
{

    return currMass;

} /* End of Aircraft::currMass */

Engine Aircraft::getEngine() const
{

    return engine;

} /* End of Aircraft::getEngine */

double Aircraft::getFuelFlow() const
{

    return engine.getFuelFlow() * engNumber;

} /* End of Aircraft::getFuelFlow() */

unsigned int Aircraft::getEngNumber() const
{

    return engNumber;

} /* End of Aircraft::getEngNumber */

double Aircraft::getVortexdeltaz1() const
{

    return vortex.getdeltaz1();

} /* End of Aircraft::getVortexDeltaz1 */

double Aircraft::getVortexdeltazw() const
{

    return vortex.getdeltazw();

} /* End of Aircraft::getVortexDeltazw */

void Aircraft::Debug( ) const
{

    std::streamsize ss = std::cout.precision();

    std::cout.precision(5);
    
    std::cout << "\n";
    std::cout << "**** Input Debugger ****" << "\n";
    std::cout << "Aircraft and Engine parameters: " << "\n";
    std::cout << "\n";
    std::cout << std::setw(15);
    std::cout << "Variable";
    std::cout << std::setw(19);
    std::cout << "Value";
    std::cout << std::setw(12);
    std::cout << "Units" << "\n";
    std::cout << "  -Aircraft Name          : " << Name << "\n";
    std::cout << "  -Aircraft Engine Number : " << engNumber << "\n";
    std::cout << "  -Aircraft Flight Speed  : " << vFlight_ms << "         [ m/s ]" << "\n";
    std::cout << "  -Aircraft Mach Number   : " << machNumber << "     [ - ]" << "\n";
    std::cout << "  -Aircraft WingSpan      : " << wingSpan << "        [ m ]" << "\n";
    std::cout << "  -Aircraft MTOW          : " << MTOW << "  [ kg ]" << "\n";
    std::cout << "  -Aircraft Current Mass  : " << currMass << "  [ kg ]" << "\n";
    std::cout << " --- \n";
    std::cout << "      +Engine Name        : " << engine.getName() << "\n";
    std::cout << "      +Engine EI NOx      : " << engine.getEI_NOx() << "      [ g/kg ]" << "\n";
    std::cout << "      +Engine EI NO       : " << engine.getEI_NO() << "      [ g/kg ]" << "\n";
    std::cout << "      +Engine EI NO2      : " << engine.getEI_NO2() << "      [ g/kg ]" << "\n";
    std::cout << "      +Engine EI HNO2     : " << engine.getEI_HNO2() << "     [ g/kg ]" << "\n";
    std::cout << "      +Engine EI CO       : " << engine.getEI_CO() << "      [ g/kg ]" << "\n";
    std::cout << "      +Engine EI HC       : " << engine.getEI_HC() << "   [ g/kg ]" << "\n";
    std::cout << "      +Engine EI Soot     : " << engine.getEI_Soot() << "        [ g/kg ]" << "\n";
    std::cout << "      +Engine SootRad     : " << engine.getSootRad() << "       [ m ]" << "\n";
    std::cout << "      +Engine FuelFlow    : " << engine.getFuelFlow() << "         [ kg/s ]" << "\n";
    std::cout << "      +Engine theta       : " << engine.getTheta() << "     [ - ]" << "\n";
    std::cout << "      +Engine delta       : " << engine.getDelta() << "     [ - ]" << "\n";
    std::cout << " --- \n";
    std::cout << "      +Wake vortex separation : " << vortex.getb() << "   [ m ]" << "\n";
    std::cout << "      +Initial circulation    : " << vortex.getgamma() << "   [ m^2/s ]" << "\n";
    std::cout << "      +Wake eff. time scale   : " << vortex.gett() << "    [ s ]" << "\n";
    std::cout << "      +Initial velocity scale : " << vortex.getw() << "   [ m/s ]" << "\n";
    std::cout << "      +Normalized dissip. rate: " << vortex.geteps() << " [ - ]" << "\n";
    std::cout << "      +Max. downwash displace.: " << vortex.getdeltazw() << "   [ m ]" << "\n";
    std::cout << "      +Mean downwash displace.: " << vortex.getdeltaz1() << "   [ m ]" << "\n";
    std::cout << "      +Initial contrail depth : " << vortex.getD1() << "   [ m ]" << "\n";
    std::cout << "\n";
    std::cout << "\n";

    std::cout.precision(ss);

} /* End of Aircraft::Debug */

/* End of Aircraft.cpp */
