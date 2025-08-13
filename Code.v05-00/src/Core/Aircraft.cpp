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
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <cmath>
#include <iomanip>
#include <iostream>
#include "Util/PhysConstant.hpp"
#include "Core/Aircraft.hpp"

Aircraft::Aircraft( )
{

    /* Default Constructor */

} /* End of Aircraft::Aircraft */

Aircraft::Aircraft( const Input& input, std::string engineFilePath, std::string engineName){
    /* Flight characteristics */
    //Sets both flight speed and mach number

    T_CA_K_ = input.temperature_K(); // From the meteorology, at the reference altitude
    RHW_CA_PC_ = input.relHumidity_w(); // From the meteorology, at the reference altitude
    RHi_CA_PC_ = input.relHumidity_i(); // From the meteorology, at the reference altitude
    nBV_Hz_ = input.nBV(); // From the yaml input
    p_CA_Pa_ = input.pressure_Pa(); // From the yaml input

    setVFlight(input.flightSpeed(), T_CA_K_);

    /* Engine characteristics */
    engine_ = Engine( engineName.c_str(), engineFilePath, T_CA_K_, p_CA_Pa_, \
                        RHW_CA_PC_, machNumber_ );

    engNumber_ = input.numEngines();
    setFuelFlow( input.fuelFlow() );
    setEI_NOx( input.EI_NOx() );
    setEI_CO( input.EI_CO() );
    setEI_HC( input.EI_HC() );
    setEI_Soot( input.EI_Soot() );
    setSootRad( input.sootRad() );
    /* Dimensions */
    wingspan_ = input.wingspan();
    currMass_ = input.aircraftMass();
}

double Aircraft::VortexLosses( const double EI_Soot,    \
                                   const double EI_SootRad, \
                                   const double wetDepth )
{
    /* Compute volume and mass of soot particles emitted */
    const double volParticle  = 4.0 / 3.0 * physConst::PI * pow( EI_SootRad, 3.0 ); //EI_SootRad in m -> volume in m3
    const double massParticle = volParticle * physConst::RHO_SOOT * 1.0E+03; //Gives mass of a particle in grams

    /* Scaling for particle number */
    const double EIice_ref = 2.8E+14; /* [#/kg_fuel] */

    /* Number of particles emitted */
    const double EIice     = EI_Soot / massParticle; /* [#/kg_fuel] */

    vortex_ = Vortex( RHi_CA_PC_, T_CA_K, p_CA_Pa, nBV_Hz, wingspan_, \
                    currMass_, vFlight_ms_, WV_exhaust, N0 );

    return icenum_survfrac;

} /* End of Aircraft::VortexLosses */

void Aircraft::setEI_NOx(const double NOx)
{

    if ( NOx > 0.0E+00 )
        engine_.setEI_NOx(NOx);

} /* End of Aircraft::setEI_NOx */

void Aircraft::setEI_CO(const double CO)
{

    if ( CO > 0.0E+00 )
        engine_.setEI_CO(CO);

} /* End of Aircraft::setEI_CO */

void Aircraft::setEI_HC(const double HC)
{

    if ( HC > 0.0E+00 )
        engine_.setEI_HC(HC);

} /* End of Aircraft::setEI_HC */

void Aircraft::setEI_Soot(const double Soot)
{

    if ( Soot > 0.0E+00 )
        engine_.setEI_Soot(Soot);

} /* End of Aircraft::setEI_Soot */

void Aircraft::setSootRad(const double sootRad)
{

    if ( sootRad > 0.0E+00 )
        engine_.setSootRad(sootRad);

} /* End of Aircraft::setSootRad */

void Aircraft::setFuelFlow(const double ff)
{

    if ( ff > 0.0E+00 )
        engine_.setFuelFlow(ff / engNumber_);


} /* End of Aircraft::setFuelFlow */

void Aircraft::setVFlight(const double Vf, double temperature_K)
{

    if ( Vf > 0.0E+00 ){
        vFlight_ms_ = Vf;
        machNumber_ = vFlight_ms_ / sqrt( physConst::GAMMA_Air * physConst::R_Air * temperature_K );
    }


} /* End of Aircraft::setVFlight */

void Aircraft::setEngNumber(const double nEng)
{

    if ( nEng > 0.0E+00 )
        engNumber_ = nEng;


} /* End of Aircraft::setEngNumber */

void Aircraft::setWingspan(const double span)
{

    if ( span > 0.0E+00 )
        wingspan_ = span;


} /* End of Aircraft::setWingspan */

void Aircraft::setMass(double mass)
{
    currMass_ = mass;
}
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
    std::cout << "  -Aircraft Name          : ";
    std::cout << std::setw(10) << Name_ << "\n";
    std::cout << "  -Aircraft Engine Number : ";
    std::cout << std::setw(10) << engNumber_ << "\n";
    std::cout << "  -Aircraft Flight Speed  : ";
    std::cout << std::setw(10) << vFlight_ms_ << " [ m/s ]" << "\n";
    std::cout << "  -Aircraft Mach Number   : ";
    std::cout << std::setw(10) << machNumber_ << " [ - ]" << "\n";
    std::cout << "  -Aircraft Wingspan      : ";
    std::cout << std::setw(10) << wingspan_ << " [ m ]" << "\n";
    std::cout << "  -Aircraft Current Mass  : ";
    std::cout << std::setw(10) << currMass_ << " [ kg ]" << "\n";
    std::cout << " --- \n";

    std::cout << "      +Engine Name        : ";
    std::cout << std::setw(10) << engine_.getName() << "\n";
    std::cout << "      +Engine EI NOx      : ";
    std::cout << std::setw(10) << engine_.getEI_NOx() << " [ g/kg ]" << "\n";
    std::cout << "      +Engine EI NO       : ";
    std::cout << std::setw(10) << engine_.getEI_NO() << " [ g/kg ]" << "\n";
    std::cout << "      +Engine EI NO2      : ";
    std::cout << std::setw(10) << engine_.getEI_NO2() << " [ g/kg ]" << "\n";
    std::cout << "      +Engine EI HNO2     : ";
    std::cout << std::setw(10) << engine_.getEI_HNO2() << " [ g/kg ]" << "\n";
    std::cout << "      +Engine EI CO       : ";
    std::cout << std::setw(10) << engine_.getEI_CO() << " [ g/kg ]" << "\n";
    std::cout << "      +Engine EI HC       : ";
    std::cout << std::setw(10) << engine_.getEI_HC() << " [ g/kg ]" << "\n";
    std::cout << "      +Engine EI Soot     : ";
    std::cout << std::setw(10) << engine_.getEI_Soot() << " [ g/kg ]" << "\n";
    std::cout << "      +Engine SootRad     : ";
    std::cout << std::setw(10) << engine_.getSootRad() << " [ m ]" << "\n";
    std::cout << "      +Engine FuelFlow    : ";
    std::cout << std::setw(10) << engine_.getFuelFlow() << " [ kg/s ]" << "\n";
    std::cout << "      +Engine theta       : ";
    std::cout << std::setw(10) << engine_.getTheta() << " [ - ]" << "\n";
    std::cout << "      +Engine delta       : ";
    std::cout << std::setw(10) << engine_.getDelta() << " [ - ]" << "\n";
    std::cout << " --- \n";

    std::cout << "      +Brunt-Vaisala frequency : ";
    std::cout << std::setw(10) << vortex_.N_BV() << " [ Hz ]" << "\n";
    std::cout << "      +Wake vortex separation  : ";
    std::cout << std::setw(10) << vortex_.b() << " [ m ]" << "\n";
    std::cout << "      +Initial circulation     : ";
    std::cout << std::setw(10) << vortex_.gamma() << " [ m^2/s ]" << "\n";
    std::cout << "      +Wake eff. time scale    : ";
    std::cout << std::setw(10) << vortex_.t() << " [ s ]" << "\n";
    std::cout << "      +Initial velocity scale  : ";
    std::cout << std::setw(10) << vortex_.w() << " [ m/s ]" << "\n";
    std::cout << "      +Normalized dissip. rate : ";
    std::cout << std::setw(10) << vortex_.eps_star() << " [ - ]" << "\n";
    std::cout << "      +Max. downwash displace. : ";
    std::cout << std::setw(10) << vortex_.z_desc() << " [ m ]" << "\n";
    std::cout << "      +Mean downwash displace. : ";
    std::cout << std::setw(10) << vortex_.z_center() << " [ m ]" << "\n";
    std::cout << "      +Initial contrail depth  : ";
    std::cout << std::setw(10) << vortex_.D1() << " [ m ]" << "\n";
    std::cout << "\n";
    std::cout << "\n";

    std::cout.precision(ss);

} /* End of Aircraft::Debug */

/* End of Aircraft.cpp */
