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

#include "Core/Aircraft.hpp"

Aircraft::Aircraft( )
{

    /* Default Constructor */

} /* End of Aircraft::Aircraft */

Aircraft::Aircraft( const char *aircraftName, std::string engineFilePath, double aircraftMass, \
                    double temperature_K, double pressure_Pa,  \
                    double relHumidity_w, double nBV )
{
    /* Constructor */

    Name_ = aircraftName;

    if ( Name_.compare( "B747" ) >= 0 ) {
        /* Flight characteristics */
        vFlight_ms_ = 250.0;
        machNumber_ = vFlight_ms_ \
                    / sqrt( physConst::GAMMA_Air * physConst::R_Air * temperature_K );

        /* Engine characteristics */
        const char *engineName = "GEnx-2B67B";
        engine_ = Engine( engineName, engineFilePath, temperature_K, pressure_Pa, \
                          relHumidity_w, machNumber_ );
        engNumber_ = 4;

        /* Dimensions */
        wingspan_ = 68.4; /* [m] */

        /* Weight */
        if ( aircraftMass > 0.0E+00 ) {
            currMass_ = aircraftMass;
        } else {
            throw std::range_error("Cannot have negative aircraft mass!");
        }

    }
    vortex_ = Vortex( temperature_K, pressure_Pa, nBV, wingspan_, \
                      currMass_, vFlight_ms_ );

} /* End of Aircraft::Aircraft */

Aircraft::Aircraft( const Input& input, std::string engineFilePath, std::string engineName){
    /* Flight characteristics */
    //Sets both flight speed and mach number
    setVFlight(input.flightSpeed(), input.temperature_K());

    /* Engine characteristics */
    engine_ = Engine( engineName.c_str(), engineFilePath, input.temperature_K(), input.pressure_Pa(), \
                        input.relHumidity_w(), machNumber_ );

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
    vortex_ = Vortex( input.temperature_K(), input.pressure_Pa(), input.nBV(), wingspan_, \
                    currMass_, vFlight_ms_ );

}

double Aircraft::VortexLosses( const double EI_Soot,    \
                                   const double EI_SootRad, \
                                   const double wetDepth )
{

    /* This function computes the fraction of contrail ice particles lost in
     * the vortex sinking phase because of turbulent sublimation.
     * This fraction is computed according to a parameterization derived from
     * large-eddy simulations and described in: 
     * Unterstrasser, Simon. "Properties of young contrails–a parametrisation based on large-eddy simulations."
     * Atmospheric Chemistry and Physics 16.4 (2016): 2059-2082.
     * */

    /* Debug flag? */
    // const bool ACDEBUG = 0;

    /* Compute volume and mass of soot particles emitted */
    const double volParticle  = 4.0 / 3.0 * physConst::PI * pow( EI_SootRad, 3.0 ); //EI_SootRad in m -> volume in m3
    const double massParticle = volParticle * physConst::RHO_SOOT * 1.0E+03; //Gives mass of a particle in grams

    /* Declare and initialize remaining fraction of ice crystals */
    double iceNumFrac = 0.0E+00;

    /* Declare and initialize each length scales */
    double z_Atm   = 0.0E+00;
    double z_Emit  = 0.0E+00;
    double z_Desc  = 0.0E+00;
    double z_Delta = 0.0E+00;

    /* Exponents */
    const double gamma_Atm  = 0.18;
    const double gamma_Emit = 0.18;

    /* Fitting coefficients */
    const double beta_0  = +0.45;
    const double beta_1  = +1.19;
    const double alpha_0 = -1.35;
    /* After reaching out to S. Unterstrasser, the beta_0 coefficient should
     * be 0.45 instead of 0.40. The paper contains a typo and equation (10a)
     * should be:
     * \beta_0 = 0.45
     */

    /* Scaling for particle number */
    const double EIice_ref = 2.8E+14; /* [#/kg_fuel] */

    /* Number of particles emitted */
    const double EIice     = EI_Soot / massParticle; /* [#/kg_fuel] */
    
    /* z_Atm = Depth of the supersaturated layer */
    z_Atm  = wetDepth;

    /* z_Desc = Vertical displacement of the wake vortex */
    if ( vortex_.N_BV() < 1.0E-05 )
        z_Desc = pow( 8.0 * vortex_.gamma() / ( physConst::PI * 1.3E-02        ), 0.5 );
    else
        z_Desc = pow( 8.0 * vortex_.gamma() / ( physConst::PI * vortex_.N_BV() ), 0.5 ); //Equation 4

    /* z_Emit = ... */
    z_Emit = 90;

    /* ==================================== */
    /* ... Parameterization starts here ... */
    /* ==================================== */

    const double alpha_Atm  = +1.70 * pow( EIice_ref / EIice, gamma_Atm  );
    const double alpha_Emit = +1.15 * pow( EIice_ref / EIice, gamma_Emit );
    const double alpha_Desc = -0.60;

    /* Combine each length scale into a single variable, zDelta, expressed in m. */
    z_Delta = + alpha_Atm  * z_Atm  \
              + alpha_Emit * z_Emit \
              + alpha_Desc * z_Desc;

    iceNumFrac = beta_0 + beta_1 / physConst::PI * atan( alpha_0 + z_Delta / 1.0E+02 );

//     if ( ( ACDEBUG ) || ( iceNumFrac < 0.1E+00 ) || ( iceNumFrac > 1.0E+00 ) ) {
//         /* If DEBUG is on or if we lose more than 90% of the initial number of 
//          * ice crystals, then print diagnostic */
// #pragma omp critical
//         {
//         if ( ACDEBUG )
//             std::cout << "----------------- DEBUG -----------------" << std::endl;
//         std::cout << "In Aircraft::VortexLosses: " << std::endl;
//         std::cout << "alpha_Atm  = " << alpha_Atm  << std::endl;
//         std::cout << "alpha_Emit = " << alpha_Emit << std::endl;
//         std::cout << "EIice      = " << EIice      << " [#/kg_fuel]" << std::endl;
//         std::cout << "EIice*     = " << EIice / EIice_ref << " [#/kg_fuel]" << std::endl;
//         std::cout << "alpha_Desc = " << alpha_Desc << std::endl;
//         std::cout << "z_Atm      = " << z_Atm  / 1.0E+02 << " [100m]" << std::endl;
//         std::cout << "z_Emit     = " << z_Emit / 1.0E+02 << " [100m]" << std::endl;
//         std::cout << "z_Desc     = " << z_Desc / 1.0E+02 << " [100m]" << std::endl;
//         std::cout << " -> Gamma  = " << vortex_.gamma() << " [m^2/s]" << std::endl;
//         std::cout << " -> N_BV   = " << vortex_.N_BV() << " [Hz]" << std::endl;
//         std::cout << "z_Delta    = " << z_Delta / 1.0E+02 << " [100m]" << std::endl;
//         std::cout << "rem. frac. = " << iceNumFrac << " [-]" << std::endl;
//         }
//     }

    iceNumFrac = std::min( std::max( iceNumFrac, 0.0E+00 ), 1.0E+00 );


    if ( iceNumFrac == 0.0E+00 )
        std::cout << "Contrail has fully melted because of vortex-sinking losses" << std::endl;

    return iceNumFrac;

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
    std::cout << std::setw(10) << vortex_.delta_zw() << " [ m ]" << "\n";
    std::cout << "      +Mean downwash displace. : ";
    std::cout << std::setw(10) << vortex_.delta_z1() << " [ m ]" << "\n";
    std::cout << "      +Initial contrail depth  : ";
    std::cout << std::setw(10) << vortex_.D1() << " [ m ]" << "\n";
    std::cout << "\n";
    std::cout << "\n";

    std::cout.precision(ss);

} /* End of Aircraft::Debug */

/* End of Aircraft.cpp */
