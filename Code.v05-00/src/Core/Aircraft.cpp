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

Aircraft::Aircraft( const char *aircraftName, double temperature_K, double pressure_Pa, double relHumidity_w )
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
        engine_ = Engine( engineName, temperature_K, pressure_Pa, \
                          relHumidity_w, machNumber_ );
        engNumber_ = 4;

        /* Dimensions */
        wingspan_ = 69.8; /* [m] */

        /* Weight */
        MTOW_ = 404.810E+03; /* [kg] */
        currMass_ = 0.9 * MTOW_; /* [kg] */

    }

    vortex_ = Vortex( temperature_K, pressure_Pa, 1.3E-02, wingspan_, \
                      currMass_, vFlight_ms_ );

} /* End of Aircraft::Aircraft */

Aircraft::Aircraft( const Aircraft &ac )
{

    Name_       = ac.Name_;
    vFlight_ms_ = ac.vFlight_ms_;
    machNumber_ = ac.machNumber_;
    wingspan_   = ac.wingspan_;
    MTOW_       = ac.MTOW_;
    currMass_   = ac.currMass_;
    engine_     = ac.engine_;
    engNumber_  = ac.engNumber_;
    vortex_     = ac.vortex_;

} /* End of Aircraft::Aircraft */

Aircraft& Aircraft::operator=( const Aircraft &ac )
{

    if ( &ac == this )
        return *this;

    Name_       = ac.Name_;
    vFlight_ms_ = ac.vFlight_ms_;
    machNumber_ = ac.machNumber_;
    wingspan_   = ac.wingspan_;
    MTOW_       = ac.MTOW_;
    currMass_   = ac.currMass_;
    engine_     = ac.engine_;
    engNumber_  = ac.engNumber_;
    vortex_     = ac.vortex_;
    return *this;

} /* End of Aircraft::operator= */

Aircraft::~Aircraft( )
{
    /* Destructor */

} /* End of Aircraft::~Aircraft */

double Aircraft::VortexLosses( const double EI_Soot, const double EI_SootRad, \
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
    const bool DEBUG = 0;

    /* Compute volume and mass of soot particles emitted */
    const double volParticle  = 4.0 / 3.0 * physConst::PI * pow( EI_SootRad, 3.0 );
    const double massParticle = volParticle * physConst::RHO_SOOT * 1.0E+03;

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
    const double beta_0  = +0.40;
    const double beta_1  = +1.19;
    const double alpha_0 = -1.35;

    /* Scaling for particle number */
    const double EIice_ref = 2.8E+14; /* [#/kg_fuel] */

    /* Number of particles emitted */
    const double EIice     = EI_Soot / massParticle; /* [#/kg_fuel] */

    /* z_Atm = Depth of the supersaturated layer */
    z_Atm  = wetDepth;

    /* z_Desc = Vertical displacement of the wake vortex */
    z_Desc = pow( 8.0 * vortex_.gamma() / ( physConst::PI * vortex_.N_BV() ), 0.5 );

    /* z_Emit = ... 
     * z_Emit is not implemented yet as z_Emit << z_Desc */
    z_Emit = 0.0E+00;

    /* ==================================== */
    /* ... Parameterization starts here ... */
    /* ==================================== */

    const double alpha_Atm  = +1.70 * pow( EIice / EIice_ref, -gamma_Atm  );
    const double alpha_Emit = +1.15 * pow( EIice / EIice_ref, -gamma_Emit );\
    const double alpha_Desc = -0.60;

    /* Combine each length scale into a single variable, zDelta, expressed in m. */
    z_Delta = + alpha_Atm  * z_Atm  \
              + alpha_Emit * z_Emit \
              + alpha_Desc * z_Desc;

    iceNumFrac = beta_0 + beta_1 / physConst::PI * atan( alpha_0 + z_Delta / 1.0E+02 );


    if ( DEBUG ) {
        std::cout << "----------------- DEBUG -----------------" << std::endl;
        std::cout << "In Aircraft::VortexLosses: " << std::endl;
        std::cout << "alpha_Atm  = " << alpha_Atm  << std::endl;
        std::cout << "alpha_Emit = " << alpha_Emit << std::endl;
        std::cout << "alpha_Desc = " << alpha_Desc << std::endl;
        std::cout << "z_Atm      = " << z_Atm   << " [m]" << std::endl;
        std::cout << "z_Emit     = " << z_Emit  << " [m]" << std::endl;
        std::cout << "z_Desc     = " << z_Desc  << " [m]" << std::endl;
        std::cout << "z_Delta    = " << z_Delta << " [m]" << std::endl;
        std::cout << "rem. frac. = " << iceNumFrac << " [-]" << std::endl;
    }

    if ( ( iceNumFrac < 0.1E+00 ) || ( iceNumFrac > 1.0E+00 ) ) {
        /* If we lose more than 90% of the initial number of ice
         * crystals, then print warning */
        std::cout << "In Aircraft::VortexLosses: ";
        std::cout << "iceNumFrac = " << iceNumFrac << std::endl;
        iceNumFrac = std::min( std::max( iceNumFrac, 0.0E+00 ), 1.0E+00 );
    }

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
    std::cout << "  -Aircraft Name          : " << Name_ << "\n";
    std::cout << "  -Aircraft Engine Number : " << engNumber_ << "\n";
    std::cout << "  -Aircraft Flight Speed  : " << vFlight_ms_ << "         [ m/s ]" << "\n";
    std::cout << "  -Aircraft Mach Number   : " << machNumber_ << "     [ - ]" << "\n";
    std::cout << "  -Aircraft Wingspan      : " << wingspan_ << "        [ m ]" << "\n";
    std::cout << "  -Aircraft MTOW          : " << MTOW_ << "  [ kg ]" << "\n";
    std::cout << "  -Aircraft Current Mass  : " << currMass_ << "  [ kg ]" << "\n";
    std::cout << " --- \n";
    std::cout << "      +Engine Name        : " << engine_.getName() << "\n";
    std::cout << "      +Engine EI NOx      : " << engine_.getEI_NOx() << "      [ g/kg ]" << "\n";
    std::cout << "      +Engine EI NO       : " << engine_.getEI_NO() << "      [ g/kg ]" << "\n";
    std::cout << "      +Engine EI NO2      : " << engine_.getEI_NO2() << "      [ g/kg ]" << "\n";
    std::cout << "      +Engine EI HNO2     : " << engine_.getEI_HNO2() << "     [ g/kg ]" << "\n";
    std::cout << "      +Engine EI CO       : " << engine_.getEI_CO() << "      [ g/kg ]" << "\n";
    std::cout << "      +Engine EI HC       : " << engine_.getEI_HC() << "   [ g/kg ]" << "\n";
    std::cout << "      +Engine EI Soot     : " << engine_.getEI_Soot() << "        [ g/kg ]" << "\n";
    std::cout << "      +Engine SootRad     : " << engine_.getSootRad() << "       [ m ]" << "\n";
    std::cout << "      +Engine FuelFlow    : " << engine_.getFuelFlow() << "         [ kg/s ]" << "\n";
    std::cout << "      +Engine theta       : " << engine_.getTheta() << "     [ - ]" << "\n";
    std::cout << "      +Engine delta       : " << engine_.getDelta() << "     [ - ]" << "\n";
    std::cout << " --- \n";
    std::cout << "      +Wake vortex separation : " << vortex_.b() << "   [ m ]" << "\n";
    std::cout << "      +Initial circulation    : " << vortex_.gamma() << "   [ m^2/s ]" << "\n";
    std::cout << "      +Wake eff. time scale   : " << vortex_.t() << "    [ s ]" << "\n";
    std::cout << "      +Initial velocity scale : " << vortex_.w() << "   [ m/s ]" << "\n";
    std::cout << "      +Normalized dissip. rate: " << vortex_.eps_star() << " [ - ]" << "\n";
    std::cout << "      +Max. downwash displace.: " << vortex_.delta_zw() << "   [ m ]" << "\n";
    std::cout << "      +Mean downwash displace.: " << vortex_.delta_z1() << "   [ m ]" << "\n";
    std::cout << "      +Initial contrail depth : " << vortex_.D1() << "   [ m ]" << "\n";
    std::cout << "\n";
    std::cout << "\n";

    std::cout.precision(ss);

} /* End of Aircraft::Debug */

/* End of Aircraft.cpp */
