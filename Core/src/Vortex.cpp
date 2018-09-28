/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Vortex Program File                                              */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Vortex.cpp                                */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Vortex.hpp"

Vortex::Vortex( ) 
{

    /* Default Constructor */

} /* End of Vortex::Vortex */

Vortex::Vortex( double temperature_K, double pressure_Pa, double n_BV, double span, double mass, double vFlight ) 
{

    /* Constructor */

    double rho = pressure_Pa / ( temperature_K * physConst::R_Air );

    /* Wake vortex separation, [ m ] */
    b = physConst::PI * span / 4;
    /* Initial circulation, [ m ^ 2 / s ] */
    gamma = 4 * mass * physConst::g / ( physConst::PI * span * rho * vFlight );
    /* Effective time scale, [ s ] */
    t = 2 * physConst::PI * b * b / gamma;
    /* Initial velocity scale, [ m / s ] */
    w = gamma / ( 2 * physConst::PI * b );
    /* Normalized dissipation rate, [ - ] */
    eps_star = pow( physConst::EPSILON * b, double(1.0/3.0) ) / w;

    if ( n_BV <= 0 ) {
        std::cout << "In Vortex::Vortex: Brunt-Vaisala frequency takes negative value, n_BV = " << n_BV << " [s^-1]\n";
        n_BV = 1.3E-02;
    }

    if ( n_BV * t >= n_BVt_threshold ) {
        delta_zw = 1.49 * w / n_BV;
    } else if ( eps_star < eps_threshold ) {
        delta_zw = b * ( 1.88 + 7.68 * ( 1.0 - 4.07 * eps_star + 5.67 * eps_star * eps_star ) * ( 0.79 - n_BV * t ));
    } else {
        std::cout << "In Vortex::Vortex:: Neither n_BV * t >= " << n_BVt_threshold << " nor eps* < " << eps_threshold << " are valid\n";
        std::cout << "Setting delta_zw to " << delta_zw_default << "\n";
        delta_zw = delta_zw_default;
    }

    delta_z1 = Cz1 * delta_zw;

    D_1 = CD_0 * delta_zw;


} /* End of Vortex::Vortex */

Vortex::~Vortex( )
{

    /* Destructor */

} /* End of Vortex::~Vortex */

Vortex::Vortex( const Vortex &v )
{

    b = v.b;
    gamma = v.gamma;
    t = v.t;
    w = v.w;
    eps_star = v.eps_star;
    delta_zw = v.delta_zw;
    delta_z1 = v.delta_z1;
    D_1 = v.D_1;

} /* End of Vortex::Vortex */

Vortex& Vortex::operator=( const Vortex &v )
{

    if ( &v == this )
        return *this;

    b = v.b;
    gamma = v.gamma;
    t = v.t;
    w = v.w;
    eps_star = v.eps_star;
    delta_zw = v.delta_zw;
    delta_z1 = v.delta_z1;
    D_1 = v.D_1;
    return *this;

} /* End of Vortex::operator= */

double Vortex::getb() const
{

    return b;

} /* End of Vortex::getb */

double Vortex::getgamma() const
{

    return gamma;

} /* End of Vortex::getgamma */

double Vortex::gett() const
{

    return t;

} /* End of Vortex::gett */

double Vortex::getw() const
{

    return w;

} /* End of Vortex::getw */

double Vortex::geteps() const
{

    return eps_star;

} /* End of Vortex::geteps */

double Vortex::getdeltazw() const
{

    return delta_zw;

} /* End of Vortex::getdeltazw */

double Vortex::getdeltaz1() const
{

    return delta_z1;

} /* End of Vortex::getdeltaz1 */

double Vortex::getD1() const
{

    return D_1;

} /* End of Vortex::getD1 */

/* End of Vortex.cpp */
