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

#include <cmath>
#include <iostream>
#include "Util/PhysConstant.hpp"
#include "Core/Vortex.hpp"

Vortex::Vortex( double temperature_K, double pressure_Pa,  \
                double N_BV, double span, double mass, \
                double vFlight ) 
{

    /* This function uses a parametric model to estimate the initial depth, 
     * width and maximum and mean downward displacements, delta_zw_ and
     * delta_z1_.
     *
     * The parametrization is taken from:
     * Schumann, U. "A contrail cirrus prediction model." Geoscientific Model Development 5.3 (2012): 543-580.
     *
     * INPUT: 
     * - temperature_K: ambient temperature expressed in K
     * - pressure_Pa  : ambient pressure expressed in Pa
     * - N_BV         : Brunt-Väisala frequency expressed in s^-1
     * - span         : aircraft span in m
     * - mass         : aircraft mass in kg
     * - vFlight      : aircraft velocity in m/s
     *
     * The function computes:
     * - b_       : wake vortex separation in m
     * - gamma_   : initial circulation in m^2/s
     * - t_       : effective time scale in s
     * - w_       : initial velocity scale in m/s
     * - eps_star_: normalized dissipation rate
     * - delta_zw_: maximum sinking in m
     * - delta_z1_: initial sinking in m
     * - D1       : initial contrail depth in m
     * */

    /* Constructor */

    double rho = pressure_Pa / ( temperature_K * physConst::R_Air );

    /* Wake vortex separation, [ m ] */
    b_ = physConst::PI * span / 4;

    /* Initial circulation, [ m ^ 2 / s ] */
    gamma_ = 4 * mass * physConst::g / ( physConst::PI * span * rho * vFlight );
    
    /* Effective time scale, [ s ] */
    t_ = 2 * physConst::PI * b_ * b_ / gamma_;
    
    /* Initial velocity scale, [ m / s ] */
    w_ = gamma_ / ( 2 * physConst::PI * b_ );
    
    /* Normalized dissipation rate, [ - ] */
    eps_star_ = pow( physConst::EPSILON * b_, 1.0 / 3.0 ) / w_;
    if ( N_BV <= 0 ) {
        std::cout << "In Vortex::Vortex: Brunt-Vaisala frequency takes negative value, N_BV = " << N_BV << " [s^-1]\n";
        N_BV = 1.3E-02;
    }

    /* Allocate input to class */
    N_BV_ = N_BV;

    if ( N_BV_ * t_ >= N_BVt_threshold ) {

        /* Strongly stratified conditions */
        delta_zw_ = 1.49 * w_ / N_BV_;

    } else if ( eps_star_ < eps_threshold ) {

        /* Weakly stratified conditions */
        delta_zw_ = b_ * ( 1.88 + 7.68 * ( 1.0 - 4.07 * eps_star_ + 5.67 * eps_star_ * eps_star_ ) * ( 0.79 - N_BV_ * t_ ));

    } else {

        std::cout << "In Vortex::Vortex:: Neither N_BV * t >= " << N_BVt_threshold << " nor eps* < " << eps_threshold << " are valid\n";
        std::cout << "Setting delta_zw to " << delta_zw_default << "\n";
        delta_zw_ = delta_zw_default;

    }

    delta_z1_ = Cz1 * delta_zw_;

    D_1_ = CD_0 * delta_zw_;

} /* End of Vortex::Vortex */


/* End of Vortex.cpp */
