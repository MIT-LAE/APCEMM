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

Vortex::Vortex( double RHi_PC, double temperature_K, double pressure_Pa,  \
                double N_BV, double wingspan, double ac_mass, \
                double vFlight, double WV_exhaust, double N0, double N0_ref);
{

    /* This function uses a parametric model to estimate the initial depth, 
     * width and maximum and mean downward displacements, z_desc_ and
     * z_center_.
     *
     * The parametrization is taken from:
     * Schumann, U. "A contrail cirrus prediction model." Geoscientific Model Development 5.3 (2012): 543-580.
     *
     * INPUT: 
     * - temperature_K: ambient temperature expressed in K
     * - pressure_Pa  : ambient pressure expressed in Pa
     * - N_BV         : Brunt-Väisala frequency expressed in s^-1
     * - wingspan         : aircraft wingspan in m
     * - ac_mass         : aircraft mass in kg
     * - vFlight      : aircraft velocity in m/s
     *
     * The function computes:
     * - b_       : wake vortex separation in m
     * - gamma_   : initial circulation in m^2/s
     * - t_       : effective time scale in s
     * - w_       : initial velocity scale in m/s
     * - eps_star_: normalized dissipation rate
     * - z_desc_  : maximum sinking in m of the contrail bottom
     * - z_center_: initial sinking in m of the contrail center
     * - D1       : initial contrail depth in m
     * */

    /* Constructor */

    double rho = pressure_Pa / ( temperature_K * physConst::R_Air );

    /* Wake vortex separation, [ m ] */
    b_ = physConst::PI * wingspan / 4;

    /* Initial circulation, [ m ^ 2 / s ] */
    gamma_ = 4 * ac_mass * physConst::g / ( physConst::PI * wingspan * rho * vFlight );
    
    /* Effective time scale, [ s ] */
    t_ = 2 * physConst::PI * b_ * b_ / gamma_;
    
    /* Initial velocity scale, [ m / s ] */
    w_ = gamma_ / ( 2 * physConst::PI * b_ );
    
    /* Normalized dissipation rate, [ - ] */
    if ( N_BV <= 1.0E-05 ) {
        std::cout << "In Vortex::Vortex: Brunt-Vaisala frequency is beneath 1.0E-5 (N_BV = " << N_BV << " [s^-1]). Taking the default value 1.3E-02.\n";
        N_BV = 1.3E-02;
    }
    N_BV_ = N_BV;

    /* Maximum downwash displacement, Eq. 5 in Lottermosser and Unterstrasser (2025) */
    z_desc_ = pow( 8.0 * gamma_ / ( physConst::PI * N_BV_ ), 0.5 );

    /* Height an air parcel has to descend until it is no longer supersaturated */
    s_ = RHi_PC / 100 - 1; // RHi in % -> excess supersaturation ratio, See S2 in U2016
    z_atm_  = 607.46 * pow(s_, 0.897) * pow(temperature_K / 205.0, 2.225); // Eq. A2 in LU2025

    /* Plume Radius */
    r_p_ = 1.5 + 0.314 * wingspan_; /* [m], from Eq. A6 in U2016 */
    
    /* Plume area before vortex breakup*/
    plume_area_0_ = 2 * physConst::PI * pow(r_p_, 2); /* [m2], see Appendix 2 in LU2025 */

    /* Temperature - 205 K*/
    T_205_ = temperature_K - 205.0; /* [K], from Eq. A3 in LU2025*/

    /* Height an air parcel has to descend until it just saturated when the emitted water vapor is added */
    rho_emit_ = WV_exhaust / plume_area_0_; // Eq. 6 in U2016
    const double rho_divisor = 10.; // 10 mg per m3
    z_emit_ = 1106.6 * pow(rho_emit * 1000. / rho_divisor, 0.678 + 0.0116 * T_205) * exp((-(0.0807+0.000428*T_205)*T_205));
    
    // /* Combine each length scale into a single variable, zDelta, expressed in m. */
    // const double N0_ref = 3.38E12; /* [#/m], hardcoded but valid for an A350 */
    n0_star_ = N0_ / N0_ref; /* [-], See Appendix A1 in LU 2025, plume area cancelled out */ 
    const double Psi = 1 / n0_star_; // Eq. 10 in LU2025
    z_delta_ = pow(Psi, gamma_exp_) * \
             (+ alpha_atm_  * z_atm_  \
              + alpha_emit_ * z_emit_) \
              - alpha_desc_ * z_desc_; // Eq. 9 in LU2025

    /* Compute the remaining fraction of ice crystals, from Eq. 12 in LU2025 */
    iceNumFrac_ = beta_0_ + beta_1_ / physConst::PI * atan( alpha_0_ + z_delta_ / 1.0E+02 );
    iceNumFrac_ = std::min( std::max( iceNumFrac_, 0.0E+00 ), 1.0E+00 );
    z_center_ = Cz1 * z_desc_;

    D_1_ = CD_0 * z_desc_;

} /* End of Vortex::Vortex */


/* End of Vortex.cpp */
