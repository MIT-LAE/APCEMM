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
     * - gamma_   : initial circulation in m^2/s
     * - z_desc_  : maximum sinking in m of the contrail bottom
     * - z_center_: initial sinking in m of the contrail center
     * - D1       : initial contrail depth in m
     * */

    /* Constructor */

    double rho = pressure_Pa / ( temperature_K * physConst::R_Air );

    /* Initial circulation, [ m ^ 2 / s ] */
    gamma_ = 4 * ac_mass * physConst::g / ( physConst::PI * wingspan * rho * vFlight );
    
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

    /* Combine each length scale into a single variable, zDelta, expressed in m. */
    const double N0_ref = 3.38E12; /* [#/m], hardcoded but valid for an A350 */
    n0_star_ = N0_ / N0_ref; /* [-], See Appendix A1 in LU 2025, plume area cancelled out */ 
    const double Psi = 1 / n0_star_; // Eq. 10 in LU2025
    z_delta_fns_ = pow(Psi, gamma_exp_) * \
             (+ alpha_atm_  * z_atm_  \
              + alpha_emit_ * z_emit_) \
              - alpha_desc_ * z_desc_; // Eq. 9 in LU2025

    /* Compute the remaining fraction of ice crystals, from Eq. 12 in LU2025 */
    icenum_survfrac_ = beta_0_ + beta_1_ / physConst::PI * atan( alpha_0_ + z_delta_fns_ / 1.0E+02 );
    icenum_survfrac_ = std::min( std::max( icenum_survfrac_, 0.0E+00 ), 1.0E+00 );

    /* Compute the parametrized contrail height by repeating the calculations for the survival fraction
     * with gamma set to 0 as in U2016 */
    z_delta_h_ =+ alpha_atm_  * z_atm_  \
                + alpha_emit_ * z_emit_ \
                - alpha_desc_ * z_desc_; // Eq. 9 in LU2025 with gamma set to 0 as in U2016
    
    icenum_survfrac_h_ = beta_0_ + beta_1_ / physConst::PI * atan( alpha_0_ + z_delta_h_ / 1.0E+02 );
    icenum_survfrac_h_ = std::min( std::max( icenum_survfrac_h_, 0.0E+00 ), 1.0E+00 );

    // Eq. 13 in U2016
    if ( icenum_survfrac_h_ <= x_s_ ) {
        bhat_ = eta_1_ * icenum_survfrac_h_;
    }
    else {
        bhat_ = eta_2_ * icenum_survfrac_h_ + (eta_1_ - eta_2_) * x_s;
    }
    height_mature_ = bhat_ * z_desc_; // Eq. 12 in U2016

    // Rectangle-equivalent width of the initial mature plume
    width_rect_mature_ = 0.63 * wingspan;

    // The height of the center is the maximum displacement minus half the contrail height
    z_center_ = z_desc_ - height_mature_ / 2.0;

} /* End of Vortex::Vortex */


/* End of Vortex.cpp */
