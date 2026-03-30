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
                double vFlight, double WV_exhaust, double N_postjet, \
                double N0_ref, double wingspan_ref)
{
    /* This constructor implements a parametric model for the vortex phase of contrail evolution,
     * estimating the initial depth, width, and the maximum and mean downward displacements (z_desc_ and z_center_).
     *
     * The formulae are based on Unterstrasser (2016), referred to as U2016; and Lottermoser & Unterstrasser (2025), referred to as LU2025.
     *
     * U2016 is available at https://doi.org/10.5194/acp-16-2059-2016
     * LU2025 is available at https://doi.org/10.5194/acp-25-7903-2025
     * 
     * INPUT:
     * - RHi_PC         : Relative humidity with respect to ice [%]
     * - temperature_K  : Ambient temperature [K]
     * - pressure_Pa    : Ambient pressure [Pa]
     * - N_BV           : Brunt-Väisälä frequency [s^-1]
     * - wingspan       : Aircraft wingspan [m]
     * - ac_mass        : Aircraft mass [kg]
     * - vFlight        : Aircraft velocity [m/s]
     * - WV_exhaust     : Water vapor exhaust [kg]
     * - N_postjet      : Ice crystal number after the jet regime [#]
     * - N0_ref         : Reference ice crystal number [#]
     * - wingspan_ref   : Reference wingspan [m]
     *
     * OUTPUT (computed member variables):
     * - gamma_              : Initial circulation [m^2/s]
     * - z_desc_             : Maximum sinking of contrail bottom [m]
     * - z_center_           : Initial sinking of contrail center [m]
     * - r_p_                : Plume radius [m]
     * - plume_area_0_       : Plume area before vortex breakup [m^2]
     * - T_205_              : Temperature offset from 205 K [K]
     * - z_atm_              : Height for air parcel to descend to lose supersaturation [m]
     * - z_emit_             : Height for air parcel to descend to just saturate with emitted WV [m]
     * - n0_star_            : Normalized initial ice crystal number [-]
     * - z_delta_fns_        : Combined length scale for survival fraction [m]
     * - icenum_survfrac_    : Fraction of surviving ice crystals [-]
     * - z_delta_h_          : Combined length scale for contrail height [m]
     * - icenum_survfrac_h_  : Fraction of surviving ice crystals for height [-]
     * - bhat_               : Height scaling factor [-]
     * - depth_mature_      : Parametrized contrail height [m]
     * - width_rect_mature_  : Rectangle-equivalent width of initial mature plume [m]
     */

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
    z_desc_ = sqrt( 8.0 * gamma_ / ( physConst::PI * N_BV_ ) );

    /* Height an air parcel has to descend until it is no longer supersaturated */
    s_ = RHi_PC / 100 - 1; // RHi in % -> excess supersaturation ratio, See S2 in U2016
    z_atm_  = 607.46 * pow(s_, 0.897) * pow(temperature_K / 205.0, 2.225); // Eq. A2 in LU2025

    /* Plume Radius */
    r_p_ = 1.5 + 0.314 * wingspan; /* [m], from Eq. A6 in U2016 */
    r_p_ref_ = 1.5 + 0.314 * wingspan_ref; /* [m], from Eq. A6 in U2016 */

    /* Plume area before vortex breakup*/
    plume_area_0_ = 2 * physConst::PI * r_p_ * r_p_; /* [m2], see Appendix 2 in LU2025 */
    plume_area_0_ref_ = 2 * physConst::PI * r_p_ref_ * r_p_ref_; /* [m2], see Appendix 2 in LU2025 */

    /* Temperature - 205 K*/
    T_205_ = temperature_K - 205.0; /* [K], from Eq. A3 in LU2025*/

    /* Height an air parcel has to descend until it just saturated when the emitted water vapor is added */
    rho_emit_ = WV_exhaust / plume_area_0_; // Eq. 6 in U2016
    std::cout << "Density of emitted water vapor: " << rho_emit_ << " [kg/m^3]" << std::endl;
    const double rho_divisor = 10.; // 10 mg per m3
    z_emit_ = 1106.6 * pow(rho_emit_ * 1000. / rho_divisor, 0.678 + 0.0116 * T_205_) * exp((-(0.0807+0.000428*T_205_)*T_205_));

    /* Combine each length scale into a single variable, zDelta, expressed in m. */
    n0_ = N_postjet / plume_area_0_;
    n0_ref_ = N0_ref / plume_area_0_ref_;
    n0_star_ = n0_ / n0_ref_; /* [-], See Appendix A1 in LU 2025 */
    const double Psi = 1 / n0_star_; // Eq. 10 in LU2025
    z_delta_fns_ = pow(Psi, gamma_exp_) * \
             (+ alpha_atm_  * z_atm_  \
              + alpha_emit_ * z_emit_) \
              - alpha_desc_ * z_desc_; // Eq. 9 in LU2025

    /* Compute the remaining fraction of ice crystals, from Eq. 12 in LU2025 */
    icenum_survfrac_ = beta_0_ + beta_1_ / physConst::PI * atan( alpha_0_ + z_delta_fns_ / 1.0E+02 );
    icenum_survfrac_ = std::min( std::max( icenum_survfrac_, 0.0E+00 ), 1.0E+00 );

    if ( icenum_survfrac_ == 0.0E+00 )
        std::cout << "Contrail has fully melted because of vortex-sinking losses" << std::endl;

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
        bhat_ = eta_2_ * icenum_survfrac_h_ + (eta_1_ - eta_2_) * x_s_;
    }
    depth_mature_ = bhat_ * z_desc_; // Eq. 12 in U2016

    // Rectangle-equivalent width of the initial mature plume
    width_rect_mature_ = 0.63 * wingspan;

    // Contrail area
    area_mature_ = width_rect_mature_ * depth_mature_; // [m^2]

    // The height of the center is the maximum displacement minus half the contrail height
    z_center_ = z_desc_ - depth_mature_ / 2.0;

} /* End of Vortex::Vortex */


/* End of Vortex.cpp */
