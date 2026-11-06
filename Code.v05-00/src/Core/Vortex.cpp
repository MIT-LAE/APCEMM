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

#include <algorithm>
#include <cmath>
#include <iostream>
#include "Util/PhysConstant.hpp"
#include "Core/Vortex.hpp"

Vortex::Vortex( double RHi_PC, double temperature_K, double pressure_Pa,  \
                double N_BV, double wingspan, double ac_mass, \
                double vFlight, double WV_exhaust, \
                double N0_ref, double wingspan_ref ) :
    RHi_PC_( RHi_PC ), temperature_K_( temperature_K ), wingspan_( wingspan ),
    vFlight_( vFlight ), WV_exhaust_( WV_exhaust ), N0_ref_( N0_ref ),
    wingspan_ref_( wingspan_ref )
{
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

    /* Height an air parcel has to descend until it is no longer supersaturated.
     * RHi in % -> excess supersaturation ratio, See S2 in U2016.
     * Subsaturated air has no supersaturation to lose, so z_atm is 0. */
    const double s = std::max( RHi_PC / 100 - 1, 0.0 );
    z_atm_  = 607.46 * pow(s, 0.897) * pow(temperature_K / 205.0, 2.225); // Eq. A2 in LU2025

    /* Plume Radius */
    const double r_p = 1.5 + 0.314 * wingspan; /* [m], from Eq. A6 in U2016 */
    const double r_p_ref = 1.5 + 0.314 * wingspan_ref; /* [m], from Eq. A6 in U2016 */

    /* Plume area before vortex breakup*/
    plume_area_0_ = 2 * physConst::PI * r_p * r_p; /* [m2], see Appendix 2 in LU2025 */
    const double plume_area_0_ref = 2 * physConst::PI * r_p_ref * r_p_ref; /* [m2], see Appendix 2 in LU2025 */

    /* Temperature - 205 K*/
    const double T_205 = temperature_K - 205.0; /* [K], from Eq. A3 in LU2025*/

    /* Height an air parcel has to descend until it just saturated when the emitted water vapor is added */
    rho_emit_ = WV_exhaust / plume_area_0_; // Eq. 6 in U2016
    const double rho_divisor = 10.; // 10 mg per m3
    z_emit_ = 1106.6 * pow(rho_emit_ * 1000. / rho_divisor, 0.678 + 0.0116 * T_205) * exp((-(0.0807+0.000428*T_205)*T_205));

    n0_ref_ = N0_ref / plume_area_0_ref;

    // Eq. 13 in U2016
    const double icenum_survfrac_h = survfrac_from_psi( 1.0 );
    double bhat;
    if ( icenum_survfrac_h <= x_s_ ) {
        bhat = eta_1_ * icenum_survfrac_h;
    }
    else {
        bhat = eta_2_ * icenum_survfrac_h + (eta_1_ - eta_2_) * x_s_;
    }
    depth_mature_ = bhat * z_desc_; // Eq. 12 in U2016

    // Rectangle-equivalent width of the initial mature plume
    width_rect_mature_ = 0.63 * wingspan;

    // The height of the center is the maximum displacement minus half the contrail height
    z_center_ = z_desc_ - depth_mature_ / 2.0;

} /* End of Vortex::Vortex */

double Vortex::z_delta( const double psi ) const
{
    return pow(psi, gamma_exp_) * \
             (+ alpha_atm_  * z_atm_  \
              + alpha_emit_ * z_emit_) \
              - alpha_desc_ * z_desc_;

} /* End of Vortex::z_delta */

double Vortex::survfrac_from_psi( const double psi ) const
{
    const double survfrac = beta_0_ + beta_1_ / physConst::PI * atan( alpha_0_ + z_delta( psi ) / 1.0E+02 );
    return std::min( std::max( survfrac, 0.0E+00 ), 1.0E+00 );

} /* End of Vortex::survfrac_from_psi */

double Vortex::psi( const double N_postjet ) const
{
    const double n0 = N_postjet / plume_area_0_;
    const double n0_star = n0 / n0_ref_; /* [-], See Appendix A1 in LU 2025 */
    return 1 / n0_star;

} /* End of Vortex::psi */

double Vortex::icenum_survfrac( const double N_postjet ) const
{
    return survfrac_from_psi( psi( N_postjet ) );

} /* End of Vortex::icenum_survfrac */

void Vortex::print_summary( const double N_postjet, const double survfrac ) const
{
    std::cout << std::endl << std::endl;
    std::cout << "***** Vortex parametrisation START *****" << std::endl;
    std::cout << std::endl;
    std::cout << "AMBIENT PARAMS" << std::endl;
    std::cout << "Cruise Temperature  = " << temperature_K_ << " [K]" << std::endl; //
    std::cout << "Cruise RHi          = " << RHi_PC_ << " [%]" << std::endl; //
    std::cout << "N_BV                = " << N_BV_ << " [1/s]" << std::endl; //
    std::cout << std::endl;
    std::cout << "VORTEX PARAMS" << std::endl;
    std::cout << "Reference ice count = " << N0_ref_ << " [#/m]" << std::endl; //
    std::cout << "Post-jet ice count  = " << N_postjet << " [#/m]" << std::endl; //
    std::cout << "Exhaust Water Vapor = " << WV_exhaust_ << " [g/m]" << std::endl; //
    std::cout << "Density of emitted water vapor: " << rho_emit_ << " [kg/m^3]" << std::endl;
    std::cout << "gamma               = " << gamma_ << " [m^2/s]" << std::endl;
    std::cout << std::endl;
    std::cout << "AIRCRAFT PARAMS" << std::endl;
    std::cout << "wingspan            = " << wingspan_ << " [m]" << std::endl; //
    std::cout << "reference wingspan  = " << wingspan_ref_ << " [m]" << std::endl; //
    std::cout << "flight speed        = " << vFlight_ << " [m/s]" << std::endl; //
    std::cout << std::endl;
    std::cout << "PARAMETRISATION RESULTS (SURVIVAL FRACTION)" << std::endl;
    std::cout << "z_desc              = " << z_desc_ << " [m]" << std::endl; //
    std::cout << "z_atm               = " << z_atm_ << " [m]" << std::endl; //
    std::cout << "z_emit              = " << z_emit_ << " [m]" << std::endl; //
    std::cout << "z_delta (survfrac)  = " << z_delta( psi( N_postjet ) ) << " [m]" << std::endl; //
    std::cout << "Survival Fraction   = " << survfrac << std::endl; //
    if ( survfrac == 0.0E+00 )
        std::cout << "Contrail has fully melted because of vortex-sinking losses" << std::endl;
    std::cout << std::endl;
    std::cout << "PARAMETRISATION RESULTS (GEOMETRY)" << std::endl;
    std::cout << "z_delta (depth)     = " << z_delta( 1.0 ) << " [m]" << std::endl; //
    std::cout << "Survfrac (depth)    = " << survfrac_from_psi( 1.0 ) << std::endl; //
    std::cout << "Contrail Depth      = " << depth_mature_ << " [m]" << std::endl; //
    std::cout << "Contrail Area Width = " << width_rect_mature_ << " [m]" << std::endl; //
    std::cout << "Contrail Area       = " << width_rect_mature_ * depth_mature_ << " [m*2]" << std::endl; //
    std::cout << "Contrail Center y   = " << z_center_ << " [m]" << std::endl; //
    std::cout << std::endl;
    std::cout << "***** Vortex parametrisation END *****" << std::endl;
    std::cout << std::endl << std::endl;

} /* End of Vortex::print_summary */


/* End of Vortex.cpp */
