/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Vortex Header File                                               */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Vortex.hpp                                */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef VORTEX_H_INCLUDED
#define VORTEX_H_INCLUDED

/* Parametric model of the vortex phase of contrail evolution.
 *
 * The formulae are based on Unterstrasser (2016), referred to as U2016;
 * and Lottermoser & Unterstrasser (2025), referred to as LU2025.
 *
 * U2016 is available at https://doi.org/10.5194/acp-16-2059-2016
 * LU2025 is available at https://doi.org/10.5194/acp-25-7903-2025
 *
 * The constructor computes all values that are known before the EPM runs.
 * This includes the geometry of the mature plume, because U2016 makes the
 * contrail depth independent of the ice crystal number.
 * Only the survival fraction needs the post-jet ice crystal number from the
 * EPM, so it is computed later via icenum_survfrac */
class Vortex
{
    public:
        Vortex( ) = delete;

        /* INPUT:
         * - RHi_PC         : Relative humidity with respect to ice [%]
         * - temperature_K  : Ambient temperature [K]
         * - pressure_Pa    : Ambient pressure [Pa]
         * - N_BV           : Brunt-Väisälä frequency [s^-1]
         * - wingspan       : Aircraft wingspan [m]
         * - ac_mass        : Aircraft mass [kg]
         * - vFlight        : Aircraft velocity [m/s]
         * - WV_exhaust     : Exhaust water vapor [g/m]
         * - N0_ref         : Reference ice crystal number [#/m]
         * - wingspan_ref   : Reference wingspan [m] */
        Vortex( double RHi_PC, double temperature_K, double pressure_Pa,  \
                double N_BV, double wingspan, double ac_mass, \
                double vFlight, double WV_exhaust, \
                double N0_ref, double wingspan_ref );

        /* Initial sinking of contrail center [m] */
        double z_center() const { return z_center_; }

        /* Initial (parametrised) contrail depth of the mature plume [m] */
        double depth_mature() const { return depth_mature_; }

        /* Rectangle-equivalent width of initial mature plume [m] */
        double width_rect_mature() const { return width_rect_mature_; }

        /* Fraction of ice crystals that survive the vortex phase [-],
         * for a post-jet ice crystal number N_postjet [#/m] */
        double icenum_survfrac( double N_postjet ) const;

        /* Print the inputs and the results of the parametrisation */
        void print_summary( double N_postjet, double survfrac ) const;

    private:
        /* Fitting coefficients,
         * Eqs. 13a to 13g in
         * Lottermoser and Unterstrasser (2025) */
        static constexpr double beta_0_  = +0.42;
        static constexpr double beta_1_  = +1.31;
        static constexpr double alpha_0_ = -1.00;
        static constexpr double alpha_atm_ = +1.27;
        static constexpr double alpha_emit_ = +0.42;
        static constexpr double alpha_desc_ = +0.49;
        static constexpr double gamma_exp_ = +0.16;

        /* Fitting coefficients,
         * Eqs. 13 in Unterstrasser (2016) */
        static constexpr double eta_1_  = +6.00;
        static constexpr double eta_2_  = +0.15;
        static constexpr double x_s_ = +0.20;

        /* Linear combination of the length scales [m], Eq. 9 in LU2025 */
        double z_delta( double psi ) const;

        /* Ice crystal survival fraction [-], Eq. 12 in LU2025.
         * The depth uses psi = 1, which is Eq. 9 with gamma set to 0 as in U2016 */
        double survfrac_from_psi( double psi ) const;

        /* psi [-] for a post-jet ice crystal number [#/m], Eq. 10 in LU2025 */
        double psi( double N_postjet ) const;

        /* Inputs that print_summary reports */
        double RHi_PC_;
        double temperature_K_;
        double wingspan_;
        double vFlight_;
        double WV_exhaust_;
        double N0_ref_;
        double wingspan_ref_;

        /* Brunt-Väisala frequency
         * Unit: s^-1 */
        double N_BV_;

        /* Initial circulation
         * Unit: m^2/s */
        double gamma_;

        /* Plume area before vortex breakup
         * Unit: m^2 */
        double plume_area_0_;

        /* Reference ice number concentration
         * Unit: #/m^3 */
        double n0_ref_;

        /* Density of the emitted water vapor
         * Unit: g / m^3 */
        double rho_emit_;

        /* Maximum downwash displacement, Eq. 5 in Lottermosser and Unterstrasser (2025)
         * Unit: m */
        double z_desc_;

        /* Height an air parcel has to descend
         * until it is no longer supersaturated
         * (due to adiabatic heating only)
         * see Lottermoser and Unterstrasser (2025)
         * Unit: m */
        double z_atm_;

        /* "Vertical displacement that corresponds
         * to an adiabatic heating such that an initially
         * saturated parcel remains at saturation when
         * the emitted water vapor is added to the parcel"
         * (Lottermoser and Unterstrasser, 2025)
         * Unit: m */
        double z_emit_;

        /* Initial (parametrised) contrail depth
         * of the mature plume
         * Unit: m */
        double depth_mature_;

        /* Initial (parametrised) contrail width
         * of the mature plume
         * Unit: m */
        double width_rect_mature_;

        /* Displacement of the early plume center
         * Unit: m */
        double z_center_;
};


#endif /* VORTEX_H_INCLUDED */

