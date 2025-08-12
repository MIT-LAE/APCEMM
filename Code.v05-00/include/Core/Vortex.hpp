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

// #include "Util/ForwardDecl.hpp"

class Vortex
{
        
    public:
        constexpr static double Cz1  = 0.25;
        constexpr static double CD_0 = 0.5;
        constexpr static double N_BVt_threshold = 0.8;  /* Unit : - */
        constexpr static double eps_threshold = 0.36;   /* Unit : - */
        constexpr static double z_desc_default = 200; /* Unit : m */

        /* Constructors */
        Vortex( ) = default;
        Vortex( double RHi_PC, double temperature_K, double pressure_Pa,  \
                double N_BV, double wingspan, double ac_mass, \
                double vFlight, double WV_exhaust, double N0, \
                double N0_ref = 3.38E12);

        /* Getters: */

        /* Fitting coefficients, 
         * Eqs. 13a to 13g in 
         * Lottermoser and Unterstrasser (2025)*/
        double beta_0() const {return beta_0_};
        double beta_1()  const {return beta_1_};
        double alpha_0() const {return alpha_0_};
        double alpha_atm() const {return alpha_atm_};
        double alpha_emit() const {return alpha_emit_};
        double alpha_desc() const {return alpha_desc_};
        double gamma_exp() const {return gamma_exp_};

        /* Brunt-Väisala frequency */
        double N_BV() const { return N_BV_; }

        /* Wake vortex separation */
        double b() const { return b_; }

        /* Initial circulation */
        double gamma() const { return gamma_; }

        /* Effective time scale */
        double t() const { return t_; }

        /* Initial velocity scale */
        double w() const { return w_; }

        /* Temperature at cruise minus 205 K */
        double T_205() const { return T_205_; }

        /* Normalised ice number */
        double n0_star() const { return n0_star_; }

        /* Plume radius */
        double r_p() const { return r_p_; }

        /* Plume area before vortex breakup */
        double plume_area_0() const { return plume_area_0_; }

        /* Degree of supersaturation*/
        double s() const { return s_; }

        /* Density of the emitted water vapor */
        double rho_emit() const { return rho_emit_; }

        /* Maximum downwash displacement */
        double z_desc() const { return z_desc_; }

        /* Height an air parcel has to descend 
         * until it is no longer supersaturated
         * (due to adiabatic heating only) */
        double z_atm() const { return z_atm_; }

        /* "Vertical displacement that corresponds 
         * to an adiabatic heating such that an initially 
         * saturated parcel remains at saturation when 
         * the emitted water vapor is added to the parcel"
         * (Lottermoser and Unterstrasser, 2025)*/
        double z_emit() const { return z_emit_; }

        /* Linear combination of the length scales
         * (for the survival fraction)
         * in the EPM */
        double z_delta_fns() const { return z_delta_fns_; }

        /* Linear combination of the length scales
         * (for the parametrised contrail height)
         * in the EPM */
        double z_delta_h() const { return z_delta_h_; }

        /* Ice number survival fraction */
        double icenum_survfrac() const { return icenum_survfrac_; }

        /* Downwash displacement of the early plume center */
        double z_center() const { return z_center_; }

        /* Initial contrail depth */
        double D1() const { return D_1_; }

    protected:

        /* Fitting coefficients, 
         * Eqs. 13a to 13g in 
         * Lottermoser and Unterstrasser (2025)*/
        const double beta_0_  = +0.42;
        const double beta_1_  = +1.31;
        const double alpha_0_ = -1.00;
        const double alpha_atm_ = +1.27;
        const double alpha_emit_ = +0.42;
        const double alpha_desc_ = +0.49;
        const double gamma_exp_ = +0.16;

        /* Brunt-Väisala frequency
         * Unit: s^-1 */
        double N_BV_;

        /* Wake vortex separation 
         * Unit: m */
        double b_;
        
        /* Initial circulation
         * Unit: m^2/s */
        double gamma_;
        
        /* Effective time scale
         * Unit: s */
        double t_;
        
        /* Initial velocity scale
         * Unit: s */
        double w_;

        /* Temperature at cruise altitude
         * minus 205 K (see Eq. A3 in
         * Lottermoser and Unterstrasser, 2025)
         * Unit: K */
        double T_205_;

        /* Normalised ice number
         * Unit: - */
        double n0_star_;

        /* Plume radius
         * Unit: m */
        double r_p_;

        /* Plume area before vortex breakup
         * Unit: m^2 */
        double plume_area_0_;

        /* Degree of supersaturation
         * see S2 in Unterstrasser (2016)
         * Unit: - */
        double s_;

        /* Density of the emitted water vapor
         * Unit: kg / m^3 */
        double rho_emit_;

        /* Maximum downwash displacement,
         * equal to the height delta between the
         * release altitude and the bottom of the
         * early plume after the vortex phase.
         * See Fig.3 in Unterstrasser and Görsch (2014) 
         * Unit: m */
        double z_desc_;
 
        /* Height an air parcel has to descend 
         * until it is no longer supersaturated
         * (due to adiabatic heating only)
         * see Lottermoser and Unterstrasser (2025)
         * Unit: - */
        double z_atm_;

        /* "Vertical displacement that corresponds 
         * to an adiabatic heating such that an initially 
         * saturated parcel remains at saturation when 
         * the emitted water vapor is added to the parcel"
         * (Lottermoser and Unterstrasser, 2025)
         * Unit: - */
        double z_emit_;
            
        /* Linear combination of the length scales
         * (for the survival fraction)
         * Unit: m */
        double z_delta_fns_;

        /* Linear combination of the length scales
         * (for the parametrised contrail height)
         * Unit: m */
        double z_delta_h_;

        /* Ice number survival fraction
         * Unit: - */
        double icenum_survfrac_;
        
        /* Displacement of the early plume center
         * Unit: m */
        double z_center_;
        
        /* Initial contrail depth
         * Unit: m */
        double D_1_;


};


#endif /* VORTEX_H_INCLUDED */

