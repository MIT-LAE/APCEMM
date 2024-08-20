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
        constexpr static double delta_zw_default = 200; /* Unit : m */

        /* Constructors */
        Vortex( ) = default;
        Vortex( double temperature_K, double pressure_Pa,  \
                double n_BV, double span, double mass, \
                double vFlight );

        /* Getters: */

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

        /* Normalized dissipation rate */
        double eps_star() const { return eps_star_; }

        /* Maximum downwash displacement */
        double delta_zw() const { return delta_zw_; }

        /* Mean downwash displacement */
        double delta_z1() const { return delta_z1_; }

        /* Initial contrail depth */
        double D1() const { return D_1_; }

    protected:

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
        
        /* Normalized dissipation rate
         * Unit: - */
        double eps_star_;
        
        /* Maximum downwash displacement
         * Unit: m */
        double delta_zw_;
        
        /* Mean downwash displacement
         * Unit: m */
        double delta_z1_;
        
        /* Initial contrail depth
         * Unit: m */
        double D_1_;


};


#endif /* VORTEX_H_INCLUDED */

