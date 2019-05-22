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

#include <iostream>
#include <cmath>

#include "Util/PhysConstant.hpp"
#include "Util/ForwardDecl.hpp"

class Vortex
{
        
    public:

        /* Constructors */
        Vortex( );
        Vortex( RealDouble temperature_K, RealDouble pressure_Pa,  \
                RealDouble n_BV, RealDouble span, RealDouble mass, \
                RealDouble vFlight );

        /* Destructor */
        ~Vortex( );

        /* Copy */
        Vortex( const Vortex &v );

        /* Copy */
        Vortex& operator=( const Vortex &v ); 

        /* Getters: */

        /* Brunt-Väisala frequency */
        RealDouble N_BV() const { return N_BV_; }

        /* Wake vortex separation */
        RealDouble b() const { return b_; }

        /* Initial circulation */
        RealDouble gamma() const { return gamma_; }

        /* Effective time scale */
        RealDouble t() const { return t_; }

        /* Initial velocity scale */
        RealDouble w() const { return w_; }

        /* Normalized dissipation rate */
        RealDouble eps_star() const { return eps_star_; }

        /* Maximum downwash displacement */
        RealDouble delta_zw() const { return delta_zw_; }

        /* Mean downwash displacement */
        RealDouble delta_z1() const { return delta_z1_; }

        /* Initial contrail depth */
        RealDouble D1() const { return D_1_; }

    protected:

        /* Brunt-Väisala frequency
         * Unit: s^-1 */
        RealDouble N_BV_;

        /* Wake vortex separation 
         * Unit: m */
        RealDouble b_;
        
        /* Initial circulation
         * Unit: m^2/s */
        RealDouble gamma_;
        
        /* Effective time scale
         * Unit: s */
        RealDouble t_;
        
        /* Initial velocity scale
         * Unit: s */
        RealDouble w_;
        
        /* Normalized dissipation rate
         * Unit: - */
        RealDouble eps_star_;
        
        /* Maximum downwash displacement
         * Unit: m */
        RealDouble delta_zw_;
        
        /* Mean downwash displacement
         * Unit: m */
        RealDouble delta_z1_;
        
        /* Initial contrail depth
         * Unit: m */
        RealDouble D_1_;

    private:

        const RealDouble Cz1  = 0.25;
        const RealDouble CD_0 = 0.5;
        const RealDouble N_BVt_threshold = 0.8;  /* Unit : - */
        const RealDouble eps_threshold = 0.36;   /* Unit : - */
        const RealDouble delta_zw_default = 200; /* Unit : m */

};


#endif /* VORTEX_H_INCLUDED */

