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

#include "PhysConstant.hpp"

class Vortex
{

    public:

        Vortex( );
        Vortex( double temperature_K, double pressure_Pa, double n_BV, double span, double mass, double vFlight );
        ~Vortex( );
        Vortex( const Vortex &v );
        Vortex& operator=( const Vortex &v ); 
        double getb() const;
        double getgamma() const;
        double gett() const;
        double getw() const;
        double geteps() const;
        double getdeltazw() const;
        double getdeltaz1() const;
        double getD1() const;

    protected:

        double b; /* wake vortex separation
                   * Unit : m */
        double gamma; /* initial circulation 
                       * Unit : m ^ 2 / s */ 
        double t; /* effective time scale
                   * Unit : s */
        double w; /* initial velocity scale
                   * Unit : m / s */
        double eps_star; /* normalized dissipation rate
                          * Unit : - */
        double delta_zw; /* maximum downwash displacement
                          * Unit : m */
        double delta_z1; /* mean downwash displacement
                          * Unit : m */
        double D_1; /* initial contrail depth 
                     * Unit : m */

    private:

        const double Cz1  = 0.25;
        const double CD_0 = 0.5;
        const double n_BVt_threshold = 0.8;  /* Unit : - */
        const double eps_threshold = 0.36;   /* Unit : - */
        const double delta_zw_default = 200; /* Unit : m */

};


#endif /* VORTEX_H_INCLUDED */

