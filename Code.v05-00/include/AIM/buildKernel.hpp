/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*                        AIrcraft Microphysics                     */
/*                              (AIM)                               */
/*                                                                  */
/* buildKernel Header File                                          */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 9/27/2018                                 */
/* File                 : buildKernel.hpp                           */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef BUILDKERNEL_H_INCLUDED
#define BUILDKERNEL_H_INCLUDED

#include "Util/ForwardDecl.hpp"

namespace AIM
{
   
    /* Brownian coagulation kernel */
    Vector_1D buildBrownianKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers, double rho_1 , double bin_R, double rho_2 );
    Vector_2D buildBrownianKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers_1, double rho_1 , Vector_1D const &bin_Centers_2, double rho_2 );

    /* Convective Brownian diffusion enhancement */
    Vector_1D buildDEKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers, double rho_1, double bin_R, double rho_2, Vector_1D const &K_Brow );
    Vector_2D buildDEKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers_1, double rho_1, Vector_1D const &bin_Centers_2, double rho_2, Vector_2D const &K_Brow );

    /* Gravitational collection kernel */
    Vector_1D buildGCKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers, double rho_1, double bin_R, double rho_2 );
    Vector_2D buildGCKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers_1, double rho_1, Vector_1D const &bin_Centers_2, double rho_2 );

    /* Turbulent initial motion */
    Vector_1D buildTIKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers, double rho_1, double bin_R, double rho_2 );
    Vector_2D buildTIKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers_1, double rho_1, Vector_1D const &bin_Centers_2, double rho_2 );

    /* Turbulent shear */
    Vector_1D buildTSKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers, double rho_1, double bin_R, double rho_2 );
    Vector_2D buildTSKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers_1, double rho_1, Vector_1D const &bin_Centers_2, double rho_2 );
    
}

#endif /* BUILDKERNEL_H_INCLUDED */
