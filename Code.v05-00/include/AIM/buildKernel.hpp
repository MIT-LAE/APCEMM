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

#include <iostream>
#include <vector>
#include <fstream>

#include "Util/ForwardDecl.hpp"
#include "Util/PhysConstant.hpp"
#include "Util/PhysFunction.hpp"

namespace AIM
{
   
    /* Brownian coagulation kernel */
    Vector_1D buildBrownianKernel( RealDouble temperature_K, RealDouble pressure_Pa, Vector_1D const &bin_Centers, RealDouble rho_1 , RealDouble bin_R, RealDouble rho_2 );
    Vector_2D buildBrownianKernel( RealDouble temperature_K, RealDouble pressure_Pa, Vector_1D const &bin_Centers_1, RealDouble rho_1 , Vector_1D const &bin_Centers_2, RealDouble rho_2 );

    /* Convective Brownian diffusion enhancement */
    Vector_1D buildDEKernel( RealDouble temperature_K, RealDouble pressure_Pa, Vector_1D const &bin_Centers, RealDouble rho_1, RealDouble bin_R, RealDouble rho_2, Vector_1D const &K_Brow );
    Vector_2D buildDEKernel( RealDouble temperature_K, RealDouble pressure_Pa, Vector_1D const &bin_Centers_1, RealDouble rho_1, Vector_1D const &bin_Centers_2, RealDouble rho_2, Vector_2D const &K_Brow );

    /* Gravitational collection kernel */
    Vector_1D buildGCKernel( RealDouble temperature_K, RealDouble pressure_Pa, Vector_1D const &bin_Centers, RealDouble rho_1, RealDouble bin_R, RealDouble rho_2 );
    Vector_2D buildGCKernel( RealDouble temperature_K, RealDouble pressure_Pa, Vector_1D const &bin_Centers_1, RealDouble rho_1, Vector_1D const &bin_Centers_2, RealDouble rho_2 );

    /* Turbulent initial motion */
    Vector_1D buildTIKernel( RealDouble temperature_K, RealDouble pressure_Pa, Vector_1D const &bin_Centers, RealDouble rho_1, RealDouble bin_R, RealDouble rho_2 );
    Vector_2D buildTIKernel( RealDouble temperature_K, RealDouble pressure_Pa, Vector_1D const &bin_Centers_1, RealDouble rho_1, Vector_1D const &bin_Centers_2, RealDouble rho_2 );

    /* Turbulent shear */
    Vector_1D buildTSKernel( RealDouble temperature_K, RealDouble pressure_Pa, Vector_1D const &bin_Centers, RealDouble rho_1, RealDouble bin_R, RealDouble rho_2 );
    Vector_2D buildTSKernel( RealDouble temperature_K, RealDouble pressure_Pa, Vector_1D const &bin_Centers_1, RealDouble rho_1, Vector_1D const &bin_Centers_2, RealDouble rho_2 );

    /* Print kernel to file for validation */
    void printKernel2File( Vector_1D Kernel, const char* fileName = "/home/fritzt/CAPCEMM/AIM/data/Kernel.out" );
    void printKernel2File( Vector_2D Kernel, const char* fileName = "/home/fritzt/CAPCEMM/AIM/data/Kernel.out" );
    
}

#endif /* BUILDKERNEL_H_INCLUDED */
