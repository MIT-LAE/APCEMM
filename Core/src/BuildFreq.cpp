/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* BuildFreq Program File                                           */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : BuildFreq.cpp                             */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <cmath>

#include "PhysConstant.hpp"

void BuildFreq( double *kx, double *ky, double *kxx, double *kyy, \
                double const xlim, double const ylim, unsigned int const nx, unsigned int const ny )
{

    unsigned int i0;
    int k;

    i0 = nx/2;
    for ( unsigned int i = 0; i < nx; i++ ) {
        k = (i0%nx) - nx/2;
        kx[i] = 2.0 * PI / xlim * k;
        kxx[i] = - kx[i] * kx[i];
        i0++;
    }

    i0 = ny/2;
    for ( unsigned int j = 0; j < ny; j++ ) {
        k = (i0%ny) - ny/2;
        ky[j] = 2.0 * PI / ylim * k;
        kyy[j] = - ky[j] * ky[j];
        i0++;
    }

} /* End of BuildFreq */


