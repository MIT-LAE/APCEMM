/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* BuildMesh Program File                                           */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : BuildMesh.cpp                             */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void BuildMesh( double *x, double *y, double const xlim, double const ylim, unsigned int const nx, unsigned int const ny )
{

    double const hx = 2 * xlim / nx;
    double const hy = 2 * ylim / ny;

    for ( unsigned int i = 0; i < nx; i++ )
        x[i] = i * hx - xlim + hx / 2;

    for ( unsigned int j = 0; j < ny; j++ )
        y[j] = j * hy - ylim + hy / 2;

} /* Enf of BuildMesh */


