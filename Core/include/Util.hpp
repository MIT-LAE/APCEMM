/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Util Header File                                                 */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Util.h                                    */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <vector>

namespace util
{
    double** EW_Multiply( double** A, double** B, unsigned int N, unsigned int M );
    double** Vect2Array( std::vector<std::vector<double> > &vals, unsigned int N, unsigned int M );
    std::vector<std::vector<double> > Array2Vect( double** A, unsigned int N, unsigned int M);
    void PrintVector( std::vector<std::vector<double> > Array );

}
