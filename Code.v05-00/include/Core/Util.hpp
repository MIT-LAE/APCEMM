/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Util Header File                                                 */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Util.hpp                                  */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED

#include <vector>

namespace util
{
    double* vect2double( const std::vector<std::vector<double>> &vals, unsigned int N, unsigned int M );

    float* vect2float( const std::vector<std::vector<std::vector<double>>> &vals, unsigned int N1, unsigned int N2, unsigned int N3 );
    float* vect2float( const std::vector<std::vector<double>> &vals, unsigned int N, unsigned int M );
    float* vect2float( const std::vector<double> &vals, unsigned int N );
}

#endif /* UTIL_H_INCLUDED */
