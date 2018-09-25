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

#include <iostream>
#include <vector>

namespace util
{
    double** EW_Multiply( double** A, double** B, unsigned int N, unsigned int M );
    double* vect2double ( std::vector<std::vector<double> > &vals, unsigned int N, unsigned int M, double scalingFactor = 1.0 );
    double* vect2double ( std::vector<double> &vals, unsigned int N, double scalingFactor = 1.0 );
    float* vect2float ( std::vector<std::vector<double> > &vals, unsigned int N, unsigned int M, double scalingFactor = 1.0 );
    float* vect2float ( std::vector<double> &vals, unsigned int N, double scalingFactor = 1.0 );
    std::vector<std::vector<double> > Array2Vect( double** A, unsigned int N, unsigned int M, double scalingFactor = 1.0 );
    void PrintVector( std::vector<std::vector<double> > Array );
    template <class T>
    void delete1D( T* temp );
    template <class T>
    void delete2D( T** temp, unsigned int N );

}

#endif /* UTIL_H_INCLUDED */
