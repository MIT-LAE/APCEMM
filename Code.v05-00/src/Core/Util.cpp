/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Util Program File                                                */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Util.cpp                                  */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Core/Util.hpp"

namespace util
{
    double* vect2double( const std::vector<std::vector<double>> &vals, unsigned int N, unsigned int M )
    {
        double* temp;
        temp = new double[N*M];

        for( unsigned int n = 0; n < N; n++ ) {
            for( unsigned int m = 0; m < M; m++ )
                temp[m + M*n] = vals[n][m];
        }

        return temp;
    }

    float* vect2float( const std::vector<std::vector<double>> &vals, unsigned int N, unsigned int M )
    {
        float* temp;
        temp = new float[N*M];

        for( unsigned int n = 0; n < N; n++ ) {
            for( unsigned int m = 0; m < M; m++ )
                temp[m + M*n] = (float) vals[n][m];
        }

        return temp;
    }

    float* vect2float( const std::vector<double> &vals, unsigned int N )
    {
        float* temp;
        temp = new float[N];

        for( unsigned int n = 0; n < N; n++ )
            temp[n] = (float) vals[n];

        return temp;
    }
}

/* End of Util.cpp */
