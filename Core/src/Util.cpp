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
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Util.hpp"

namespace util
{
    double** EW_Multiply( double** A, double** B, unsigned int N, unsigned int M )
    {
        double** C;
        C = new double*[N];

        for ( unsigned int i = 0; i < N; i++ ) {
            C[i] = new double[M];
            for ( unsigned int j = 0; j < M; j++ )
                C[i][j] = A[i][j]*B[i][j];
        }

        return C;
    }

    double** Vect2Array( std::vector<std::vector<double> > &vals, unsigned int N, unsigned int M )
    {
        double** temp;
        temp = new double*[N];

        for( unsigned int i=0; i < N; i++ ) {
            temp[i] = new double[M];
            for( unsigned int j=0; j < M; j++ )
                temp[i][j] = vals[i][j];
        }

        return temp;
    }

    std::vector<std::vector<double> > Array2Vect( double** A, unsigned int N, unsigned int M)
    {
        std::vector<std::vector<double> > vect;

        for ( unsigned int i = 0; i < N; i++ ) {
            vect.push_back(std::vector<double>( M ) );
            for ( unsigned int j = 0; j < M; j++ )
                vect[i][j] = A[i][j];
        }

        return vect;
    }

//    std::complex<double>** Vect2Array_Complex( std::vector<std::vector<std::complex<double>> > &vals, unsigned int N, unsigned int M )
//    {
//        double** temp;
//        temp = new double*[N];
//
//        for( unsigned i=0; i < N; i++ ) {
//            temp[i] = new double[M];
//            for( unsigned j=0; j < M; j++ )
//                temp[i][j] = vals[i][j];
//        }
//
//        return temp;
//    }

    void Print2DVector( std::vector<std::vector<double> > Array )
    {
        for ( unsigned int i = 0; i < Array.size(); i++ ) {
            for ( unsigned int j = 0; j < Array.size(); j++ ) {
                std::cout << Array[i][j] << " ";
            }
            std::cout << "" << std::endl;
        }
    }


}
