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

#include <algorithm>
#include <iostream>
#include "Core/Util.hpp"

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

    double* vect2double( const std::vector<std::vector<std::vector<std::vector<double>>>> &vals, unsigned int N, unsigned int M, \
                         unsigned int L, unsigned int K, double scalingFactor )
    {
        double* temp;
        temp = new double[N*M*L*K];

        for( unsigned int n = 0; n < N; n++ ) {
            for( unsigned int m = 0; m < M; m++ ) {
                for ( unsigned int l = 0; l < L; l++ ) {
                    for ( unsigned int k = 0; k < K; k++ ) 
                        temp[k + K*l + K*L*m + K*L*M*n] = vals[n][m][l][k] * scalingFactor;
                }
            }
        }

        return temp;
    }
    
    double* vect2double( const std::vector<std::vector<std::vector<double>>> &vals, unsigned int N, unsigned int M, \
                         unsigned int L, double scalingFactor )
    {
        double* temp;
        temp = new double[N*M*L];

        for( unsigned int n = 0; n < N; n++ ) {
            for( unsigned int m = 0; m < M; m++ ) {
                for ( unsigned int l = 0; l < L; l++ ) 
                    temp[l + L*m + L*M*n] = vals[n][m][l] * scalingFactor;
            }
        }

        return temp;
    }
    
    double* vect2double( const std::vector<std::vector<double>> &vals, unsigned int N, unsigned int M, \
                         double scalingFactor )
    {
        double* temp;
        temp = new double[N*M];

        for( unsigned int n = 0; n < N; n++ ) {
            for( unsigned int m = 0; m < M; m++ )
                temp[m + M*n] = vals[n][m] * scalingFactor;
        }

        return temp;
    }
    
    double* vect2double( const std::vector<double> &vals, unsigned int N, double scalingFactor )
    {
        double* temp;
        temp = new double[N];

        for( unsigned int n = 0; n < N; n++ )
            temp[n] = vals[n] * scalingFactor; 

        return temp;
    }
    
//    double** vect2double( std::vector<std::vector<double> > &vals, unsigned int N, unsigned int M, double scalingFactor )
//    {
//        double** temp;
//        temp = new double*[N];
//
//        for( unsigned int i = 0; i < N; i++ ) {
//            temp[i] = new double[M];
//            for( unsigned int j = 0; j < M; j++ )
//                temp[i][j] = vals[i][j] * scalingFactor;
//        }
//
//        return temp;
//    }
    
    float* vect2float( const std::vector<std::vector<std::vector<std::vector<double>>>> &vals, unsigned int N, unsigned int M, \
                       unsigned int L, unsigned int K, double scalingFactor )
    {
        float* temp;
        temp = new float[N*M*L*K];

        for( unsigned int n = 0; n < N; n++ ) {
            for( unsigned int m = 0; m < M; m++ ) {
                for ( unsigned int l = 0; l < L; l++ ) {
                    for ( unsigned int k = 0; k < K; k++ ) 
                        temp[k + K*l + K*L*m + K*L*M*n] = (float) vals[n][m][l][k] * scalingFactor;
                }
            }
        }

        return temp;
    }
    
    float* vect2float( const std::vector<std::vector<std::vector<double>>> &vals, unsigned int N, unsigned int M, \
                       unsigned int L, double scalingFactor )
    {
        float* temp;
        temp = new float[N*M*L];

        for( unsigned int n = 0; n < N; n++ ) {
            for( unsigned int m = 0; m < M; m++ ) {
                for ( unsigned int l = 0; l < L; l++ ) 
                    temp[l + L*m + L*M*n] = vals[n][m][l] * scalingFactor;
            }
        }

        return temp;
    }
    
    float* vect2float( const std::vector<std::vector<double>> &vals, unsigned int N, unsigned int M, \
                       double scalingFactor )
    {
        float* temp;
        temp = new float[N*M];

        for( unsigned int n = 0; n < N; n++ ) {
            for( unsigned int m = 0; m < M; m++ )
                temp[m + M*n] = (float) vals[n][m] * scalingFactor;
        }

        return temp;
    }
    
    float* vect2float( const std::vector<double> &vals, unsigned int N, double scalingFactor )
    {
        float* temp;
        temp = new float[N];

        for( unsigned int n = 0; n < N; n++ )
            temp[n] = (float) vals[n] * scalingFactor; 

        return temp;
    }

//    float** vect2float( std::vector<std::vector<double> > &vals, unsigned int N, unsigned int M, double scalingFactor )
//    {
//        float** temp;
//        temp = new float*[N];
//
//        for( unsigned int i = 0; i < N; i++ ) {
//            temp[i] = new float[M];
//            for( unsigned int j = 0; j < M; j++ )
//                temp[i][j] = float( vals[i][j] * scalingFactor );
//        }
//
//        return temp;
//    }

    template <class T>
    void delete1D( T* temp )
    {
        
        delete[] temp;
        temp = NULL;

    }
    
    template <class T>
    void delete2D( T** temp, unsigned int N )
    {
        for ( unsigned int i = 0; i < N; i++ ) {
            delete[] temp[i];
        }
        delete[] temp;
        temp = NULL;

    }

    std::vector<std::vector<double> > Array2Vect( const double** A, unsigned int N, unsigned int M, double scalingFactor )
    {
        std::vector<std::vector<double> > vect;

        for ( unsigned int i = 0; i < N; i++ ) {
            vect.push_back(std::vector<double>( M ) );
            for ( unsigned int j = 0; j < M; j++ )
                vect[i][j] = A[i][j] * scalingFactor;
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

    template <typename T>
    std::vector<T> add1D(const std::vector<T>&a, const std::vector<T>& b) 
    {
    
        std::vector<T> result;

        if ( a.size() != b.size() )
        {
            std::cout << "In Util.cpp: operator+ couldn't be performed because vectors don't have the same size\n";
            return result;
        }

        std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::plus<T>());

        return result;

    }

    template <typename T>
    std::vector<std::vector<T>> add2D(const std::vector<std::vector<T>>&a, const std::vector<std::vector<T>>& b) 
    {
    
        std::vector<std::vector<T>> result;

        if ( ( a.size() != b.size() ) )
        {
            std::cout << "In Util.cpp: operator+ couldn't be performed because vectors don't have the same size\n";
            return result;
        }

        for ( unsigned int i = 0; i < a.size(); i++ ) {
            result.push_back( std::vector<T> ( a[i].size() ) );
            for ( unsigned int j = 0; j < a[i].size(); j++ ) {
                result[i][j] = a[i][j] + b[i][j];
            }
        }

        return result;

    }

template void delete1D<double>( double* temp );
template void delete1D<float>( float* temp );
template std::vector<double> add1D(const std::vector<double>&a, const std::vector<double>& b);
template std::vector<float> add1D(const std::vector<float>&a, const std::vector<float>& b);
template std::vector<std::vector<double>> add2D(const std::vector<std::vector<double>>&a, const std::vector<std::vector<double>>& b);
template std::vector<std::vector<float>> add2D(const std::vector<std::vector<float>>&a, const std::vector<std::vector<float>>& b);

}

/* End of Util.cpp */
