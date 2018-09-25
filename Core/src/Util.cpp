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

    double* vect2double( std::vector<std::vector<double> > &vals, unsigned int N, unsigned int M, double scalingFactor )
    {
        double* temp;
        temp = new double[N*M];

        for( unsigned int i = 0; i < N; i++ ) {
            for( unsigned int j = 0; j < M; j++ )
                temp[j + M*i] = vals[i][j] * scalingFactor;
        }

        return temp;
    }
    
    double* vect2double( std::vector<double> &vals, unsigned int N, double scalingFactor )
    {
        double* temp;
        temp = new double[N];

        for( unsigned int i = 0; i < N; i++ ) {
            temp[i] = vals[i] * scalingFactor;
        }

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
    
    float* vect2float( std::vector<std::vector<double> > &vals, unsigned int N, unsigned int M, double scalingFactor )
    {
        float* temp;
        temp = new float[N*M];

        for( unsigned int i = 0; i < N; i++ ) {
            for( unsigned int j = 0; j < M; j++ )
                temp[j + M*i] = (float) vals[i][j] * scalingFactor;
        }

        return temp;
    }
    
    float* vect2float( std::vector<double> &vals, unsigned int N, double scalingFactor )
    {
        float* temp;
        temp = new float[N];

        for( unsigned int i = 0; i < N; i++ ) {
            temp[i] = (float) vals[i] * scalingFactor;
        }

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

    std::vector<std::vector<double> > Array2Vect( double** A, unsigned int N, unsigned int M, double scalingFactor )
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

template void delete1D<double>( double* temp );
template void delete1D<float>( float* temp );

}

/* End of Util.cpp */
