/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*              Spectral Advection aNd Diffusion Solver             */
/*                             (SANDS)                              */
/*                                                                  */
/* FFT Program File                                                 */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 11/14/2018                                */
/* File                 : FFT.cpp                                   */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "SANDS/FFT.hpp"

FourierTransform_1D<float>::FourierTransform_1D( const unsigned int rows_ )
    :   rows( rows_ ),
        rowsC( rows_/2 + 1 ),
        fftScaling( rows_ )
{

    /* Allocate the ins and outs */
    in_FFT   = (scalar_type*)  fftwf_malloc( sizeof(scalar_type)  * rows  );
    out_FFT  = (complex_type*) fftwf_malloc( sizeof(complex_type) * rowsC );
    in_IFFT  = (complex_type*) fftwf_malloc( sizeof(complex_type) * rowsC );
    out_IFFT = (scalar_type*)  fftwf_malloc( sizeof(scalar_type)  * rows  );

    /* Create FFT plan and compute the best FFT algorithm */
    plan_FFT = fftwf_plan_dft_r2c_1d( rows, in_FFT, out_FFT, FFTW_EXHAUSTIVE );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_FFT == NULL ) {
        std::cout << " In FourierTransform_1D: Plan creation failed!\n";
        exit(-1);
    }

    /* Create IFFT plan and compute the best IFFT algorithm */
    plan_IFFT = fftwf_plan_dft_c2r_1d( rows, in_IFFT, out_IFFT, FFTW_EXHAUSTIVE );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_IFFT == NULL ) {
        std::cout << " In FourierTransform_1D: Plan creation failed!\n";
        exit(-1);
    }

} /* End of FourierTransform_1D<float>::FourierTransform_1D */

FourierTransform_1D<float>::~FourierTransform_1D()
{

    /* Free plan memory */
    fftwf_destroy_plan( plan_FFT );
    fftwf_destroy_plan( plan_IFFT );

    /* Free arrays */
    fftwf_free( in_FFT ); 
    in_FFT = NULL;
    fftwf_free( out_FFT ); 
    out_FFT = NULL;
    fftwf_free( in_IFFT ); 
    in_IFFT = NULL;
    fftwf_free( out_IFFT ); 
    out_IFFT = NULL;

} /* End of FourierTransform_1D<float>::~FourierTransform_1D */

void FourierTransform_1D<float>::Forward( scalar_type* in, complex_type* out ) const
{

    /* Computes forward DFT */
    fftwf_execute_dft_r2c( plan_FFT, in, out );

} /* End of FourierTransform_1D<float>::Forward */

void FourierTransform_1D<float>::Forward( const Vector_1Df &V, complex_type* out ) const
{

    for ( unsigned int i = 0; i < rows; i++ )
        in_FFT[i] = (scalar_type) V[i];

    /* Computes forward DFT */
    fftwf_execute_dft_r2c( plan_FFT, in_FFT, out );

} /* End of FourierTransform_1D<float>::Forward */

void FourierTransform_1D<float>::Backward( complex_type* in, scalar_type* out ) const
{

    /* Computes backward DFT */
    fftwf_execute_dft_c2r( plan_IFFT, in, out );

} /* End of FourierTransform_1D<float>::Backward */

void FourierTransform_1D<float>::Backward( complex_type* in, Vector_1Df &V ) const
{

    /* Computes backward DFT */
    fftwf_execute_dft_c2r( plan_IFFT, in, out_IFFT );

    Scale( out_IFFT );

    for ( unsigned int i = 0; i < rows; i++ )
        V[i] = out_IFFT[i];


} /* End of FourierTransform_1D<float>::Backward */

void FourierTransform_1D<float>::Scale( scalar_type* in ) const
{

    /* Scale output */
    for ( unsigned int i = 0; i < rows; i++ ) 
        in[i] /= fftScaling;

} /* End of FourierTransform_1D<float>::Scale */

void FourierTransform_1D<float>::ApplyShear( const Vector_2Dcf &shearFactor, \
                                             Vector_2Df &V ) const
{

    /* Declare array to store each horizontal layer of V */
    Vector_1Df array;
    array.assign( V[0].size(), 0.0E+00 );

    /* Dealing with each layer one at a time */
    for ( unsigned int jNy = 0; jNy < V.size(); jNy++ ) {

        /* Storing horizontal layer */
        for ( unsigned int iNx = 0; iNx < V[0].size(); iNx++ )
            array[iNx] = V[jNy][iNx];
        
        /* Apply forward transformation */
        Forward( array, out_FFT );

        /* Convolve and scale the frequencies */
        for ( unsigned int i = 0; i < rowsC; i++ ) {
            in_IFFT[i][REAL] = ( out_FFT[i][REAL] * shearFactor[jNy][i].real() \
                               - out_FFT[i][IMAG] * shearFactor[jNy][i].imag() );
            in_IFFT[i][IMAG] = ( out_FFT[i][REAL] * shearFactor[jNy][i].imag() \
                               + out_FFT[i][IMAG] * shearFactor[jNy][i].real() );
        }

        /* Apply backward transformation + scaling */
        Backward( in_IFFT, array );

        /* Apply results back into original 2D vector */
        for ( unsigned int iNx = 0; iNx < V[0].size(); iNx++ )
            V[jNy][iNx] = array[iNx];

    }

} /* End of FourierTransform_1D<float>::ApplyShear */

FourierTransform_1D<double>::FourierTransform_1D( const unsigned int rows_ )
    :   rows( rows_ ),
        rowsC( rows_/2 + 1 ),
        fftScaling( rows_ )
{

    /* Allocate the ins and outs */
    in_FFT   = (scalar_type*)  fftw_malloc( sizeof(scalar_type)  * rows  );
    out_FFT  = (complex_type*) fftw_malloc( sizeof(complex_type) * rowsC );
    in_IFFT  = (complex_type*) fftw_malloc( sizeof(complex_type) * rowsC );
    out_IFFT = (scalar_type*)  fftw_malloc( sizeof(scalar_type)  * rows  );

    /* Create FFT plan and compute the best FFT algorithm */
    plan_FFT = fftw_plan_dft_r2c_1d( rows, in_FFT, out_FFT, FFTW_EXHAUSTIVE );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_FFT == NULL ) {
        std::cout << " In FourierTransform_1D: Plan creation failed!\n";
        exit(-1);
    }

    /* Create IFFT plan and compute the best IFFT algorithm */
    plan_IFFT = fftw_plan_dft_c2r_1d( rows, in_IFFT, out_IFFT, FFTW_EXHAUSTIVE );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_IFFT == NULL ) {
        std::cout << " In FourierTransform_1D: Plan creation failed!\n";
        exit(-1);
    }

} /* End of FourierTransform_1D<double>::FourierTransform_1D */

FourierTransform_1D<double>::~FourierTransform_1D()
{

    /* Free plan memory */
    fftw_destroy_plan( plan_FFT );
    fftw_destroy_plan( plan_IFFT );

    /* Free arrays */
    fftw_free( in_FFT ); 
    in_FFT = NULL;
    fftw_free( out_FFT ); 
    out_FFT = NULL;
    fftw_free( in_IFFT ); 
    in_IFFT = NULL;
    fftw_free( out_IFFT ); 
    out_IFFT = NULL;

} /* End of FourierTransform_1D<double>::~FourierTransform_1D */

void FourierTransform_1D<double>::Forward( scalar_type* in, complex_type* out ) const
{

    /* Computes forward DFT */
    fftw_execute_dft_r2c( plan_FFT, in, out );

} /* End of FourierTransform_1D<double>::Forward */

void FourierTransform_1D<double>::Forward( const Vector_1D &V, complex_type* out ) const
{

    for ( unsigned int i = 0; i < rows; i++ )
        in_FFT[i] = (scalar_type) V[i];

    /* Computes forward DFT */
    fftw_execute_dft_r2c( plan_FFT, in_FFT, out );

} /* End of FourierTransform_1D<double>::Forward */

void FourierTransform_1D<double>::Backward( complex_type* in, scalar_type* out ) const
{

    /* Computes backward DFT */
    fftw_execute_dft_c2r( plan_IFFT, in, out );

} /* End of FourierTransform_1D<double>::Backward */

void FourierTransform_1D<double>::Backward( complex_type* in, Vector_1D &V ) const
{

    /* Computes backward DFT */
    fftw_execute_dft_c2r( plan_IFFT, in, out_IFFT );

    Scale( out_IFFT );

    for ( unsigned int i = 0; i < rows; i++ )
        V[i] = out_IFFT[i];


} /* End of FourierTransform_1D<double>::Backward */

void FourierTransform_1D<double>::Scale( scalar_type* in ) const
{

    /* Scale output */
    for ( unsigned int i = 0; i < rows; i++ ) 
        in[i] /= fftScaling;

} /* End of FourierTransform_1D<double>::Scale */

void FourierTransform_1D<double>::ApplyShear( const Vector_2Dc &shearFactor, \
                                              Vector_2D &V ) const
{

    /* Declare array to store each horizontal layer of V */
    Vector_1D array;
    array.assign( V[0].size(), 0.0E+00 );

    /* Dealing with each layer one at a time */
    for ( unsigned int jNy = 0; jNy < V.size(); jNy++ ) {

        /* Storing horizontal layer */
        for ( unsigned int iNx = 0; iNx < V[0].size(); iNx++ )
            array[iNx] = V[jNy][iNx];
        
        /* Apply forward transformation */
        Forward( array, out_FFT );

        /* Convolve and scale the frequencies */
        for ( unsigned int i = 0; i < rowsC; i++ ) {
            in_IFFT[i][REAL] = ( out_FFT[i][REAL] * shearFactor[jNy][i].real() \
                               - out_FFT[i][IMAG] * shearFactor[jNy][i].imag() );
            in_IFFT[i][IMAG] = ( out_FFT[i][REAL] * shearFactor[jNy][i].imag() \
                               + out_FFT[i][IMAG] * shearFactor[jNy][i].real() );
        }

        /* Apply backward transformation + scaling */
        Backward( in_IFFT, array );

        /* Apply results back into original 2D vector */
        for ( unsigned int iNx = 0; iNx < V[0].size(); iNx++ )
            V[jNy][iNx] = array[iNx];

    }

} /* End of FourierTransform_1D<double>::ApplyShear */

FourierTransform_1D<long double>::FourierTransform_1D( const unsigned int rows_ )
    :   rows( rows_ ),
        rowsC( rows_/2 + 1 ),
        fftScaling( rows_ )
{

    /* Allocate the ins and outs */
    in_FFT   = (scalar_type*)  fftwl_malloc( sizeof(scalar_type)  * rows  );
    out_FFT  = (complex_type*) fftwl_malloc( sizeof(complex_type) * rowsC );
    in_IFFT  = (complex_type*) fftwl_malloc( sizeof(complex_type) * rowsC );
    out_IFFT = (scalar_type*)  fftwl_malloc( sizeof(scalar_type)  * rows  );

    /* Create FFT plan and compute the best FFT algorithm */
    plan_FFT = fftwl_plan_dft_r2c_1d( rows, in_FFT, out_FFT, FFTW_EXHAUSTIVE );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_FFT == NULL ) {
        std::cout << " In FourierTransform_1D: Plan creation failed!\n";
        exit(-1);
    }

    /* Create IFFT plan and compute the best IFFT algorithm */
    plan_IFFT = fftwl_plan_dft_c2r_1d( rows, in_IFFT, out_IFFT, FFTW_EXHAUSTIVE );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_IFFT == NULL ) {
        std::cout << " In FourierTransform_1D: Plan creation failed!\n";
        exit(-1);
    }

} /* End of FourierTransform_1D<long double>::FourierTransform_1D */

FourierTransform_1D<long double>::~FourierTransform_1D()
{

    /* Free plan memory */
    fftwl_destroy_plan( plan_FFT );
    fftwl_destroy_plan( plan_IFFT );

    /* Free arrays */
    fftwl_free( in_FFT ); 
    in_FFT = NULL;
    fftwl_free( out_FFT ); 
    out_FFT = NULL;
    fftwl_free( in_IFFT ); 
    in_IFFT = NULL;
    fftwl_free( out_IFFT ); 
    out_IFFT = NULL;

} /* End of FourierTransform_1D<long double>::~FourierTransform_1D */

void FourierTransform_1D<long double>::Forward( scalar_type* in, complex_type* out ) const
{

    /* Computes forward DFT */
    fftwl_execute_dft_r2c( plan_FFT, in, out );

} /* End of FourierTransform_1D<long double>::Forward */

void FourierTransform_1D<long double>::Forward( const Vector_1Dl &V, complex_type* out ) const
{

    for ( unsigned int i = 0; i < rows; i++ )
        in_FFT[i] = (scalar_type) V[i];

    /* Computes forward DFT */
    fftwl_execute_dft_r2c( plan_FFT, in_FFT, out );

} /* End of FourierTransform_1D<long double>::Forward */

void FourierTransform_1D<long double>::Backward( complex_type* in, scalar_type* out ) const
{

    /* Computes backward DFT */
    fftwl_execute_dft_c2r( plan_IFFT, in, out );

} /* End of FourierTransform_1D<long double>::Backward */

void FourierTransform_1D<long double>::Backward( complex_type* in, Vector_1Dl &V ) const
{

    /* Computes backward DFT */
    fftwl_execute_dft_c2r( plan_IFFT, in, out_IFFT );

    Scale( out_IFFT );

    for ( unsigned int i = 0; i < rows; i++ )
        V[i] = out_IFFT[i];


} /* End of FourierTransform_1D<long double>::Backward */

void FourierTransform_1D<long double>::Scale( scalar_type* in ) const
{

    /* Scale output */
    for ( unsigned int i = 0; i < rows; i++ ) 
        in[i] /= fftScaling;

} /* End of FourierTransform_1D<long double>::Scale */

void FourierTransform_1D<long double>::ApplyShear( const Vector_2Dcl &shearFactor, \
                                                   Vector_2Dl &V ) const
{

    /* Declare array to store each horizontal layer of V */
    Vector_1Dl array;
    array.assign( V[0].size(), 0.0E+00 );

    /* Dealing with each layer one at a time */
    for ( unsigned int jNy = 0; jNy < V.size(); jNy++ ) {

        /* Storing horizontal layer */
        for ( unsigned int iNx = 0; iNx < V[0].size(); iNx++ )
            array[iNx] = V[jNy][iNx];
        
        /* Apply forward transformation */
        Forward( array, out_FFT );

        /* Convolve and scale the frequencies */
        for ( unsigned int i = 0; i < rowsC; i++ ) {
            in_IFFT[i][REAL] = ( out_FFT[i][REAL] * shearFactor[jNy][i].real() \
                               - out_FFT[i][IMAG] * shearFactor[jNy][i].imag() );
            in_IFFT[i][IMAG] = ( out_FFT[i][REAL] * shearFactor[jNy][i].imag() \
                               + out_FFT[i][IMAG] * shearFactor[jNy][i].real() );
        }

        /* Apply backward transformation + scaling */
        Backward( in_IFFT, array );

        /* Apply results back into original 2D vector */
        for ( unsigned int iNx = 0; iNx < V[0].size(); iNx++ )
            V[jNy][iNx] = array[iNx];

    }

} /* End of FourierTransform_1D<long double>::ApplyShear */

FourierTransform_2D<float>::FourierTransform_2D( const unsigned int rows_, \
                                                 const unsigned int cols_ )
    :   rows( rows_ ),
        cols( cols_ ),
        colsC( cols_/2 + 1 ),
        fftScaling( cols_ * rows_ )
{

    /* Allocate the ins and outs */
    in_FFT   = (scalar_type*)  fftwf_malloc( sizeof(scalar_type)  * rows * cols  );
    out_FFT  = (complex_type*) fftwf_malloc( sizeof(complex_type) * rows * colsC );
    in_IFFT  = (complex_type*) fftwf_malloc( sizeof(complex_type) * rows * colsC );
    out_IFFT = (scalar_type*)  fftwf_malloc( sizeof(scalar_type)  * rows * cols  );

    /* Create FFT plan and compute the best FFT algorithm */
    plan_FFT = fftwf_plan_dft_r2c_2d( rows, cols, in_FFT, out_FFT, FFTW_EXHAUSTIVE );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_FFT == NULL ) {
        std::cout << " In FourierTransform_2D: Plan creation failed!\n";
        exit(-1);
    }

    /* Create IFFT plan and compute the best IFFT algorithm */
    plan_IFFT = fftwf_plan_dft_c2r_2d( rows, cols, in_IFFT, out_IFFT, FFTW_EXHAUSTIVE );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_IFFT == NULL ) {
        std::cout << " In FourierTransform_2D: Plan creation failed!\n";
        exit(-1);
    }

} /* End of FourierTransform_2D<float>::FourierTransform_2D */

FourierTransform_2D<float>::~FourierTransform_2D()
{

    /* Free plan memory */
    fftwf_destroy_plan( plan_FFT );
    fftwf_destroy_plan( plan_IFFT );

    /* Free arrays */
    fftwf_free( in_FFT ); 
    in_FFT = NULL;
    fftwf_free( out_FFT ); 
    out_FFT = NULL;
    fftwf_free( in_IFFT ); 
    in_IFFT = NULL;
    fftwf_free( out_IFFT ); 
    out_IFFT = NULL;

} /* End of FourierTransform_2D<float>::~FourierTransform_2D */

void FourierTransform_2D<float>::Forward( scalar_type* in, complex_type* out ) const
{

    /* Computes forward DFT */
    fftwf_execute_dft_r2c( plan_FFT, in, out );

} /* End of FourierTransform_2D<float>::Forward */

void FourierTransform_2D<float>::Forward( const Vector_2Df &V, complex_type* out ) const
{

    for ( unsigned int i = 0; i < rows; i++ ) {
        for ( unsigned int j = 0; j < cols; j++ )
            in_FFT[i * cols + j] = (scalar_type) V[j][i];
    }

    /* Computes forward DFT */
    fftwf_execute_dft_r2c( plan_FFT, in_FFT, out );

} /* End of FourierTransform_2D<float>::Forward */

void FourierTransform_2D<float>::Backward( complex_type* in, scalar_type* out ) const
{

    /* Computes backward DFT */
    fftwf_execute_dft_c2r( plan_IFFT, in, out );

} /* End of FourierTransform_2D<float>::Backward */

void FourierTransform_2D<float>::Backward( complex_type* in, Vector_2Df &V ) const
{

    /* Computes backward DFT */
    fftwf_execute_dft_c2r( plan_IFFT, in, out_IFFT );

    Scale( out_IFFT );

    for ( unsigned int i = 0; i < rows; i++ ) {
        for ( unsigned int j = 0; j < cols; j++ )
            V[j][i] = out_IFFT[i * cols + j];
    }


} /* End of FourierTransform_2D<float>::Backward */

void FourierTransform_2D<float>::Scale( scalar_type* in ) const
{

    /* Scale output */
    for ( unsigned int i = 0; i < cols * rows; i++ ) 
        in[i] /= fftScaling;

} /* End of FourierTransform_2D<float>::Scale */

void FourierTransform_2D<float>::SANDS( const Vector_2Df &diffFactor, \
                                        const Vector_2Dcf &advFactor, \
                                        Vector_2Df &V ) const
{

    /* Transform input into frequency domain */
    Forward( V, out_FFT );

    /* Convolve and scale the frequencies */
    for ( unsigned int i = 0; i < rows; i++ ) {
        for ( unsigned int j = 0; j < colsC; j++ ) {
            in_IFFT[i * colsC + j][REAL] = ( out_FFT[i * colsC + j][REAL] * advFactor[j][i].real() \
                                           - out_FFT[i * colsC + j][IMAG] * advFactor[j][i].imag() ) * diffFactor[j][i] ;
            in_IFFT[i * colsC + j][IMAG] = ( out_FFT[i * colsC + j][REAL] * advFactor[j][i].imag() \
                                           + out_FFT[i * colsC + j][IMAG] * advFactor[j][i].real() ) * diffFactor[j][i] ;
        }
    }

    /* Convert back to real domain
     * Backward returns the scaled IFFT such that fft(ifft) == Id */
    Backward( in_IFFT, V );

} /* End of FourierTransform_2D<float>::SANDS */

FourierTransform_2D<double>::FourierTransform_2D( const unsigned int rows_, \
                                                  const unsigned int cols_ )
    :   rows( rows_ ),
        cols( cols_ ),
        colsC( cols_/2 + 1 ),
        fftScaling( cols_ * rows_ )
{

    /* Allocate the ins and outs */
    in_FFT   = (scalar_type*)  fftw_malloc( sizeof(scalar_type)  * rows * cols  );
    out_FFT  = (complex_type*) fftw_malloc( sizeof(complex_type) * rows * colsC );
    in_IFFT  = (complex_type*) fftw_malloc( sizeof(complex_type) * rows * colsC );
    out_IFFT = (scalar_type*)  fftw_malloc( sizeof(scalar_type)  * rows * cols  );

    /* Create FFT plan and compute the best FFT algorithm */
    plan_FFT = fftw_plan_dft_r2c_2d( rows, cols, in_FFT, out_FFT, FFTW_EXHAUSTIVE );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_FFT == NULL ) {
        std::cout << " In FourierTransform_2D: Plan creation failed!\n";
        exit(-1);
    }

    /* Create IFFT plan and compute the best IFFT algorithm */
    plan_IFFT = fftw_plan_dft_c2r_2d( rows, cols, in_IFFT, out_IFFT, FFTW_EXHAUSTIVE );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_IFFT == NULL ) {
        std::cout << " In FourierTransform_2D: Plan creation failed!\n";
        exit(-1);
    }

} /* End of FourierTransform_2D<double>::FourierTransform_2D */

FourierTransform_2D<double>::~FourierTransform_2D()
{

    /* Free plan memory */
    fftw_destroy_plan( plan_FFT );
    fftw_destroy_plan( plan_IFFT );

    /* Free arrays */
    fftw_free( in_FFT ); 
    in_FFT = NULL;
    fftw_free( out_FFT ); 
    out_FFT = NULL;
    fftw_free( in_IFFT ); 
    in_IFFT = NULL;
    fftw_free( out_IFFT ); 
    out_IFFT = NULL;

} /* End of FourierTransform_2D<double>::~FourierTransform_2D */

void FourierTransform_2D<double>::Forward( scalar_type* in, complex_type* out ) const
{

    /* Computes forward DFT */
    fftw_execute_dft_r2c( plan_FFT, in, out );

} /* End of FourierTransform_2D<double>::Forward */

void FourierTransform_2D<double>::Forward( const Vector_2D &V, complex_type* out ) const
{

    for ( unsigned int i = 0; i < rows; i++ ) {
        for ( unsigned int j = 0; j < cols; j++ )
            in_FFT[i * cols + j] = (scalar_type) V[j][i];
    }

    /* Computes forward DFT */
    fftw_execute_dft_r2c( plan_FFT, in_FFT, out );

} /* End of FourierTransform_2D<double>::Forward */

void FourierTransform_2D<double>::Backward( complex_type* in, scalar_type* out ) const
{

    /* Computes backward DFT */
    fftw_execute_dft_c2r( plan_IFFT, in, out );

} /* End of FourierTransform_2D<double>::Backward */

void FourierTransform_2D<double>::Backward( complex_type* in, Vector_2D &V ) const
{

    /* Computes backward DFT */
    fftw_execute_dft_c2r( plan_IFFT, in, out_IFFT );

    Scale( out_IFFT );

    for ( unsigned int i = 0; i < rows; i++ ) {
        for ( unsigned int j = 0; j < cols; j++ )
            V[j][i] = out_IFFT[i * cols + j];
    }


} /* End of FourierTransform_2D<double>::Backward */

void FourierTransform_2D<double>::Scale( scalar_type* in ) const
{

    /* Scale output */
    for ( unsigned int i = 0; i < cols * rows; i++ ) 
        in[i] /= fftScaling;

} /* End of FourierTransform_2D<double>::Scale */

void FourierTransform_2D<double>::SANDS( const Vector_2D &diffFactor, \
                                         const Vector_2Dc &advFactor, \
                                         Vector_2D &V ) const
{

    /* Transform input into frequency domain */
    Forward( V, out_FFT );

    /* Convolve and scale the frequencies */
    for ( unsigned int i = 0; i < rows; i++ ) {
        for ( unsigned int j = 0; j < colsC; j++ ) {
            in_IFFT[i * colsC + j][REAL] = ( out_FFT[i * colsC + j][REAL] * advFactor[j][i].real() \
                                           - out_FFT[i * colsC + j][IMAG] * advFactor[j][i].imag() ) * diffFactor[j][i] ;
            in_IFFT[i * colsC + j][IMAG] = ( out_FFT[i * colsC + j][REAL] * advFactor[j][i].imag() \
                                           + out_FFT[i * colsC + j][IMAG] * advFactor[j][i].real() ) * diffFactor[j][i] ;
        }
    }

    /* Convert back to real domain
     * Backward returns the scaled IFFT such that fft(ifft) == Id */
    Backward( in_IFFT, V );

} /* End of FourierTransform_2D<double>::SANDS */


FourierTransform_2D<long double>::FourierTransform_2D( const unsigned int rows_, \
                                                       const unsigned int cols_ )
    :   rows( rows_ ),
        cols( cols_ ),
        colsC( cols_/2 + 1 ),
        fftScaling( cols_ * rows_ )
{

    /* Allocate the ins and outs */
    in_FFT   = (scalar_type*)  fftwl_malloc( sizeof(scalar_type)  * rows * cols  );
    out_FFT  = (complex_type*) fftwl_malloc( sizeof(complex_type) * rows * colsC );
    in_IFFT  = (complex_type*) fftwl_malloc( sizeof(complex_type) * rows * colsC );
    out_IFFT = (scalar_type*)  fftwl_malloc( sizeof(scalar_type)  * rows * cols  );

    /* Create FFT plan and compute the best FFT algorithm */
    plan_FFT = fftwl_plan_dft_r2c_2d( rows, cols, in_FFT, out_FFT, FFTW_EXHAUSTIVE );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_FFT == NULL ) {
        std::cout << " In FourierTransform_2D: Plan creation failed!\n";
        exit(-1);
    }

    /* Create IFFT plan and compute the best IFFT algorithm */
    plan_IFFT = fftwl_plan_dft_c2r_2d( rows, cols, in_IFFT, out_IFFT, FFTW_EXHAUSTIVE );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_IFFT == NULL ) {
        std::cout << " In FourierTransform_2D: Plan creation failed!\n";
        exit(-1);
    }

} /* End of FourierTransform_2D<long double>::FourierTransform_2D */

FourierTransform_2D<long double>::~FourierTransform_2D()
{

    /* Free plan memory */
    fftwl_destroy_plan( plan_FFT );
    fftwl_destroy_plan( plan_IFFT );

    /* Free arrays */
    fftwl_free( in_FFT ); 
    in_FFT = NULL;
    fftwl_free( out_FFT ); 
    out_FFT = NULL;
    fftwl_free( in_IFFT ); 
    in_IFFT = NULL;
    fftwl_free( out_IFFT ); 
    out_IFFT = NULL;

} /* End of FourierTransform_2D<long double>::~FourierTransform_2D */

void FourierTransform_2D<long double>::Forward( scalar_type* in, complex_type* out ) const
{

    /* Computes forward DFT */
    fftwl_execute_dft_r2c( plan_FFT, in, out );

} /* End of FourierTransform_2D<long double>::Forward */

void FourierTransform_2D<long double>::Forward( const Vector_2Dl &V, complex_type* out ) const
{

    for ( unsigned int i = 0; i < rows; i++ ) {
        for ( unsigned int j = 0; j < cols; j++ )
            in_FFT[i * cols + j] = (scalar_type) V[j][i];
    }

    /* Computes forward DFT */
    fftwl_execute_dft_r2c( plan_FFT, in_FFT, out );

} /* End of FourierTransform_2D<long double>::Forward */

void FourierTransform_2D<long double>::Backward( complex_type* in, scalar_type* out ) const
{

    /* Computes backward DFT */
    fftwl_execute_dft_c2r( plan_IFFT, in, out );

} /* End of FourierTransform_2D<long double>::Backward */

void FourierTransform_2D<long double>::Backward( complex_type* in, Vector_2Dl &V ) const
{

    /* Computes backward DFT */
    fftwl_execute_dft_c2r( plan_IFFT, in, out_IFFT );

    Scale( out_IFFT );

    for ( unsigned int i = 0; i < rows; i++ ) {
        for ( unsigned int j = 0; j < cols; j++ )
            V[j][i] = out_IFFT[i * cols + j];
    }


} /* End of FourierTransform_2D<long double>::Backward */

void FourierTransform_2D<long double>::Scale( scalar_type* in ) const
{

    /* Scale output */
    for ( unsigned int i = 0; i < cols * rows; i++ ) 
        in[i] /= fftScaling;

} /* End of FourierTransform_2D<long double>::Scale */

void FourierTransform_2D<long double>::SANDS( const Vector_2Dl &diffFactor, \
                                              const Vector_2Dcl &advFactor, \
                                              Vector_2Dl &V ) const
{

    /* Transform input into frequency domain */
    Forward( V, out_FFT );

    /* Convolve and scale the frequencies */
    for ( unsigned int i = 0; i < rows; i++ ) {
        for ( unsigned int j = 0; j < colsC; j++ ) {
            in_IFFT[i * colsC + j][REAL] = ( out_FFT[i * colsC + j][REAL] * advFactor[j][i].real() \
                                           - out_FFT[i * colsC + j][IMAG] * advFactor[j][i].imag() ) * diffFactor[j][i] ;
            in_IFFT[i * colsC + j][IMAG] = ( out_FFT[i * colsC + j][REAL] * advFactor[j][i].imag() \
                                           + out_FFT[i * colsC + j][IMAG] * advFactor[j][i].real() ) * diffFactor[j][i] ;
        }
    }

    /* Convert back to real domain
     * Backward returns the scaled IFFT such that fft(ifft) == Id */
    Backward( in_IFFT, V );

} /* End of FourierTransform_2D<long double>::SANDS */




