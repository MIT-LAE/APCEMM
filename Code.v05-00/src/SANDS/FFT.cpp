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

FourierTransform<float>::FourierTransform( const unsigned int rows_, \
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
        std::cout << " In FourierTransform: Plan creation failed!\n";
        exit(-1);
    }

    /* Create IFFT plan and compute the best IFFT algorithm */
    plan_IFFT = fftwf_plan_dft_c2r_2d( rows, cols, in_IFFT, out_IFFT, FFTW_EXHAUSTIVE );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_IFFT == NULL ) {
        std::cout << " In FourierTransform: Plan creation failed!\n";
        exit(-1);
    }

} /* End of FourierTransform<float>::FourierTransform */

FourierTransform<float>::~FourierTransform()
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

} /* End of FourierTransform<float>::~FourierTransform */

void FourierTransform<float>::Forward( scalar_type* in, complex_type* out ) const
{

    /* Computes forward DFT */
    fftwf_execute_dft_r2c( plan_FFT, in, out );

} /* End of FourierTransform<float>::Forward */

void FourierTransform<float>::Forward( const Vector_2Df &V, complex_type* out ) const
{

    for ( unsigned int i = 0; i < rows; i++ ) {
        for ( unsigned int j = 0; j < cols; j++ )
            in_FFT[i * cols + j] = (scalar_type) V[j][i];
    }

    /* Computes forward DFT */
    fftwf_execute_dft_r2c( plan_FFT, in_FFT, out );

} /* End of FourierTransform<float>::Forward */

void FourierTransform<float>::Backward( complex_type* in, scalar_type* out ) const
{

    /* Computes backward DFT */
    fftwf_execute_dft_c2r( plan_IFFT, in, out );

} /* End of FourierTransform<float>::Backward */

void FourierTransform<float>::Backward( complex_type* in, Vector_2Df &V ) const
{

    /* Computes backward DFT */
    fftwf_execute_dft_c2r( plan_IFFT, in, out_IFFT );

    Scale( out_IFFT );

    for ( unsigned int i = 0; i < rows; i++ ) {
        for ( unsigned int j = 0; j < cols; j++ )
            V[j][i] = out_IFFT[i * cols + j];
    }


} /* End of FourierTransform<float>::Backward */

void FourierTransform<float>::Scale( scalar_type* in ) const
{

    /* Scale output */
    for ( unsigned int i = 0; i < cols * rows; i++ ) 
        in[i] /= fftScaling;

} /* End of FourierTransform<float>::Scale */

void FourierTransform<float>::SANDS( const Vector_2Df &diffFactor, \
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

} /* End of FourierTransform<float>::SANDS */

FourierTransform<double>::FourierTransform( const unsigned int rows_, \
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
        std::cout << " In FourierTransform: Plan creation failed!\n";
        exit(-1);
    }

    /* Create IFFT plan and compute the best IFFT algorithm */
    plan_IFFT = fftw_plan_dft_c2r_2d( rows, cols, in_IFFT, out_IFFT, FFTW_EXHAUSTIVE );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_IFFT == NULL ) {
        std::cout << " In FourierTransform: Plan creation failed!\n";
        exit(-1);
    }

} /* End of FourierTransform<double>::FourierTransform */

FourierTransform<double>::~FourierTransform()
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

} /* End of FourierTransform<double>::~FourierTransform */

void FourierTransform<double>::Forward( scalar_type* in, complex_type* out ) const
{

    /* Computes forward DFT */
    fftw_execute_dft_r2c( plan_FFT, in, out );

} /* End of FourierTransform<double>::Forward */

void FourierTransform<double>::Forward( const Vector_2D &V, complex_type* out ) const
{

    for ( unsigned int i = 0; i < rows; i++ ) {
        for ( unsigned int j = 0; j < cols; j++ )
            in_FFT[i * cols + j] = (scalar_type) V[j][i];
    }

    /* Computes forward DFT */
    fftw_execute_dft_r2c( plan_FFT, in_FFT, out );

} /* End of FourierTransform<double>::Forward */

void FourierTransform<double>::Backward( complex_type* in, scalar_type* out ) const
{

    /* Computes backward DFT */
    fftw_execute_dft_c2r( plan_IFFT, in, out );

} /* End of FourierTransform<double>::Backward */

void FourierTransform<double>::Backward( complex_type* in, Vector_2D &V ) const
{

    /* Computes backward DFT */
    fftw_execute_dft_c2r( plan_IFFT, in, out_IFFT );

    Scale( out_IFFT );

    for ( unsigned int i = 0; i < rows; i++ ) {
        for ( unsigned int j = 0; j < cols; j++ )
            V[j][i] = out_IFFT[i * cols + j];
    }


} /* End of FourierTransform<double>::Backward */

void FourierTransform<double>::Scale( scalar_type* in ) const
{

    /* Scale output */
    for ( unsigned int i = 0; i < cols * rows; i++ ) 
        in[i] /= fftScaling;

} /* End of FourierTransform<double>::Scale */

void FourierTransform<double>::SANDS( const Vector_2D &diffFactor, \
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

} /* End of FourierTransform<double>::SANDS */


FourierTransform<long double>::FourierTransform( const unsigned int rows_, \
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
        std::cout << " In FourierTransform: Plan creation failed!\n";
        exit(-1);
    }

    /* Create IFFT plan and compute the best IFFT algorithm */
    plan_IFFT = fftwl_plan_dft_c2r_2d( rows, cols, in_IFFT, out_IFFT, FFTW_EXHAUSTIVE );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_IFFT == NULL ) {
        std::cout << " In FourierTransform: Plan creation failed!\n";
        exit(-1);
    }

} /* End of FourierTransform<long double>::FourierTransform */

FourierTransform<long double>::~FourierTransform()
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

} /* End of FourierTransform<long double>::~FourierTransform */

void FourierTransform<long double>::Forward( scalar_type* in, complex_type* out ) const
{

    /* Computes forward DFT */
    fftwl_execute_dft_r2c( plan_FFT, in, out );

} /* End of FourierTransform<long double>::Forward */

void FourierTransform<long double>::Forward( const Vector_2Dl &V, complex_type* out ) const
{

    for ( unsigned int i = 0; i < rows; i++ ) {
        for ( unsigned int j = 0; j < cols; j++ )
            in_FFT[i * cols + j] = (scalar_type) V[j][i];
    }

    /* Computes forward DFT */
    fftwl_execute_dft_r2c( plan_FFT, in_FFT, out );

} /* End of FourierTransform<long double>::Forward */

void FourierTransform<long double>::Backward( complex_type* in, scalar_type* out ) const
{

    /* Computes backward DFT */
    fftwl_execute_dft_c2r( plan_IFFT, in, out );

} /* End of FourierTransform<long double>::Backward */

void FourierTransform<long double>::Backward( complex_type* in, Vector_2Dl &V ) const
{

    /* Computes backward DFT */
    fftwl_execute_dft_c2r( plan_IFFT, in, out_IFFT );

    Scale( out_IFFT );

    for ( unsigned int i = 0; i < rows; i++ ) {
        for ( unsigned int j = 0; j < cols; j++ )
            V[j][i] = out_IFFT[i * cols + j];
    }


} /* End of FourierTransform<long double>::Backward */

void FourierTransform<long double>::Scale( scalar_type* in ) const
{

    /* Scale output */
    for ( unsigned int i = 0; i < cols * rows; i++ ) 
        in[i] /= fftScaling;

} /* End of FourierTransform<long double>::Scale */

void FourierTransform<long double>::SANDS( const Vector_2Dl &diffFactor, \
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

} /* End of FourierTransform<long double>::SANDS */




