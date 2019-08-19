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

FourierTransform_1D<float>::FourierTransform_1D( const bool WISDOM,    \
                                                 const char* FFTW_DIR, \
                                                 const UInt rows_ )
    :   rows( rows_ ),
        rowsC( rows_/2 + 1 ),
        fftScaling( rows_ )
{

    UInt nThreads = 0;
    int wisdomExists = 0;

    std::string fileName = FFTW_DIR;
    std::string fileName_FFT, fileName_IFFT;

    fileName += "FFTW_1Dplan_" + std::to_string(rows);

    if ( !PARALLEL_CASES ) {
        /* Performing the one-time initialization required to use threads with
         * FFTW */
        int Success = fftw_init_threads();
        if ( Success == 0 ) {
            std::cout << " Could not perform the initialization required by FFTW when using multiple threads" << std::endl;
            exit(-1);
        }
    }

    /* Allocate the ins and outs */
    in_FFT   = (scalar_type*)  fftwf_malloc( sizeof(scalar_type)  * rows  );
    out_FFT  = (complex_type*) fftwf_malloc( sizeof(complex_type) * rowsC );
    in_IFFT  = (complex_type*) fftwf_malloc( sizeof(complex_type) * rowsC );
    out_IFFT = (scalar_type*)  fftwf_malloc( sizeof(scalar_type)  * rows  );

    if ( !PARALLEL_CASES ) {
        int nThreads = omp_get_max_threads();

        /* All plans subsequently created with any planner routine will use 
         * that many threads */
        fftw_plan_with_nthreads( nThreads );

        fileName += "_" + std::to_string(nThreads);
    }

    fileName_FFT  = fileName + "_FFT.pl";
    fileName_IFFT = fileName + "_IFFT.pl";

    if ( WISDOM )
        wisdomExists = fftw_import_wisdom_from_filename( fileName_FFT.c_str() );

    /* Create FFT plan and compute the best FFT algorithm */
    if ( wisdomExists )
        plan_FFT = fftwf_plan_dft_r2c_1d( rows, in_FFT, out_FFT, FFTW_WISDOM_ONLY );
    else
        plan_FFT = fftwf_plan_dft_r2c_1d( rows, in_FFT, out_FFT, FFTW_PATIENT );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_FFT == NULL ) {
        std::cout << " In FourierTransform_1D: Plan creation failed!\n";
        exit(-1);
    }

    if ( WISDOM && ( wisdomExists == 0 ) ) {
        /* Export wisdom from FFTW plan */
        fftw_export_wisdom_to_filename( fileName_FFT.c_str() );
    }

    if ( WISDOM )
        wisdomExists = fftw_import_wisdom_from_filename( fileName_IFFT.c_str() );

    /* Create IFFT plan and compute the best IFFT algorithm */
    if ( wisdomExists )
        plan_IFFT = fftwf_plan_dft_c2r_1d( rows, in_IFFT, out_IFFT, FFTW_WISDOM_ONLY );
    else
        plan_IFFT = fftwf_plan_dft_c2r_1d( rows, in_IFFT, out_IFFT, FFTW_PATIENT );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_IFFT == NULL ) {
        std::cout << " In FourierTransform_1D: Plan creation failed!\n";
        exit(-1);
    }

    if ( WISDOM && ( wisdomExists == 0 ) ) {
        /* Export wisdom from FFTW plan */
        fftw_export_wisdom_to_filename( fileName_IFFT.c_str() );
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

    if ( !PARALLEL_CASES ) {
        /* Cleanup and get rid of all memory allocated by FFTW */
        fftw_cleanup_threads();
    }

} /* End of FourierTransform_1D<float>::~FourierTransform_1D */

//void FourierTransform_1D<float>::Forward( scalar_type* in, complex_type* out ) const
//{
//
//    /* Computes forward DFT */
//    fftwf_execute_dft_r2c( plan_FFT, in, out );
//
//} /* End of FourierTransform_1D<float>::Forward */
//
//void FourierTransform_1D<float>::Forward( const Vector_1Df &V, complex_type* out ) const
//{
//
//    UInt i = 0;
//
//#pragma omp parallel for   \
//    private ( i          ) \
//    schedule( dynamic, 1 )
//    for ( i = 0; i < rows; i++ )
//        in_FFT[i] = (scalar_type) V[i];
//
//    /* Computes forward DFT */
//    fftwf_execute_dft_r2c( plan_FFT, in_FFT, out );
//
//} /* End of FourierTransform_1D<float>::Forward */
//
//void FourierTransform_1D<float>::Backward( complex_type* in, scalar_type* out ) const
//{
//
//    /* Computes backward DFT */
//    fftwf_execute_dft_c2r( plan_IFFT, in, out );
//
//} /* End of FourierTransform_1D<float>::Backward */
//
//void FourierTransform_1D<float>::Backward( complex_type* in, Vector_1Df &V ) const
//{
//
//    UInt i = 0;
//
//    /* Computes backward DFT */
//    fftwf_execute_dft_c2r( plan_IFFT, in, out_IFFT );
//
//    Scale( out_IFFT );
//
//#pragma omp parallel for   \
//    private ( i          ) \
//    schedule( dynamic, 1 )
//    for ( i = 0; i < rows; i++ )
//        V[i] = out_IFFT[i];
//
//
//} /* End of FourierTransform_1D<float>::Backward */
//
//void FourierTransform_1D<float>::Scale( scalar_type* in ) const
//{
//
//    UInt i = 0;
//
//    /* Scale output */
//
//#pragma omp parallel for   \
//    private ( i          ) \
//    schedule( dynamic, 1 )
//    for ( i = 0; i < rows; i++ ) 
//        in[i] /= fftScaling;
//
//} /* End of FourierTransform_1D<float>::Scale */

void FourierTransform_1D<float>::ApplyShear( const Vector_2Dcf &shearFactor, \
                                             Vector_2Df &V ) const
{

    UInt i = 0;
    UInt j = 0;

    /* Dealing with each layer one at a time */
    for ( j = 0; j < V.size(); j++ ) {

        /* Storing horizontal layer */
#pragma omp parallel for       \
        private ( i          ) \
        schedule( dynamic, 1 )
        for ( i = 0; i < V[0].size(); i++ )
            in_FFT[i] = (scalar_type) V[j][i];

        /* Computes forward DFT */
        fftwf_execute_dft_r2c( plan_FFT, in_FFT, out_FFT );

        /* Convolve and scale the frequencies */
#pragma omp parallel for       \
        private ( i          ) \
        schedule( dynamic, 1 )
        for ( i = 0; i < rowsC; i++ ) {
            in_IFFT[i][REAL] = ( out_FFT[i][REAL] * shearFactor[j][i].real() \
                               - out_FFT[i][IMAG] * shearFactor[j][i].imag() );
            in_IFFT[i][IMAG] = ( out_FFT[i][REAL] * shearFactor[j][i].imag() \
                               + out_FFT[i][IMAG] * shearFactor[j][i].real() );
        }

        /* Computes backward DFT */
        fftwf_execute_dft_c2r( plan_IFFT, in_IFFT, out_IFFT );

        /* Apply results back into original 2D vector */
#pragma omp parallel for       \
        private ( i          ) \
        schedule( dynamic, 1 )
        for ( i = 0; i < V[0].size(); i++ )
            V[j][i] = out_IFFT[i] / fftScaling;

    }

} /* End of FourierTransform_1D<float>::ApplyShear */

FourierTransform_1D<double>::FourierTransform_1D( const bool WISDOM,    \
                                                  const char* FFTW_DIR, \
                                                  const UInt rows_ )
    :   rows( rows_ ),
        rowsC( rows_/2 + 1 ),
        fftScaling( rows_ )
{

    UInt nThreads = 0;
    int wisdomExists = 0;

    std::string fileName = FFTW_DIR;
    std::string fileName_FFT, fileName_IFFT; 

    fileName += "FFTW_1Dplan_" + std::to_string(rows);

    if ( !PARALLEL_CASES ) {
        /* Performing the one-time initialization required to use threads with
         * FFTW */
        int Success = fftw_init_threads();
        if ( Success == 0 ) {
            std::cout << " Could not perform the initialization required by FFTW when using multiple threads" << std::endl;
            exit(-1);
        }
    }

    /* Allocate the ins and outs */
    in_FFT   = (scalar_type*)  fftw_malloc( sizeof(scalar_type)  * rows  );
    out_FFT  = (complex_type*) fftw_malloc( sizeof(complex_type) * rowsC );
    in_IFFT  = (complex_type*) fftw_malloc( sizeof(complex_type) * rowsC );
    out_IFFT = (scalar_type*)  fftw_malloc( sizeof(scalar_type)  * rows  );

    if ( !PARALLEL_CASES ) {
        int nThreads = omp_get_max_threads();

        /* All plans subsequently created with any planner routine will use 
         * that many threads */
        fftw_plan_with_nthreads( nThreads );

        fileName += "_" + std::to_string(nThreads);
    }

    fileName_FFT  = fileName + "_FFT.pl";
    fileName_IFFT = fileName + "_IFFT.pl";

    if ( WISDOM )
        wisdomExists = fftw_import_wisdom_from_filename( fileName_FFT.c_str() );

    /* Create FFT plan and compute the best FFT algorithm */
    if ( wisdomExists )
        plan_FFT = fftw_plan_dft_r2c_1d( rows, in_FFT, out_FFT, FFTW_WISDOM_ONLY );
    else
        plan_FFT = fftw_plan_dft_r2c_1d( rows, in_FFT, out_FFT, FFTW_PATIENT );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_FFT == NULL ) {
        std::cout << " In FourierTransform_1D: Plan creation failed!\n";
        exit(-1);
    }

    if ( WISDOM && ( wisdomExists == 0 ) ) {
        /* Export wisdom from FFTW plan */
        fftw_export_wisdom_to_filename( fileName_FFT.c_str() );
    }

    if ( WISDOM )
        wisdomExists = fftw_import_wisdom_from_filename( fileName_IFFT.c_str() );

    /* Create IFFT plan and compute the best IFFT algorithm */
    if ( wisdomExists )
        plan_IFFT = fftw_plan_dft_c2r_1d( rows, in_IFFT, out_IFFT, FFTW_WISDOM_ONLY );
    else
        plan_IFFT = fftw_plan_dft_c2r_1d( rows, in_IFFT, out_IFFT, FFTW_PATIENT );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_IFFT == NULL ) {
        std::cout << " In FourierTransform_1D: Plan creation failed!\n";
        exit(-1);
    }

    if ( WISDOM && ( wisdomExists == 0 ) ) {
        /* Export wisdom from FFTW plan */
        fftw_export_wisdom_to_filename( fileName_IFFT.c_str() );
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

    if ( !PARALLEL_CASES ) {
        /* Cleanup and get rid of all memory allocated by FFTW */
        fftw_cleanup_threads();
    }

} /* End of FourierTransform_1D<double>::~FourierTransform_1D */

//void FourierTransform_1D<double>::Forward( scalar_type* in, complex_type* out ) const
//{
//
//        /* Computes forward DFT */
//        fftw_execute_dft_r2c( plan_FFT, in, out );
//
//} /* End of FourierTransform_1D<double>::Forward */
//
//void FourierTransform_1D<double>::Forward( const Vector_1D &V, complex_type* out ) const
//{
//
//    UInt i = 0;
//
//#pragma omp parallel for   \
//    private ( i          ) \
//    schedule( dynamic, 1 )
//    for ( i = 0; i < rows; i++ )
//        in_FFT[i] = (scalar_type) V[i];
//
//    /* Computes forward DFT */
//    fftw_execute_dft_r2c( plan_FFT, in_FFT, out );
//
//} /* End of FourierTransform_1D<double>::Forward */
//
//void FourierTransform_1D<double>::Backward( complex_type* in, scalar_type* out ) const
//{
//
//    /* Computes backward DFT */
//    fftw_execute_dft_c2r( plan_IFFT, in, out );
//
//} /* End of FourierTransform_1D<double>::Backward */
//
//void FourierTransform_1D<double>::Backward( complex_type* in, Vector_1D &V ) const
//{
//
//    UInt i = 0;
//
//    /* Computes backward DFT */
//    fftw_execute_dft_c2r( plan_IFFT, in, out_IFFT );
//
//    Scale( out_IFFT );
//
//#pragma omp parallel for   \
//    private ( i          ) \
//    schedule( dynamic, 1 )
//    for ( i = 0; i < rows; i++ )
//        V[i] = out_IFFT[i];
//
//
//} /* End of FourierTransform_1D<double>::Backward */
//
//void FourierTransform_1D<double>::Scale( scalar_type* in ) const
//{
//
//    UInt i = 0;
//
//    /* Scale output */
//
//#pragma omp parallel for   \
//    private ( i          ) \
//    schedule( dynamic, 1 )
//    for ( i = 0; i < rows; i++ ) 
//        in[i] /= fftScaling;
//
//} /* End of FourierTransform_1D<double>::Scale */

void FourierTransform_1D<double>::ApplyShear( const Vector_2Dc &shearFactor, \
                                              Vector_2D &V ) const
{

    UInt i = 0;
    UInt j = 0;

    /* Dealing with each layer one at a time */
    for ( j = 0; j < V.size(); j++ ) {

        /* Storing horizontal layer */
#pragma omp parallel for       \
        private ( i          ) \
        schedule( dynamic, 1 )
        for ( i = 0; i < V[0].size(); i++ )
            in_FFT[i] = (scalar_type) V[j][i];

        /* Computes forward DFT */
        fftw_execute_dft_r2c( plan_FFT, in_FFT, out_FFT );

        /* Convolve and scale the frequencies */
#pragma omp parallel for       \
        private ( i          ) \
        schedule( dynamic, 1 )
        for ( i = 0; i < rowsC; i++ ) {
            in_IFFT[i][REAL] = ( out_FFT[i][REAL] * shearFactor[j][i].real() \
                               - out_FFT[i][IMAG] * shearFactor[j][i].imag() );
            in_IFFT[i][IMAG] = ( out_FFT[i][REAL] * shearFactor[j][i].imag() \
                               + out_FFT[i][IMAG] * shearFactor[j][i].real() );
        }

        /* Computes backward DFT */
        fftw_execute_dft_c2r( plan_IFFT, in_IFFT, out_IFFT );

        /* Apply results back into original 2D vector */
#pragma omp parallel for       \
        private ( i          ) \
        schedule( dynamic, 1 )
        for ( i = 0; i < V[0].size(); i++ )
            V[j][i] = out_IFFT[i] / fftScaling;

    }

} /* End of FourierTransform_1D<double>::ApplyShear */

FourierTransform_1D<long double>::FourierTransform_1D( const bool WISDOM,    \
                                                       const char* FFTW_DIR, \
                                                       const UInt rows_ )
    :   rows( rows_ ),
        rowsC( rows_/2 + 1 ),
        fftScaling( rows_ )
{

    UInt nThreads = 0;
    int wisdomExists = 0;

    std::string fileName = FFTW_DIR;
    std::string fileName_FFT, fileName_IFFT; 

    fileName += "FFTW_1Dplan_" + std::to_string(rows);

    if ( !PARALLEL_CASES ) {
        /* Performing the one-time initialization required to use threads with
         * FFTW */
        int Success = fftw_init_threads();
        if ( Success == 0 ) {
            std::cout << " Could not perform the initialization required by FFTW when using multiple threads" << std::endl;
            exit(-1);
        }
    }

    /* Allocate the ins and outs */
    in_FFT   = (scalar_type*)  fftwl_malloc( sizeof(scalar_type)  * rows  );
    out_FFT  = (complex_type*) fftwl_malloc( sizeof(complex_type) * rowsC );
    in_IFFT  = (complex_type*) fftwl_malloc( sizeof(complex_type) * rowsC );
    out_IFFT = (scalar_type*)  fftwl_malloc( sizeof(scalar_type)  * rows  );

    if ( !PARALLEL_CASES ) {
        int nThreads = omp_get_max_threads();

        /* All plans subsequently created with any planner routine will use 
         * that many threads */
        fftw_plan_with_nthreads( nThreads );

        fileName += "_" + std::to_string(nThreads);
    }

    fileName_FFT  = fileName + "_FFT.pl";
    fileName_IFFT = fileName + "_IFFT.pl";

    if ( WISDOM )
        wisdomExists = fftw_import_wisdom_from_filename( fileName_FFT.c_str() );

    /* Create FFT plan and compute the best FFT algorithm */
    if ( wisdomExists )
        plan_FFT = fftwl_plan_dft_r2c_1d( rows, in_FFT, out_FFT, FFTW_WISDOM_ONLY );
    else
        plan_FFT = fftwl_plan_dft_r2c_1d( rows, in_FFT, out_FFT, FFTW_PATIENT );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_FFT == NULL ) {
        std::cout << " In FourierTransform_1D: Plan creation failed!\n";
        exit(-1);
    }

    if ( WISDOM && ( wisdomExists == 0 ) ) {
        /* Export wisdom from FFTW plan */
        fftw_export_wisdom_to_filename( fileName_FFT.c_str() );
    }

    if ( WISDOM )
        wisdomExists = fftw_import_wisdom_from_filename( fileName_IFFT.c_str() );

    /* Create IFFT plan and compute the best IFFT algorithm */
    if ( wisdomExists )
        plan_IFFT = fftwl_plan_dft_c2r_1d( rows, in_IFFT, out_IFFT, FFTW_PATIENT );
    else
        plan_IFFT = fftwl_plan_dft_c2r_1d( rows, in_IFFT, out_IFFT, FFTW_WISDOM_ONLY );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_IFFT == NULL ) {
        std::cout << " In FourierTransform_1D: Plan creation failed!\n";
        exit(-1);
    }

    if ( WISDOM && ( wisdomExists == 0 ) ) {
        /* Export wisdom from FFTW plan */
        fftw_export_wisdom_to_filename( fileName_IFFT.c_str() );
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

    if ( !PARALLEL_CASES ) {
        /* Cleanup and get rid of all memory allocated by FFTW */
        fftw_cleanup_threads();
    }

} /* End of FourierTransform_1D<long double>::~FourierTransform_1D */

//void FourierTransform_1D<long double>::Forward( scalar_type* in, complex_type* out ) const
//{
//
//    /* Computes forward DFT */
//    fftwl_execute_dft_r2c( plan_FFT, in, out );
//
//} /* End of FourierTransform_1D<long double>::Forward */
//
//void FourierTransform_1D<long double>::Forward( const Vector_1Dl &V, complex_type* out ) const
//{
//
//    UInt i = 0;
//
//#pragma omp parallel for   \
//    private ( i          ) \
//    schedule( dynamic, 1 )
//    for ( i = 0; i < rows; i++ )
//        in_FFT[i] = (scalar_type) V[i];
//
//    /* Computes forward DFT */
//    fftwl_execute_dft_r2c( plan_FFT, in_FFT, out );
//
//} /* End of FourierTransform_1D<long double>::Forward */
//
//void FourierTransform_1D<long double>::Backward( complex_type* in, scalar_type* out ) const
//{
//
//    /* Computes backward DFT */
//    fftwl_execute_dft_c2r( plan_IFFT, in, out );
//
//} /* End of FourierTransform_1D<long double>::Backward */
//
//void FourierTransform_1D<long double>::Backward( complex_type* in, Vector_1Dl &V ) const
//{
//
//    UInt i = 0;
//
//    /* Computes backward DFT */
//    fftwl_execute_dft_c2r( plan_IFFT, in, out_IFFT );
//
//    Scale( out_IFFT );
//
//#pragma omp parallel for   \
//    private ( i          ) \
//    schedule( dynamic, 1 )
//    for ( i = 0; i < rows; i++ )
//        V[i] = out_IFFT[i];
//
//
//} /* End of FourierTransform_1D<long double>::Backward */
//
//void FourierTransform_1D<long double>::Scale( scalar_type* in ) const
//{
//
//    UInt i = 0;
//
//    /* Scale output */
//
//#pragma omp parallel for   \
//    private ( i          ) \
//    schedule( dynamic, 1 )
//    for ( i = 0; i < rows; i++ ) 
//        in[i] /= fftScaling;
//
//} /* End of FourierTransform_1D<long double>::Scale */

void FourierTransform_1D<long double>::ApplyShear( const Vector_2Dcl &shearFactor, \
                                                   Vector_2Dl &V ) const
{

    UInt i = 0;
    UInt j = 0;

    /* Dealing with each layer one at a time */
    for ( j = 0; j < V.size(); j++ ) {

        /* Storing horizontal layer */
#pragma omp parallel for       \
        private ( i          ) \
        schedule( dynamic, 1 )
        for ( i = 0; i < V[0].size(); i++ )
            in_FFT[i] = (scalar_type) V[j][i];

        /* Computes forward DFT */
        fftwl_execute_dft_r2c( plan_FFT, in_FFT, out_FFT );

        /* Convolve and scale the frequencies */
#pragma omp parallel for       \
        private ( i          ) \
        schedule( dynamic, 1 )
        for ( i = 0; i < rowsC; i++ ) {
            in_IFFT[i][REAL] = ( out_FFT[i][REAL] * shearFactor[j][i].real() \
                               - out_FFT[i][IMAG] * shearFactor[j][i].imag() );
            in_IFFT[i][IMAG] = ( out_FFT[i][REAL] * shearFactor[j][i].imag() \
                               + out_FFT[i][IMAG] * shearFactor[j][i].real() );
        }

        /* Computes backward DFT */
        fftwl_execute_dft_c2r( plan_IFFT, in_IFFT, out_IFFT );

        /* Apply results back into original 2D vector */
#pragma omp parallel for       \
        private ( i          ) \
        schedule( dynamic, 1 )
        for ( i = 0; i < V[0].size(); i++ )
            V[j][i] = out_IFFT[i] / fftScaling;

    }

} /* End of FourierTransform_1D<long double>::ApplyShear */

FourierTransform_2D<float>::FourierTransform_2D( const bool WISDOM,    \
                                                 const char* FFTW_DIR, \
                                                 const UInt rows_,     \
                                                 const UInt cols_)
    :   rows( rows_ ),
        cols( cols_ ),
        colsC( cols_/2 + 1 ),
        fftScaling( cols_ * rows_ )
{

    UInt nThreads = 0;
    int wisdomExists = 0;

    std::string fileName = FFTW_DIR;
    std::string fileName_FFT, fileName_IFFT; 

    fileName += "FFTW_2Dplan_" + std::to_string(rows) \
                             + "_" + std::to_string(cols);

    if ( !PARALLEL_CASES ) {
        /* Performing the one-time initialization required to use threads with
         * FFTW */
        int Success = fftw_init_threads();
        if ( Success == 0 ) {
            std::cout << " Could not perform the initialization required by FFTW when using multiple threads" << std::endl;
            exit(-1);
        }
    }

    /* Allocate the ins and outs */
    in_FFT   = (scalar_type*)  fftwf_malloc( sizeof(scalar_type)  * rows * cols  );
    out_FFT  = (complex_type*) fftwf_malloc( sizeof(complex_type) * rows * colsC );
    in_IFFT  = (complex_type*) fftwf_malloc( sizeof(complex_type) * rows * colsC );
    out_IFFT = (scalar_type*)  fftwf_malloc( sizeof(scalar_type)  * rows * cols  );

    if ( !PARALLEL_CASES ) {
        int nThreads = omp_get_max_threads();

        /* All plans subsequently created with any planner routine will use 
         * that many threads */
        fftw_plan_with_nthreads( nThreads );

        fileName += "_" + std::to_string(nThreads);
    }

    fileName_FFT  = fileName + "_FFT.pl";
    fileName_IFFT = fileName + "_IFFT.pl";

    if ( WISDOM )
        wisdomExists = fftw_import_wisdom_from_filename( fileName_FFT.c_str() );

    /* Create FFT plan and compute the best FFT algorithm */
    if ( wisdomExists )
        plan_FFT = fftwf_plan_dft_r2c_2d( rows, cols, in_FFT, out_FFT, FFTW_WISDOM_ONLY );
    else
        plan_FFT = fftwf_plan_dft_r2c_2d( rows, cols, in_FFT, out_FFT, FFTW_PATIENT );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_FFT == NULL ) {
        std::cout << " In FourierTransform_2D: Plan creation failed!\n";
        exit(-1);
    }

    if ( WISDOM && ( wisdomExists == 0 ) ) {
        /* Export wisdom from FFTW plan */
        fftw_export_wisdom_to_filename( fileName_FFT.c_str() );
    }

    if ( WISDOM )
        wisdomExists = fftw_import_wisdom_from_filename( fileName_IFFT.c_str() );

    /* Create IFFT plan and compute the best IFFT algorithm */
    if ( wisdomExists )
        plan_IFFT = fftwf_plan_dft_c2r_2d( rows, cols, in_IFFT, out_IFFT, FFTW_WISDOM_ONLY );
    else
        plan_IFFT = fftwf_plan_dft_c2r_2d( rows, cols, in_IFFT, out_IFFT, FFTW_PATIENT );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_IFFT == NULL ) {
        std::cout << " In FourierTransform_2D: Plan creation failed!\n";
        exit(-1);
    }

    if ( WISDOM && ( wisdomExists == 0 ) ) {
        /* Export wisdom from FFTW plan */
        fftw_export_wisdom_to_filename( fileName_IFFT.c_str() );
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

    if ( !PARALLEL_CASES ) {
        /* Cleanup and get rid of all memory allocated by FFTW */
        fftw_cleanup_threads();
    }

} /* End of FourierTransform_2D<float>::~FourierTransform_2D */

//void FourierTransform_2D<float>::Forward( scalar_type* in, complex_type* out ) const
//{
//
//    /* Computes forward DFT */
//    fftwf_execute_dft_r2c( plan_FFT, in, out );
//
//} /* End of FourierTransform_2D<float>::Forward */
//
//void FourierTransform_2D<float>::Forward( const Vector_2Df &V, complex_type* out ) const
//{
//
//    UInt i = 0;
//    UInt j = 0;
//
//#pragma omp parallel for   \
//    private ( i, j       ) \
//    schedule( dynamic, 1 )
//    for ( i = 0; i < rows; i++ ) {
//        for ( j = 0; j < cols; j++ )
//            in_FFT[i * cols + j] = (scalar_type) V[j][i];
//    }
//
//    /* Computes forward DFT */
//    fftwf_execute_dft_r2c( plan_FFT, in_FFT, out );
//
//} /* End of FourierTransform_2D<float>::Forward */
//
//void FourierTransform_2D<float>::Backward( complex_type* in, scalar_type* out ) const
//{
//
//    /* Computes backward DFT */
//    fftwf_execute_dft_c2r( plan_IFFT, in, out );
//
//} /* End of FourierTransform_2D<float>::Backward */
//
//void FourierTransform_2D<float>::Backward( complex_type* in, Vector_2Df &V ) const
//{
//
//    UInt i = 0;
//    UInt j = 0;
//
//    /* Computes backward DFT */
//    fftwf_execute_dft_c2r( plan_IFFT, in, out_IFFT );
//
//    Scale( out_IFFT );
//
//#pragma omp parallel for   \
//    private ( i, j       ) \
//    schedule( dynamic, 1 )
//    for ( i = 0; i < rows; i++ ) {
//        for ( j = 0; j < cols; j++ )
//            V[j][i] = out_IFFT[i * cols + j];
//    }
//
//
//} /* End of FourierTransform_2D<float>::Backward */
//
//void FourierTransform_2D<float>::Scale( scalar_type* in ) const
//{
//
//    UInt i = 0;
//
//    /* Scale output */
//
//#pragma omp parallel for   \
//    private ( i          ) \
//    schedule( dynamic, 1 )
//    for ( i = 0; i < cols * rows; i++ ) 
//        in[i] /= fftScaling;
//
//} /* End of FourierTransform_2D<float>::Scale */

void FourierTransform_2D<float>::SANDS( const Vector_2Df &diffFactor, \
                                        const Vector_2Dcf &advFactor, \
                                        Vector_2Df &V ) const
{

    UInt i = 0;
    UInt j = 0;

#pragma omp parallel for   \
    private ( i, j       ) \
    schedule( dynamic, 1 )
    for ( i = 0; i < rows; i++ ) {
        for ( j = 0; j < cols; j++ )
            in_FFT[i * cols + j] = (scalar_type) V[j][i];
    }

    /* Computes forward DFT */
    fftwf_execute_dft_r2c( plan_FFT, in_FFT, out_FFT );

    /* Convolve and scale the frequencies */
#pragma omp parallel for   \
    private ( i, j       ) \
    schedule( dynamic, 1 )
    for ( i = 0; i < rows; i++ ) {
        for ( j = 0; j < colsC; j++ ) {
            in_IFFT[i * colsC + j][REAL] = ( out_FFT[i * colsC + j][REAL] * advFactor[j][i].real() \
                                           - out_FFT[i * colsC + j][IMAG] * advFactor[j][i].imag() ) * diffFactor[j][i] ;
            in_IFFT[i * colsC + j][IMAG] = ( out_FFT[i * colsC + j][REAL] * advFactor[j][i].imag() \
                                           + out_FFT[i * colsC + j][IMAG] * advFactor[j][i].real() ) * diffFactor[j][i] ;
        }
    }

    /* Computes backward DFT */
    fftwf_execute_dft_c2r( plan_IFFT, in_IFFT, out_IFFT );

#pragma omp parallel for   \
    private ( i, j       ) \
    schedule( dynamic, 1 )
    for ( i = 0; i < rows; i++ ) {
        for ( j = 0; j < cols; j++ ) {
            V[j][i] = out_IFFT[i * cols + j] / fftScaling;
        }
    }

} /* End of FourierTransform_2D<float>::SANDS */

FourierTransform_2D<double>::FourierTransform_2D( const bool WISDOM,    \
                                                  const char* FFTW_DIR, \
                                                  const UInt rows_,     \
                                                  const UInt cols_ )
    :   rows( rows_ ),
        cols( cols_ ),
        colsC( cols_/2 + 1 ),
        fftScaling( cols_ * rows_ )
{

    UInt nThreads = 0;
    int wisdomExists = 0;

    std::string fileName = FFTW_DIR;
    std::string fileName_FFT, fileName_IFFT; 

    fileName += "FFTW_2Dplan_" + std::to_string(rows) \
                             + "_" + std::to_string(cols);

    if ( !PARALLEL_CASES ) {
        /* Performing the one-time initialization required to use threads with
         * FFTW */
        int Success = fftw_init_threads();
        if ( Success == 0 ) {
            std::cout << " Could not perform the initialization required by FFTW when using multiple threads" << std::endl;
            exit(-1);
        }
    }

    /* Allocate the ins and outs */
    in_FFT   = (scalar_type*)  fftw_malloc( sizeof(scalar_type)  * rows * cols  );
    out_FFT  = (complex_type*) fftw_malloc( sizeof(complex_type) * rows * colsC );
    in_IFFT  = (complex_type*) fftw_malloc( sizeof(complex_type) * rows * colsC );
    out_IFFT = (scalar_type*)  fftw_malloc( sizeof(scalar_type)  * rows * cols  );

    if ( !PARALLEL_CASES ) {
        nThreads = omp_get_max_threads();

        /* All plans subsequently created with any planner routine will use 
         * that many threads */
        fftw_plan_with_nthreads( nThreads );

        fileName += "_" + std::to_string(nThreads);
    }

    fileName_FFT  = fileName + "_FFT.pl";
    fileName_IFFT = fileName + "_IFFT.pl";

    if ( WISDOM )
        wisdomExists = fftw_import_wisdom_from_filename( fileName_FFT.c_str() );

    /* Create FFT plan and compute the best FFT algorithm */
    if ( wisdomExists )
        plan_FFT = fftw_plan_dft_r2c_2d( rows, cols, in_FFT, out_FFT, FFTW_WISDOM_ONLY );
    else
        plan_FFT = fftw_plan_dft_r2c_2d( rows, cols, in_FFT, out_FFT, FFTW_PATIENT );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_FFT == NULL ) {
        std::cout << " In FourierTransform_2D: Plan creation failed!\n";
        exit(-1);
    }

    if ( WISDOM && ( wisdomExists == 0 ) ) {
        /* Export wisdom from FFTW plan */
        fftw_export_wisdom_to_filename( fileName_FFT.c_str() );
    }

    if ( WISDOM )
        wisdomExists = fftw_import_wisdom_from_filename( fileName_IFFT.c_str() );

    /* Create IFFT plan and compute the best IFFT algorithm */
    if ( wisdomExists )
        plan_IFFT = fftw_plan_dft_c2r_2d( rows, cols, in_IFFT, out_IFFT, FFTW_WISDOM_ONLY );
    else
        plan_IFFT = fftw_plan_dft_c2r_2d( rows, cols, in_IFFT, out_IFFT, FFTW_PATIENT );

    /* Check that plan was successfully created. 
     * Otherwise exit. */
    if ( plan_IFFT == NULL ) {
        std::cout << " In FourierTransform_2D: Plan creation failed!\n";
        exit(-1);
    }

    if ( WISDOM && ( wisdomExists == 0 ) ) {
        /* Export wisdom from FFTW plan */
        fftw_export_wisdom_to_filename( fileName_IFFT.c_str() );
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

    if ( !PARALLEL_CASES ) {
        /* Cleanup and get rid of all memory allocated by FFTW */
        fftw_cleanup_threads();
    }

} /* End of FourierTransform_2D<double>::~FourierTransform_2D */

//void FourierTransform_2D<double>::Forward( scalar_type* in, complex_type* out ) const
//{
//
//    /* Computes forward DFT */
//    fftw_execute_dft_r2c( plan_FFT, in, out );
//
//} /* End of FourierTransform_2D<double>::Forward */
//
//void FourierTransform_2D<double>::Forward( const Vector_2D &V, complex_type* out ) const
//{
//
//    UInt i = 0;
//    UInt j = 0;
//
//#pragma omp parallel for   \
//    private ( i, j       ) \
//    schedule( dynamic, 1 )
//    for ( i = 0; i < rows; i++ ) {
//        for ( j = 0; j < cols; j++ )
//            in_FFT[i * cols + j] = (scalar_type) V[j][i];
//    }
//
//    /* Computes forward DFT */
//    fftw_execute_dft_r2c( plan_FFT, in_FFT, out );
//
//} /* End of FourierTransform_2D<double>::Forward */
//
//void FourierTransform_2D<double>::Backward( complex_type* in, scalar_type* out ) const
//{
//
//    /* Computes backward DFT */
//    fftw_execute_dft_c2r( plan_IFFT, in, out );
//
//} /* End of FourierTransform_2D<double>::Backward */
//
//void FourierTransform_2D<double>::Backward( complex_type* in, Vector_2D &V ) const
//{
//
//    UInt i = 0;
//    UInt j = 0;
//
//    /* Computes backward DFT */
//    fftw_execute_dft_c2r( plan_IFFT, in, out_IFFT );
//
//    Scale( out_IFFT );
//
//#pragma omp parallel for   \
//    private ( i, j       ) \
//    schedule( dynamic, 1 )
//    for ( i = 0; i < rows; i++ ) {
//        for ( j = 0; j < cols; j++ )
//            V[j][i] = out_IFFT[i * cols + j];
//    }
//
//
//} /* End of FourierTransform_2D<double>::Backward */
//
//void FourierTransform_2D<double>::Scale( scalar_type* in ) const
//{
//
//    UInt i = 0;
//
//    /* Scale output */
//
//#pragma omp parallel for   \
//    private ( i          ) \
//    schedule( dynamic, 1 )
//    for ( i = 0; i < cols * rows; i++ ) 
//        in[i] /= fftScaling;
//
//} /* End of FourierTransform_2D<double>::Scale */

void FourierTransform_2D<double>::SANDS( const Vector_2D &diffFactor, \
                                         const Vector_2Dc &advFactor, \
                                         Vector_2D &V ) const
{

    UInt i = 0;
    UInt j = 0;

#pragma omp parallel for   \
    private ( i, j       ) \
    schedule( dynamic, 1 )
    for ( i = 0; i < rows; i++ ) {
        for ( j = 0; j < cols; j++ )
            in_FFT[i * cols + j] = (scalar_type) V[j][i];
    }

    /* Computes forward DFT */
    fftw_execute_dft_r2c( plan_FFT, in_FFT, out_FFT );

    /* Convolve and scale the frequencies */
#pragma omp parallel for   \
    private ( i, j       ) \
    schedule( dynamic, 1 )
    for ( i = 0; i < rows; i++ ) {
        for ( j = 0; j < colsC; j++ ) {
            in_IFFT[i * colsC + j][REAL] = ( out_FFT[i * colsC + j][REAL] * advFactor[j][i].real() \
                                           - out_FFT[i * colsC + j][IMAG] * advFactor[j][i].imag() ) * diffFactor[j][i] ;
            in_IFFT[i * colsC + j][IMAG] = ( out_FFT[i * colsC + j][REAL] * advFactor[j][i].imag() \
                                           + out_FFT[i * colsC + j][IMAG] * advFactor[j][i].real() ) * diffFactor[j][i] ;
        }
    }

    /* Computes backward DFT */
    fftw_execute_dft_c2r( plan_IFFT, in_IFFT, out_IFFT );

#pragma omp parallel for   \
    private ( i, j       ) \
    schedule( dynamic, 1 )
    for ( i = 0; i < rows; i++ ) {
        for ( j = 0; j < cols; j++ ) {
            V[j][i] = out_IFFT[i * cols + j] / fftScaling;
        }
    }

} /* End of FourierTransform_2D<double>::SANDS */


FourierTransform_2D<long double>::FourierTransform_2D( const bool WISDOM,    \
                                                       const char* FFTW_DIR, \
                                                       const UInt rows_,     \
                                                       const UInt cols_ )
    :   rows( rows_ ),
        cols( cols_ ),
        colsC( cols_/2 + 1 ),
        fftScaling( cols_ * rows_ )
{

    UInt nThreads = 0;
    int wisdomExists = 0;

    std::string fileName = FFTW_DIR;
    std::string fileName_FFT, fileName_IFFT; 

    fileName += "FFTW_2Dplan_" + std::to_string(rows) \
                             + "_" + std::to_string(cols);

    if ( !PARALLEL_CASES ) {
        /* Performing the one-time initialization required to use threads with
         * FFTW */
        int Success = fftw_init_threads();
        if ( Success == 0 ) {
            std::cout << " Could not perform the initialization required by FFTW when using multiple threads" << std::endl;
            exit(-1);
        }
    }

    /* Allocate the ins and outs */
    in_FFT   = (scalar_type*)  fftwl_malloc( sizeof(scalar_type)  * rows * cols  );
    out_FFT  = (complex_type*) fftwl_malloc( sizeof(complex_type) * rows * colsC );
    in_IFFT  = (complex_type*) fftwl_malloc( sizeof(complex_type) * rows * colsC );
    out_IFFT = (scalar_type*)  fftwl_malloc( sizeof(scalar_type)  * rows * cols  );

    if ( !PARALLEL_CASES ) {
        int nThreads = omp_get_max_threads();

        /* All plans subsequently created with any planner routine will use 
         * that many threads */
        fftw_plan_with_nthreads( nThreads );

        fileName += "_" + std::to_string(nThreads);
    }

    fileName_FFT  = fileName + "_FFT.pl";
    fileName_IFFT = fileName + "_IFFT.pl";

    if ( WISDOM )
        wisdomExists = fftw_import_wisdom_from_filename( fileName_FFT.c_str() );

    /* Create FFT plan and compute the best FFT algorithm */
    if ( wisdomExists )
        plan_FFT = fftwl_plan_dft_r2c_2d( rows, cols, in_FFT, out_FFT, FFTW_WISDOM_ONLY );
    else
        plan_FFT = fftwl_plan_dft_r2c_2d( rows, cols, in_FFT, out_FFT, FFTW_PATIENT );

    /* Check that plan was successfully created.
     * Otherwise exit. */
    if ( plan_FFT == NULL ) {
        std::cout << " In FourierTransform_2D: Plan creation failed!\n";
        exit(-1);
    }

    if ( WISDOM && ( wisdomExists == 0 ) ) {
        /* Export wisdom from FFTW plan */
        fftw_export_wisdom_to_filename( fileName_FFT.c_str() );
    }

    if ( WISDOM )
        wisdomExists = fftw_import_wisdom_from_filename( fileName_IFFT.c_str() );

    /* Create IFFT plan and compute the best IFFT algorithm */
    if ( wisdomExists )
        plan_IFFT = fftwl_plan_dft_c2r_2d( rows, cols, in_IFFT, out_IFFT, FFTW_WISDOM_ONLY );
    else
        plan_IFFT = fftwl_plan_dft_c2r_2d( rows, cols, in_IFFT, out_IFFT, FFTW_PATIENT );

    /* Check that plan was successfully created.
     * Otherwise exit. */
    if ( plan_IFFT == NULL ) {
        std::cout << " In FourierTransform_2D: Plan creation failed!\n";
        exit(-1);
    }

    if ( WISDOM && ( wisdomExists == 0 ) ) {
        /* Export wisdom from FFTW plan */
        fftw_export_wisdom_to_filename( fileName_IFFT.c_str() );
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

    if ( !PARALLEL_CASES ) {
        /* Cleanup and get rid of all memory allocated by FFTW */
        fftw_cleanup_threads();
    }

} /* End of FourierTransform_2D<long double>::~FourierTransform_2D */

//void FourierTransform_2D<long double>::Forward( scalar_type* in, complex_type* out ) const
//{
//
//    /* Computes forward DFT */
//    fftwl_execute_dft_r2c( plan_FFT, in, out );
//
//} /* End of FourierTransform_2D<long double>::Forward */
//
//void FourierTransform_2D<long double>::Forward( const Vector_2Dl &V, complex_type* out ) const
//{
//
//    UInt i = 0;
//    UInt j = 0;
//
//#pragma omp parallel for   \
//    private ( i, j       ) \
//    schedule( dynamic, 1 )
//    for ( i = 0; i < rows; i++ ) {
//        for ( j = 0; j < cols; j++ )
//            in_FFT[i * cols + j] = (scalar_type) V[j][i];
//    }
//
//    /* Computes forward DFT */
//    fftwl_execute_dft_r2c( plan_FFT, in_FFT, out );
//
//} /* End of FourierTransform_2D<long double>::Forward */
//
//void FourierTransform_2D<long double>::Backward( complex_type* in, scalar_type* out ) const
//{
//
//    /* Computes backward DFT */
//    fftwl_execute_dft_c2r( plan_IFFT, in, out );
//
//} /* End of FourierTransform_2D<long double>::Backward */
//
//void FourierTransform_2D<long double>::Backward( complex_type* in, Vector_2Dl &V ) const
//{
//
//    UInt i = 0;
//    UInt j = 0;
//
//    /* Computes backward DFT */
//    fftwl_execute_dft_c2r( plan_IFFT, in, out_IFFT );
//
//    Scale( out_IFFT );
//
//#pragma omp parallel for   \
//    private ( i, j       ) \
//    schedule( dynamic, 1 )
//    for ( i = 0; i < rows; i++ ) {
//        for ( j = 0; j < cols; j++ )
//            V[j][i] = out_IFFT[i * cols + j];
//    }
//
//
//} /* End of FourierTransform_2D<long double>::Backward */
//
//void FourierTransform_2D<long double>::Scale( scalar_type* in ) const
//{
//
//    UInt i = 0;
//
//    /* Scale output */
//
//#pragma omp parallel for   \
//    private ( i          ) \
//    schedule( dynamic, 1 )
//    for ( i = 0; i < cols * rows; i++ ) 
//        in[i] /= fftScaling;
//
//} /* End of FourierTransform_2D<long double>::Scale */

void FourierTransform_2D<long double>::SANDS( const Vector_2Dl &diffFactor, \
                                              const Vector_2Dcl &advFactor, \
                                              Vector_2Dl &V ) const
{

    UInt i = 0;
    UInt j = 0;

#pragma omp parallel for   \
    private ( i, j       ) \
    schedule( dynamic, 1 )
    for ( i = 0; i < rows; i++ ) {
        for ( j = 0; j < cols; j++ )
            in_FFT[i * cols + j] = (scalar_type) V[j][i];
    }

    /* Computes forward DFT */
    fftwl_execute_dft_r2c( plan_FFT, in_FFT, out_FFT );

    /* Convolve and scale the frequencies */
#pragma omp parallel for   \
    private ( i, j       ) \
    schedule( dynamic, 1 )
    for ( i = 0; i < rows; i++ ) {
        for ( j = 0; j < colsC; j++ ) {
            in_IFFT[i * colsC + j][REAL] = ( out_FFT[i * colsC + j][REAL] * advFactor[j][i].real() \
                                           - out_FFT[i * colsC + j][IMAG] * advFactor[j][i].imag() ) * diffFactor[j][i] ;
            in_IFFT[i * colsC + j][IMAG] = ( out_FFT[i * colsC + j][REAL] * advFactor[j][i].imag() \
                                           + out_FFT[i * colsC + j][IMAG] * advFactor[j][i].real() ) * diffFactor[j][i] ;
        }
    }

    /* Computes backward DFT */
    fftwl_execute_dft_c2r( plan_IFFT, in_IFFT, out_IFFT );

#pragma omp parallel for   \
    private ( i, j       ) \
    schedule( dynamic, 1 )
    for ( i = 0; i < rows; i++ ) {
        for ( j = 0; j < cols; j++ ) {
            V[j][i] = out_IFFT[i * cols + j] / fftScaling;
        }
    }

} /* End of FourierTransform_2D<long double>::SANDS */



