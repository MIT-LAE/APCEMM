/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*              Spectral Advection aNd Diffusion Solver             */
/*                             (SANDS)                              */
/*                                                                  */
/* SANDS Program File                                               */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : SANDS.cpp                                 */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* 11/5/2018 : Thibaud M. Fritz
 *             Added pragma critical statement as the only thread-safe
 *             (re-entrant) routine in FFTW is fftw_execute. All the 
 *             other routines should only be called from one thread at
 *             a time. */

#include "SANDS/Solver.hpp"

void SANDS( Real_2DVector &vect, \
            Real_2DVector &diffFactor, \
            Complex_2DVector &advFactor, \
            const char* fileName_FFTW, \
            const bool realInput = 1 )
{

    unsigned int i, j;

    /* If dealing with real input, set realInput to 1.
     * Output will be the same, the algorithm is roughly
     * 3 times faster after and requires ~75% of the 
     * original memory requirement. */
        
    unsigned flags;
    /* Restoring plans from disk.
     * See following link for more information:
     * http://www.fftw.org/fftw3_doc/Wisdom.html#Wisdom */
    if ( FFTW_WISDOM && 0 ) { /* Under development ... */
        fftw_import_wisdom_from_filename( fileName_FFTW );
        flags = FFTW_WISDOM_ONLY | FFTW_PATIENT;
    }
    else
        flags = FFTW_ESTIMATE;


    /* Is data real? */
    if ( realInput ) { 
        /* After an r2c transform, the output is an ( NX * (NY/2 + 1) ) array */
        /* See the following link for more information:
         * http://www.fftw.org/fftw3_doc/Multi_002dDimensional-DFTs-of-Real-Data.html */
        const unsigned int NYH = NY/2 + 1;

        FFTW_ComplexDouble *in_IFFT, *out_FFT;
        RealDouble *in_FFT, *out_IFFT;
        in_FFT   = (RealDouble*) fftw_malloc(sizeof(RealDouble) * NX * NY);
        out_FFT  = (FFTW_ComplexDouble*) fftw_malloc(sizeof(FFTW_ComplexDouble) * NX * NYH);
    //  in_IFFT  = out_FFT .* diffFactor .* advFactor;
        in_IFFT  = (FFTW_ComplexDouble*) fftw_malloc(sizeof(FFTW_ComplexDouble) * NX * NYH);
        out_IFFT = (RealDouble*) fftw_malloc(sizeof(RealDouble) * NX * NY);

        /* Allocate FFT plans */
        fftw_plan plan_FFT, plan_IFFT;

        /* Fill input */
        for ( i = 0; i < NX; i++ ) {
            for ( j = 0; j < NY; j++ )
                in_FFT[i*NY + j] = (RealDouble) vect[j][i];
        }

        /* Define plan FFT */
        #pragma omp critical (make_plan) 
            plan_FFT = fftw_plan_dft_r2c_2d( NX, NY, in_FFT, out_FFT, flags);

        /* Execute FFT */
        fftw_execute( plan_FFT );
    
        /* Fill */
        for ( i = 0; i < NX; i++ ) {
            for ( j = 0; j < NYH; j++ ) {
                in_IFFT[i*NYH + j][0] = ( out_FFT[i*NYH + j][0] * advFactor[j][i].real() \
                                        - out_FFT[i*NYH + j][1] * advFactor[j][i].imag() ) * ( diffFactor[j][i] );
                in_IFFT[i*NYH + j][1] = ( out_FFT[i*NYH + j][0] * advFactor[j][i].imag() \
                                        + out_FFT[i*NYH + j][1] * advFactor[j][i].real() ) * ( diffFactor[j][i] );
            }
        }

        /* Define plan IFFT */ 
        #pragma omp critical (make_plan) 
            plan_IFFT = fftw_plan_dft_c2r_2d( NX, NY, in_IFFT, out_IFFT, flags);

        /* Execute IFFT */
        fftw_execute( plan_IFFT );

        /* Fill output */
        for ( i = 0; i < NX; i++ ) {
            for ( j = 0; j < NY; j++ ) 
                vect[j][i] = ( out_IFFT[i*NY+j] ) / ( RealDouble( NX * NY ) );
        }

        /* Destroy FFT plans */
        #pragma omp critical (destroy_plan)
        {
            fftw_destroy_plan( plan_FFT  );
            fftw_destroy_plan( plan_IFFT );
        }

        /* Free */
        fftw_free( in_FFT   ); 
        fftw_free( out_FFT  ); 
        fftw_free( in_IFFT  );
        fftw_free( out_IFFT );

    }
    /* Is data complex? */
    else {
    
        /* Dynamic allocation */
        FFTW_ComplexDouble *in_FFT, *out_FFT, *in_IFFT, *out_IFFT;
        in_FFT   = (FFTW_ComplexDouble*) fftw_malloc(sizeof(FFTW_ComplexDouble) * NX * NY);
        out_FFT  = (FFTW_ComplexDouble*) fftw_malloc(sizeof(FFTW_ComplexDouble) * NX * NY);
    //  in_IFFT  = out_FFT .* diffFactor .* advFactor;
        in_IFFT  = (FFTW_ComplexDouble*) fftw_malloc(sizeof(FFTW_ComplexDouble) * NX * NY);
        out_IFFT = (FFTW_ComplexDouble*) fftw_malloc(sizeof(FFTW_ComplexDouble) * NX * NY);
    
        /* Allocate FFT plans */
        fftw_plan plan_FFT, plan_IFFT;

        /* Fill input */
        for ( i = 0; i < NX; i++ ) {
            for ( j = 0; j < NY; j++ )
                in_FFT[i*NY + j][0] = (RealDouble) vect[j][i];
        }

        /* Define plan FFT */
        #pragma omp critical (make_plan) 
            plan_FFT  = fftw_plan_dft_2d( NX, NY, in_FFT, out_FFT, FFTW_FORWARD, flags);

        /* Execute FFT */
        fftw_execute( plan_FFT );
    
        /* Fill */
        for ( i = 0; i < NX; i++ ) {
            for ( j = 0; j < NY; j++ ) {
              in_IFFT[i*NY + j][0] = ( out_FFT[i*NY + j][0] * advFactor[j][i].real() \
                                     - out_FFT[i*NY + j][1] * advFactor[j][i].imag() ) * ( diffFactor[j][i] );
              in_IFFT[i*NY + j][1] = ( out_FFT[i*NY + j][0] * advFactor[j][i].imag() \
                                     + out_FFT[i*NY + j][1] * advFactor[j][i].real() ) * ( diffFactor[j][i] );
            }
         }

        /* Define plan IFFT */ 
        #pragma omp critical (make_plan) 
            plan_IFFT = fftw_plan_dft_2d( NX, NY, in_IFFT, out_IFFT, FFTW_BACKWARD, flags);

        /* Execute IFFT */
        fftw_execute( plan_IFFT );

        /* Fill output */
        for ( i = 0; i < NX; i++ ) {
            for ( j = 0; j < NY; j++ ) 
                vect[j][i] = ( out_IFFT[i*NY+j][0] ) / ( RealDouble( NX * NY ) );
        }

        /* Destroy FFT plans */
        #pragma omp critical (destroy_plan)
        {
            fftw_destroy_plan( plan_FFT  );
            fftw_destroy_plan( plan_IFFT );
        }

        /* Free */
        fftw_free( in_FFT   ); 
        fftw_free( out_FFT  ); 
        fftw_free( in_IFFT  );
        fftw_free( out_IFFT );
    
    }

/*    1D EXAMPLE    */

//    fftw_complex *in_d, *out_d;
//    
//    in_d = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
//    out_d = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
//    for ( unsigned int i = 0; i < N; i++ ) {
//        in_d[i][0] = (double) i;
//        std::cout << "in_d[" << i << "] = " << *in_d[i] << std::endl;
//    }
//    
//    fftw_plan p;
//    
//    p = fftw_plan_dft_1d( N, &in_d[0], &out_d[0], FFTW_FORWARD, FFTW_ESTIMATE );
//    fftw_execute( p );
//
//    fftw_destroy_plan( p );
//
//    std::cout << "in = [ ";
//    for ( int i = 0; i < N; i++ )
//        std::cout << *in_d[i] << ", " ;
//    std::cout << "]" << std::endl;
//
//    std::cout << "out = [ ";
//    for ( int i = 0; i < N; i++ )
//        std::cout << *out_d[i] << ", " ;
//    std::cout << "]" << std::endl;
//
//    fftw_free(in_d);
//    fftw_free(out_d);

   
/*    2D EXAMPLE    */

//    fftw_complex *in, *out;

//    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
//    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
//    for ( unsigned int i = 0; i < N; i++ ) {
//        for ( unsigned int j = 0; j < N; j++ ) {
//            in[i*N + j][0] = (double) i+j;
//            in[i*N + j][1] = (double) 0.0;
//        }
//    }
//
//    fftw_plan p;
//
//    p = fftw_plan_dft_2d( N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
//    fftw_execute( p );
//
//    fftw_destroy_plan( p );

//    std::cout << "in = \n";
//    for ( unsigned int i = 0; i < N; i++ ) {
//        std::cout << "[ ";
//        for ( unsigned int j = 0; j < N; j++ )
//            std::cout << in[i*N+j][0] << ", " ;
//        std::cout << "]" << std::endl;
//    }
//
//    std::cout << "out = \n";
//    for ( unsigned int i = 0; i < N; i++ ) {
//        std::cout << "[ ";
//        for ( unsigned int j = 0; j < N; j++ )
//            std::cout << out[i*N+j][0] << ", " ;
//        std::cout << "]" << std::endl;
//    }
//
//    fftw_free(in);
//    fftw_free(out);

}
