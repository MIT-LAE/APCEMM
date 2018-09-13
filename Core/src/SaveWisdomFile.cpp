/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* SaveWisdomFile Program File                                      */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : SaveWisdomFile.cpp                        */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <fftw3.h>

#include "Parameters.hpp"

typedef fftw_complex T;

void SaveWisdomFile( std::vector<std::vector<double> > &Data, const char* fileName )
{

    unsigned int i, j;
    int success;

    /* Dynamic allocation */
    T *in_FFT, *out_FFT;
    in_FFT   = (T*) fftw_malloc(sizeof(T) * NCELL);
    out_FFT  = (T*) fftw_malloc(sizeof(T) * NCELL);

    unsigned flags = FFTW_EXHAUSTIVE; 
//    unsigned flags = FFTW_MEASURE;
    
    /* Allocate FFT plan */
    fftw_plan plan_FFT;

    plan_FFT  = fftw_plan_dft_2d( NX, NY, in_FFT, out_FFT, FFTW_FORWARD, flags );

    /* Fill */
    for ( i = 0; i < NX; i++ ) {
        for ( j = 0; j < NY; j++ ) {
            in_FFT[i*NY + j][0] = (double) Data[j][i];
            in_FFT[i*NY + j][1] = (double) 0.0;
        }
    }

    fftw_execute( plan_FFT );

    success = fftw_export_wisdom_to_filename( fileName );

    if ( success == 0 )
        std::cout << "In SaveWisdomFile.cpp: Couldn't save FFTW3's wisdom to disk (fileName: " << fileName << ")" << std::endl;

    fftw_destroy_plan( plan_FFT );
}
