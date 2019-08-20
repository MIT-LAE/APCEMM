/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*              Spectral Advection aNd Diffusion Solver             */
/*                             (SANDS)                              */
/*                                                                  */
/* FFT Header File                                                  */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 11/14/2018                                */
/* File                 : FFT.hpp                                   */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef FFT_H_INCLUDED
#define FFT_H_INCLUDED

#include <iostream>
#include <cstdlib>
#ifdef OMP
    #include "omp.h"
#endif /* OMP */

/* Include Parameters.hpp for multithreading option */
#include <Core/Parameters.hpp>

#include <Util/ForwardDecl.hpp>
#include <fftw3.h>

#define REAL      0
#define IMAG      1


template<typename T>
inline T* fftw_cast( const T* p )
{
    return const_cast<T*>( p);
}

inline fftw_complex* fftw_cast( const std::complex<double> *p )
{
    return const_cast<fftw_complex*>( reinterpret_cast<const fftw_complex*>(p) );
}

inline fftwf_complex* fftw_cast( const std::complex<float> *p )
{
    return const_cast<fftwf_complex*>( reinterpret_cast<const fftwf_complex*>(p) );
}

inline fftwl_complex* fftw_cast( const std::complex<long double> *p )
{
    return const_cast<fftwl_complex*>( reinterpret_cast<const fftwl_complex*>(p) );
}

/* Fourier Transform Interface */

template <typename T>
class FourierTransform_1D {};
template <typename T>
class FourierTransform_2D {};

template <>
class FourierTransform_1D<float>
{

    /* Typedefs */
    typedef float scalar_type;
    typedef fftwf_complex complex_type;

    public:

        /** 
         * Constructor
         *
         * @param MULTITHREAD_FFT  : Use threaded FFT?
         * @param WISDOM (bool)    : Use FFTW_WISDOM?
         * @param FFTW_DIR (char*) : Path to storage for FFTW plans?
         * @param rows_ (UInt)     : number of rows 
         */

        FourierTransform_1D( const bool MULTITHREADED_FFT, \
                             const bool WISDOM,            \
                             const char* FFTW_DIR,         \
                             const UInt rows_ );

        /**
         * Destructor 
         *
         * Destroy the FFT and IFFT plans 
         */

        ~FourierTransform_1D();

        /**
         * Performs forward 1D Fourier transform of real-valued data
         * Takes scalar array as an input and outputs a complex array
         *
         * Reuses the FFT plan.
         *
         * @param in  (1D scalar)  : input to the FFT
         * @param out (1D complex) : output to the FFT
         */

//        virtual void Forward( scalar_type* in, complex_type* out ) const;
//
//        /**
//         * Performs forward 1D Fourier transform of real-valued data
//         * Takes scalar 1D vector as an input and outputs a complex array
//         *
//         * Reuses the FFT plan.
//         *
//         * @param V   (1D scalar)  : input to the FFT
//         * @param out (1D complex) : output to the FFT
//         */
//
//        virtual void Forward( const Vector_1Df &V, complex_type* out ) const;
//
//        /**
//         * Performs backward 1D Fourier transform of real-valued data
//         * Takes complex array as an input and outputs a real-valued array
//         *
//         * Reuses the IFFT plan.
//         *
//         * Results are scaled by rows
//         *
//         * @param in  (1D complex) : input to the IFFT
//         * @param out (1D scalar)  : output to the IFFT
//         */
//
//        virtual void Backward( complex_type* in, scalar_type* out ) const;
//
//        /**
//         * Performs backward 1D Fourier transform of real-valued data
//         * Takes complex array as an input and outputs a 1D real-valued vector
//         *
//         * Reuses the IFFT plan.
//         *
//         * Results are scaled by rows
//         *
//         * @param in (1D complex) : input to the IFFT
//         * @param V  (1D scalar)  : output to the IFFT
//         */
//
//        virtual void Backward( complex_type* in, Vector_1Df &V ) const;
//
//        /** 
//         * Scales the input array by 1 / rows 
//         *
//         * @param in (1D scalar) : array to be scaled
//         */
//
//        virtual void Scale( scalar_type* in ) const;
//
//        /**
//         * Performs backward 1D Fourier transform of real-valued data
//         * Takes complex array as an input and outputs a 1D real-valued vector
//         *
//         * Reuses the IFFT plan.
//         *
//         * Results are scaled by rows
//         *
//         * @param in (2D complex) : Advection array corresponding to shear
//         * @param V  (2D scalar)  : vector to be "sheared"
//         */

        virtual void ApplyShear( const Vector_2Dcf &shearFactor, \
                                 Vector_2Df &V ) const;

        /* Rows of real and complex data */
        const UInt rows;
        /* Rows in the complex spectrum */
        const UInt rowsC;

        /* Scaling factor */
        const UInt fftScaling;

    private:

        /* Switch for threaded FFT? */
        const bool THREADED_FFT;

        /* Pointer to arrays for fftw_plans */
        scalar_type  *in_FFT , *out_IFFT;
        complex_type *in_IFFT, *out_FFT;

        /* Reusable plan for forward transformation */
        fftwf_plan plan_FFT;
        /* Reusable plan for backward transformation */
        fftwf_plan plan_IFFT;

};

template <>
class FourierTransform_1D<double>
{

    /* Typedefs */
    typedef double scalar_type;
    typedef fftw_complex complex_type;

    public:

        /** 
         * Constructor
         *
         * @param MULTITHREAD_FFT  : Use threaded FFT?
         * @param WISDOM (bool)    : Use FFTW_WISDOM?
         * @param FFTW_DIR (char*) : Path to storage for FFTW plans?
         * @param rows_ (UInt)     : number of rows 
         */

        FourierTransform_1D( const bool MULTITHREADED_FFT, \
                             const bool WISDOM,            \
                             const char* FFTW_DIR,         \
                             const UInt rows_ );

        /**
         * Destructor 
         *
         * Destroy the FFT and IFFT plans 
         */

        ~FourierTransform_1D();

        /**
         * Performs forward 1D Fourier transform of real-valued data
         * Takes scalar array as an input and outputs a complex array
         *
         * Reuses the FFT plan.
         *
         * @param in  (1D scalar)  : input to the FFT
         * @param out (1D complex) : output to the FFT
         */

//        virtual void Forward( scalar_type* in, complex_type* out ) const;
//
//        /**
//         * Performs forward 1D Fourier transform of real-valued data
//         * Takes scalar 1D vector as an input and outputs a complex array
//         *
//         * Reuses the FFT plan.
//         *
//         * @param V   (1D scalar)  : input to the FFT
//         * @param out (1D complex) : output to the FFT
//         */
//
//        virtual void Forward( const Vector_1D &V, complex_type* out ) const;
//
//        /**
//         * Performs backward 1D Fourier transform of real-valued data
//         * Takes complex array as an input and outputs a real-valued array
//         *
//         * Reuses the IFFT plan.
//         *
//         * Results are scaled by rows
//         *
//         * @param in  (1D complex) : input to the IFFT
//         * @param out (1D scalar)  : output to the IFFT
//         */
//
//        virtual void Backward( complex_type* in, scalar_type* out ) const;
//
//        /**
//         * Performs backward 1D Fourier transform of real-valued data
//         * Takes complex array as an input and outputs a 1D real-valued vector
//         *
//         * Reuses the IFFT plan.
//         *
//         * Results are scaled by rows
//         *
//         * @param in (1D complex) : input to the IFFT
//         * @param V  (1D scalar)  : output to the IFFT
//         */
//
//        virtual void Backward( complex_type* in, Vector_1D &V ) const;
//
//        /** 
//         * Scales the input array by 1 / rows 
//         *
//         * @param in (1D scalar) : array to be scaled
//         */
//
//        virtual void Scale( scalar_type* in ) const;
//
//        /**
//         * Performs backward 1D Fourier transform of real-valued data
//         * Takes complex array as an input and outputs a 1D real-valued vector
//         *
//         * Reuses the IFFT plan.
//         *
//         * Results are scaled by rows
//         *
//         * @param in (2D complex) : Advection array corresponding to shear
//         * @param V  (2D scalar)  : vector to be "sheared"
//         */

        virtual void ApplyShear( const Vector_2Dc &shearFactor, \
                                 Vector_2D &V ) const;

        /* Rows of real and complex data */
        const UInt rows;
        /* Rows in the complex spectrum */
        const UInt rowsC;

        /* Scaling factor */
        const UInt fftScaling;

    private:

        /* Switch for threaded FFT? */
        const bool THREADED_FFT;

        /* Pointer to arrays for fftw_plans */
        scalar_type  *in_FFT , *out_IFFT;
        complex_type *in_IFFT, *out_FFT;

        /* Reusable plan for forward transformation */
        fftw_plan plan_FFT;
        /* Reusable plan for backward transformation */
        fftw_plan plan_IFFT;

};

template <>
class FourierTransform_1D<long double>
{

    /* Typedefs */
    typedef long double scalar_type;
    typedef fftwl_complex complex_type;

    public:

        /** 
         * Constructor
         *
         * @param MULTITHREAD_FFT  : Use threaded FFT?
         * @param WISDOM (bool)    : Use FFTW_WISDOM?
         * @param FFTW_DIR (char*) : Path to storage for FFTW plans?
         * @param rows_ (UInt)     : number of rows 
         */

        FourierTransform_1D( const bool MULTITHREADED_FFT, \
                             const bool WISDOM,            \
                             const char* FFTW_DIR,         \
                             const UInt rows_ );

        /**
         * Destructor 
         *
         * Destroy the FFT and IFFT plans 
         */

        ~FourierTransform_1D();

        /**
         * Performs forward 1D Fourier transform of real-valued data
         * Takes scalar array as an input and outputs a complex array
         *
         * Reuses the FFT plan.
         *
         * @param in  (1D scalar)  : input to the FFT
         * @param out (1D complex) : output to the FFT
         */

//        virtual void Forward( scalar_type* in, complex_type* out ) const;
//
//        /**
//         * Performs forward 1D Fourier transform of real-valued data
//         * Takes scalar 1D vector as an input and outputs a complex array
//         *
//         * Reuses the FFT plan.
//         *
//         * @param V   (1D scalar)  : input to the FFT
//         * @param out (1D complex) : output to the FFT
//         */
//
//        virtual void Forward( const Vector_1Dl &V, complex_type* out ) const;
//
//        /**
//         * Performs backward 1D Fourier transform of real-valued data
//         * Takes complex array as an input and outputs a real-valued array
//         *
//         * Reuses the IFFT plan.
//         *
//         * Results are scaled by rows
//         *
//         * @param in  (1D complex) : input to the IFFT
//         * @param out (1D scalar)  : output to the IFFT
//         */
//
//        virtual void Backward( complex_type* in, scalar_type* out ) const;
//
//        /**
//         * Performs backward 1D Fourier transform of real-valued data
//         * Takes complex array as an input and outputs a 1D real-valued vector
//         *
//         * Reuses the IFFT plan.
//         *
//         * Results are scaled by rows
//         *
//         * @param in (1D complex) : input to the IFFT
//         * @param V  (1D scalar)  : output to the IFFT
//         */
//
//        virtual void Backward( complex_type* in, Vector_1Dl &Vl ) const;
//
//        /** 
//         * Scales the input array by 1 / rows 
//         *
//         * @param in (1D scalar) : array to be scaled
//         */
//
//        virtual void Scale( scalar_type* in ) const;
//
//        /**
//         * Performs backward 1D Fourier transform of real-valued data
//         * Takes complex array as an input and outputs a 1D real-valued vector
//         *
//         * Reuses the IFFT plan.
//         *
//         * Results are scaled by rows
//         *
//         * @param in (2D complex) : Advection array corresponding to shear
//         * @param V  (2D scalar)  : vector to be "sheared"
//         */

        virtual void ApplyShear( const Vector_2Dcl &shearFactor, \
                                 Vector_2Dl &V ) const;

        /* Rows of real and complex data */
        const UInt rows;
        /* Rows in the complex spectrum */
        const UInt rowsC;

        /* Scaling factor */
        const UInt fftScaling;

    private:

        /* Switch for threaded FFT? */
        const bool THREADED_FFT;

        /* Pointer to arrays for fftw_plans */
        scalar_type  *in_FFT , *out_IFFT;
        complex_type *in_IFFT, *out_FFT;

        /* Reusable plan for forward transformation */
        fftwl_plan plan_FFT;
        /* Reusable plan for backward transformation */
        fftwl_plan plan_IFFT;

};

template <>
class FourierTransform_2D<float>
{

    /* Typedefs */
    typedef float scalar_type;
    typedef fftwf_complex complex_type;

    public:

        /** 
         * Constructor
         *
         * @param MULTITHREAD_FFT  : Use threaded FFT?
         * @param WISDOM (bool)    : Use FFTW_WISDOM?
         * @param FFTW_DIR (char*) : Path to storage for FFTW plans?
         * @param rows_ (UInt)     : number of rows 
         * @param cols_ (UInt)     : number of columns 
         */

        FourierTransform_2D( const bool MULTITHREADED_FFT, \
                             const bool WISDOM,            \
                             const char* FFTW_DIR,         \
                             const UInt rows_,             \
                             const UInt cols_ );

        /**
         * Destructor 
         *
         * Destroy the FFT and IFFT plans 
         */

        ~FourierTransform_2D();

        /**
         * Performs forward 2D Fourier transform of real-valued data
         * Takes scalar array as an input and outputs a complex array
         *
         * Reuses the FFT plan.
         *
         * @param in  (1D scalar)  : input to the FFT
         * @param out (1D complex) : output to the FFT
         */

//        virtual void Forward( scalar_type* in, complex_type* out ) const;
//
//        /**
//         * Performs forward 2D Fourier transform of real-valued data
//         * Takes scalar 2D vector as an input and outputs a complex array
//         *
//         * Reuses the FFT plan.
//         *
//         * @param V   (2D scalar)  : input to the FFT
//         * @param out (1D complex) : output to the FFT
//         */
//
//        virtual void Forward( const Vector_2Df &V, complex_type* out ) const;
//
//        /**
//         * Performs backward 2D Fourier transform of real-valued data
//         * Takes complex array as an input and outputs a real-valued array
//         *
//         * Reuses the IFFT plan.
//         *
//         * Results are scaled by rows * cols
//         *
//         * @param in  (1D complex) : input to the IFFT
//         * @param out (1D scalar)  : output to the IFFT
//         */
//
//        virtual void Backward( complex_type* in, scalar_type* out ) const;
//
//        /**
//         * Performs backward 2D Fourier transform of real-valued data
//         * Takes complex array as an input and outputs a 2D real-valued vector
//         *
//         * Reuses the IFFT plan.
//         *
//         * Results are scaled by rows * cols
//         *
//         * @param in (1D complex) : input to the IFFT
//         * @param V  (2D scalar)  : output to the IFFT
//         */
//
//        virtual void Backward( complex_type* in, Vector_2Df &V ) const;
//
//        /** 
//         * Scales the input array by 1 / (rows * cols) 
//         *
//         * @param in (1D scalar) : array to be scaled
//         */
//
//        virtual void Scale( scalar_type* in ) const;
//
//        /** 
//         * Solves the 2D diffusion-advection equation
//         *
//         * @param DiffFactor (2D scalar)  : diffusion factor for the 2D spectral solver 
//         * @param AdvFactor  (2D complex) : advection factor for the 2D spectral solver 
//         * @param V          (2D scalar)  : vector to be transported
//         */

        void SANDS( const Vector_2Df &DiffFactor, \
                    const Vector_2Dcf &AdvFactor, \
                    Vector_2Df &V ) const;

        /* Rows of real and complex data */
        const UInt rows;
        /* Columns of real data only */
        const UInt cols;
        /* Columns in the complex spectrum */
        const UInt colsC;

        /* Scaling factor */
        const UInt fftScaling;

    private:

        /* Switch for threaded FFT? */
        const bool THREADED_FFT;

        /* Pointer to arrays for fftw_plans */
        scalar_type  *in_FFT , *out_IFFT;
        complex_type *in_IFFT, *out_FFT;

        /* Reusable plan for forward transformation */
        fftwf_plan plan_FFT;
        /* Reusable plan for backward transformation */
        fftwf_plan plan_IFFT;

};

template <>
class FourierTransform_2D<double>
{

    /* Typedefs */
    typedef double scalar_type;
    typedef fftw_complex complex_type;

    public:

        /** 
         * Constructor
         *
         * @param MULTITHREAD_FFT  : Use threaded FFT?
         * @param WISDOM (bool)    : Use FFTW_WISDOM?
         * @param FFTW_DIR (char*) : Path to storage for FFTW plans?
         * @param rows_ (UInt)     : number of rows 
         * @param cols_ (UInt)     : number of columns 
         */

        FourierTransform_2D( const bool MULTITHREADED_FFT, \
                             const bool WISDOM,            \
                             const char* FFTW_DIR,         \
                             const UInt rows_,             \
                             const UInt cols_ );

        /**
         * Destructor 
         *
         * Destroy the FFT and IFFT plans 
         */

        ~FourierTransform_2D();

        /**
         * Performs forward 2D Fourier transform of real-valued data
         * Takes scalar array as an input and outputs a complex array
         *
         * Reuses the FFT plan.
         *
         * @param in  (1D scalar)  : input to the FFT
         * @param out (1D complex) : output to the FFT
         */

//        virtual void Forward( scalar_type* in, complex_type* out ) const;
//
//        /**
//         * Performs forward 2D Fourier transform of real-valued data
//         * Takes scalar 2D vector as an input and outputs a complex array
//         *
//         * Reuses the FFT plan.
//         *
//         * @param V   (2D scalar)  : input to the FFT
//         * @param out (1D complex) : output to the FFT
//         */
//
//        virtual void Forward( const Vector_2D &V, complex_type* out ) const;
//
//        /**
//         * Performs backward 2D Fourier transform of real-valued data
//         * Takes complex array as an input and outputs a real-valued array
//         *
//         * Reuses the IFFT plan.
//         *
//         * Results are scaled by rows * cols
//         *
//         * @param in  (1D complex) : input to the IFFT
//         * @param out (1D scalar)  : output to the IFFT
//         */
//
//        virtual void Backward( complex_type* in, scalar_type* out ) const;
//
//        /**
//         * Performs backward 2D Fourier transform of real-valued data
//         * Takes complex array as an input and outputs a 2D real-valued vector
//         *
//         * Reuses the IFFT plan.
//         *
//         * Results are scaled by rows * cols
//         *
//         * @param in (1D complex) : input to the IFFT
//         * @param V  (2D scalar)  : output to the IFFT
//         */
//
//        virtual void Backward( complex_type* in, Vector_2D &V ) const;
//
//        /** 
//         * Scales the input array by 1 / (rows * cols) 
//         *
//         * @param in (1D scalar) : array to be scaled
//         */
//
//        virtual void Scale( scalar_type* in ) const;
//
//        /** 
//         * Solves the 2D diffusion-advection equation
//         *
//         * @param DiffFactor (2D scalar)  : diffusion factor for the 2D spectral solver 
//         * @param AdvFactor  (2D complex) : advection factor for the 2D spectral solver 
//         * @param V          (2D scalar)  : vector to be transported
//         */

        void SANDS( const Vector_2D &DiffFactor, \
                    const Vector_2Dc &AdvFactor, \
                    Vector_2D &V ) const;

        /* Rows of real and complex data */
        const UInt rows;
        /* Columns of real data only */
        const UInt cols;
        /* Columns in the complex spectrum */
        const UInt colsC;

        /* Scaling factor */
        const UInt fftScaling;

    private:

        /* Switch for threaded FFT? */
        const bool THREADED_FFT;

        /* Pointer to arrays for fftw_plans */
        scalar_type  *in_FFT , *out_IFFT;
        complex_type *in_IFFT, *out_FFT;

        /* Reusable plan for forward transformation */
        fftw_plan plan_FFT;
        /* Reusable plan for backward transformation */
        fftw_plan plan_IFFT;


};

template <>
class FourierTransform_2D<long double>
{

    /* Typedefs */
    typedef long double scalar_type;
    typedef fftwl_complex complex_type;

    public:

        /** 
         * Constructor
         *
         * @param MULTITHREAD_FFT  : Use threaded FFT?
         * @param WISDOM (bool)    : Use FFTW_WISDOM?
         * @param FFTW_DIR (char*) : Path to storage for FFTW plans?
         * @param rows_ (UInt)     : number of rows 
         * @param cols_ (UInt)     : number of columns 
         */

        FourierTransform_2D( const bool MULTITHREADED_FFT, \
                             const bool WISDOM,            \
                             const char* FFTW_DIR,         \
                             const UInt rows_,             \
                             const UInt cols_ );

        /**
         * Destructor 
         *
         * Destroy the FFT and IFFT plans 
         */

        ~FourierTransform_2D();

        /**
         * Performs forward 2D Fourier transform of real-valued data
         * Takes scalar array as an input and outputs a complex array
         *
         * Reuses the FFT plan.
         *
         * @param in  (1D scalar)  : input to the FFT
         * @param out (1D complex) : output to the FFT
         */

//        virtual void Forward( scalar_type* in, complex_type* out ) const;
//
//        /**
//         * Performs forward 2D Fourier transform of real-valued data
//         * Takes scalar 2D vector as an input and outputs a complex array
//         *
//         * Reuses the FFT plan.
//         *
//         * @param V   (2D scalar)  : input to the FFT
//         * @param out (1D complex) : output to the FFT
//         */
//
//        virtual void Forward( const Vector_2Dl &V, complex_type* out ) const;
//
//        /**
//         * Performs backward 2D Fourier transform of real-valued data
//         * Takes complex array as an input and outputs a real-valued array
//         *
//         * Reuses the IFFT plan.
//         *
//         * Results are scaled by rows * cols
//         *
//         * @param in  (1D complex) : input to the IFFT
//         * @param out (1D scalar)  : output to the IFFT
//         */
//
//        virtual void Backward( complex_type* in, scalar_type* out ) const;
//
//        /**
//         * Performs backward 2D Fourier transform of real-valued data
//         * Takes complex array as an input and outputs a 2D real-valued vector
//         *
//         * Reuses the IFFT plan.
//         *
//         * Results are scaled by rows * cols
//         *
//         * @param in (1D complex) : input to the IFFT
//         * @param V  (2D scalar)  : output to the IFFT
//         */
//
//        virtual void Backward( complex_type* in, Vector_2Dl &V ) const;
//
//        /** 
//         * Scales the input array by 1 / (rows * cols) 
//         *
//         * @param in (1D scalar) : array to be scaled
//         */
//
//        virtual void Scale( scalar_type* in ) const;
//
//        /** 
//         * Solves the 2D diffusion-advection equation
//         *
//         * @param DiffFactor (2D scalar)  : diffusion factor for the 2D spectral solver 
//         * @param AdvFactor  (2D complex) : advection factor for the 2D spectral solver 
//         * @param V          (2D scalar)  : vector to be transported
//         */

        void SANDS( const Vector_2Dl &DiffFactor, \
                    const Vector_2Dcl &AdvFactor, \
                    Vector_2Dl &V ) const;

        /* Rows of real and complex data */
        const UInt rows;
        /* Columns of real data only */
        const UInt cols;
        /* Columns in the complex spectrum */
        const UInt colsC;

        /* Scaling factor */
        const UInt fftScaling;

    private:

        /* Switch for threaded FFT? */
        const bool THREADED_FFT;

        /* Pointer to arrays for fftw_plans */
        scalar_type  *in_FFT , *out_IFFT;
        complex_type *in_IFFT, *out_FFT;

        /* Reusable plan for forward transformation */
        fftwl_plan plan_FFT;
        /* Reusable plan for backward transformation */
        fftwl_plan plan_IFFT;


};


#endif /* FFT_H_INCLUDED */
