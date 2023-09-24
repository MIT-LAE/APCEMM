/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*              Spectral Advection aNd Diffusion Solver             */
/*                             (SANDS)                              */
/*                                                                  */
/* Solver Header File                                               */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 8/12/2018                                 */
/* File                 : Solver.hpp                                */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef SOLVER_H_INCLUDED
#define SOLVER_H_INCLUDED

#include <iostream>
#include <vector>
#include <complex>
#include <fftw3.h>
#include <fstream>
#ifdef OMP
    #include "omp.h"
#endif /* OMP */

#include "Core/Parameters.hpp"
#include "Core/Interface.hpp"
#include "Util/PhysConstant.hpp"
#include "Util/ForwardDecl.hpp"
#include "SANDS/FFT.hpp"

namespace SANDS 
{

    /* A class for solving the 2D-Avection Diffusion equation */

    class Solver
    {

        public:
            
            /**
             * Constructor
             */

            Solver();

            Solver(int nx, int ny,
                   double xlim_left, double xlim_right,
                   double ylim_up, double ylim_down,
                   bool doFill = true,
                   int fillOpt = 1,
                   double fillVal = 0);

            /** 
             * Initializer 
             *
             * Initialize 2D solver
             *
             * @param MULTITHREADED_FFT: Use threaded FFT?
             * @param WISDOM (bool)    : Use FFTW_WISDOM? (default = 0)
             * @param FFTW_DIR (char*) : Path to storage for FFTW (default = "")
             * @param fill (bool)      : Fill negative values? (default = 1)
             * @param fillVal (double) : Fill with? (default = 0.0E+00) 
             * @param fillOpt (int)    : Fill option
             */

            void Initialize( const bool MULTITHREADED_FFT,         \
                             const bool WISDOM = 0,                \
                             const char* FFTW_DIR = "",            \
                             const bool fill = 1,                  \
                             const double fillValue = 0.0E+00, \
                             const UInt fillOpt_ = 1 );

            /**
             * Destructor 
             *
             * Destroys 2D solver
             */

            ~Solver( );

            /**
             * Assign frequencies 
             */

            void AssignFreq( );

            /**
             * Update time step 
             *
             * @param T (double) : New timestep [s]
             *
             * Returns an error if T <= 0.0
             */

            void UpdateTimeStep( const double T );

            /**
             * Update diffusion field 
             *
             * @param dH (double) : Horizontal diffusion coefficient [m2/s]
             * @param dV (double) : Vertical diffusion coefficient [m2/s]
             *
             * Returns an error if dH or dV is non positive 
             */

            void UpdateDiff( const double dH, const double dV );

            /**
             * Update advection field 
             *
             * @param vH (double) : Horizontal advection velocity [m/s]
             * @param vV (double) : Vertical advection velocity [m/s]
             */

            void UpdateAdv( const double vH, const double vV, const bool artFilt = 0 );
 
            /**
             * Update shear value
             *
             * @param shear (double) : Shear [1/s]
             * @param m (Mesh)
             * */
 
            void UpdateShear( const double shear, const Vector_1D &y, const bool artFilt = 0 );
            void UpdateShear( const Vector_1D shear, const Vector_1D &y, const bool artFilt = 0 );

            /**
             * Solves the 2D advection-diffusion equation over dt using
             * the diffusion and advection fields 
             *
             * @param V (2D vector)         : Field to be diffused 
             * @param cellAreas (2D vector) : Cell areas in m^2
             * @param fillOpt_ (int)        : Fill option
             */

            void Run( Vector_2D &V, const Vector_2D &cellAreas, \
                      const int fillOpt_ = 0 );

            /**
             * Fill value below threshold with value
             *
             * @param V (2D vector)      : Field to be filled
             * @param val (double)       : Filling value
             * @param threshold (double) : Threshold (default = 0.0)
             */

            void Fill( Vector_2D &V, const double val, \
                       const double threshold = 0.0 );
            
            /**
             * Apply correction scheme to get rid of Gibbs oscillations
             *
             * @param V (2D vector)         : Field to be filled
             * @param mass0 (double)        : Mass pre-diffusion/advection
             * @param cellAreas (2D vector) : Cell areas in m^2
             */
    
            void ScinoccaCorr( Vector_2D &V, const double mass0, \
                               const Vector_2D &cellAreas );

            /** 
             * Returns the 2D diffusion field 
             */

            inline const Vector_2D& getDiffFactor( ) const{ return DiffFactor; }
            
            /** 
             * Returns the 2D advection field 
             */

            inline const Vector_2Dc& getAdvFactor( ) const{ return AdvFactor; }

        protected:

            unsigned int n_x, n_y;
            double xlim_right, xlim_left, ylim_up, ylim_down;
            bool doFill;
            UInt fillOpt;
            double fillVal;
            double dt;
            double dH, dV;
            double vH, vV;
            double shear;

        private:

            /* Diffusion and advection field */
            Vector_2D DiffFactor;
            Vector_2Dc AdvFactor;
            Vector_2Dc ShearFactor;

            /* Frequencies */
            Vector_1D kx, ky;
            Vector_1D kxx, kyy;
            
            /* FFT Solver */
            FourierTransform_1D<double> *FFT_1D;
            FourierTransform_2D<double> *FFT_2D;



    };

} /* SANDS */

#endif /* SOLVER_H_INCLUDED */

