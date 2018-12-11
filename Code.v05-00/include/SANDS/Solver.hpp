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

            /** 
             * Initializer 
             *
             * Initialize 2D solver
             *
             * @param fill (bool)        : Fill negative values?
             * @param fillValue (double) : Fill with? (default = 0.0E+00) 
             */

            void Initialize( const bool fill = 1, const RealDouble fillValue = 0.0E+00 );

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

            void UpdateTimeStep( const RealDouble T );

            /**
             * Update diffusion field 
             *
             * @param dH (double) : Horizontal diffusion coefficient [m2/s]
             * @param dV (double) : Vertical diffusion coefficient [m2/s]
             *
             * Returns an error if dH or dV is non positive 
             */

            void UpdateDiff( const RealDouble dH, const RealDouble dV );

            /**
             * Update advection field 
             *
             * @param vH (double) : Horizontal advection velocity [m/s]
             * @param vV (double) : Vertical advection velocity [m/s]
             */

            void UpdateAdv( const RealDouble vH, const RealDouble vV );
           
            /**
             * Solves the 2D advection-diffusion equation over dt using
             * the diffusion and advection fields 
             *
             * @param V (2D vector) : Field to be diffused 
             */

            void Run( Vector_2D &V );

            /**
             * Fill value below threshold with value
             *
             * @param V (2D vector)      : Field to be filled
             * @param val (double)       : Filling value
             * @param threshold (double) : Threshold (default = 0.0)
             */

            void Fill( Vector_2D &V, const RealDouble val, const RealDouble threshold = 0.0 );

            /** 
             * Returns the 2D diffusion field 
             */

            Vector_2D getDiffFactor( ) const;
            
            /** 
             * Returns the 2D advection field 
             */

            Vector_2Dc getAdvFactor( ) const;

        protected:

            unsigned int n_x, n_y;
            RealDouble xlim, ylim;
            bool doFill;
            RealDouble fillVal;
            double dt;

        private:

            /* Diffusion and advectoin field */
            Vector_2D DiffFactor;
            Vector_2Dc AdvFactor;
            
            /* Frequencies */
            Vector_1D kx, ky;
            Vector_1D kxx, kyy;
            
            /* FFT Solver */
            FourierTransform<double> *FFT;



    };

} /* SANDS */

#endif /* SOLVER_H_INCLUDED */

