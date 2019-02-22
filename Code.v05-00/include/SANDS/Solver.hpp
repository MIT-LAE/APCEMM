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
             * @param fill (bool)      : Fill negative values?
             * @param fillVal (double) : Fill with? (default = 0.0E+00) 
             * @param fillOpt (int)    : Fill option
             */

            void Initialize( const bool fill = 1, \
                             const RealDouble fillValue = 0.0E+00, \
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
             * Update shear value
             *
             * @param shear (double) : Shear [1/s]
             * @param m (Mesh)
             * */
 
            void UpdateShear( const RealDouble shear, const Vector_1D &y );

            /**
             * Solves the 2D advection-diffusion equation over dt using
             * the diffusion and advection fields 
             *
             * @param V (2D vector)         : Field to be diffused 
             * @param cellAreas (2D vector) : Cell areas in m^2
             * @param fillOpt_ (UInt)       : Fill option
             */

            void Run( Vector_2D &V, const Vector_2D &cellAreas, \
                      const UInt fillOpt_ = 0 );

            /**
             * Fill value below threshold with value
             *
             * @param V (2D vector)      : Field to be filled
             * @param val (double)       : Filling value
             * @param threshold (double) : Threshold (default = 0.0)
             */

            void Fill( Vector_2D &V, const RealDouble val, const RealDouble threshold = 0.0 );
            
            /**
             * Apply correction scheme to get rid of Gibbs oscillations
             *
             * @param V (2D vector)         : Field to be filled
             * @param mass0 (double)        : Mass pre-diffusion/advection
             * @param cellAreas (2D vector) : Cell areas in m^2
             */
    
            void ScinoccaCorr( Vector_2D &V, const RealDouble mass0, const Vector_2D &cellAreas );

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
            UInt fillOpt;
            RealDouble fillVal;
            RealDouble dt;
            RealDouble shear;

        private:

            /* Diffusion and advection field */
            Vector_2D DiffFactor;
            Vector_2Dc AdvFactor;
            Vector_2Dc ShearFactor;

            /* Frequencies */
            Vector_1D kx, ky;
            Vector_1D kxx, kyy;
            
            /* FFT Solver */
            FourierTransform_1D<RealDouble> *FFT_1D;
            FourierTransform_2D<RealDouble> *FFT_2D;



    };

} /* SANDS */

#endif /* SOLVER_H_INCLUDED */

