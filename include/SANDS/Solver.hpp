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
#include "Util/ForwardsDecl.hpp"

    typedef fftw_complex FFTW_ComplexDouble;

    const ComplexDouble _1j ( 0.0, 1.0 );          /* j^2 = -1 */

    void SANDS( Real_2DVector &V, Real_2DVector &diff, Complex_2DVector &adv, \
                const char* fileName_FFTW, const bool realInput );
    void SaveWisdom( Real_2DVector &V, const char* fileName_FFTW );

    /* A class for solving the 2D-Avection Diffusion equation */

    class Solver
    {

        public:
            
            /* Base constructor  */
            Solver( );

            /* Constructor */
            Solver( const bool fill = 0, const RealDouble fillValue = 0.0, const unsigned flag = FFTW_ESTIMATE );

            /* Copy constructor */
            /*
             * Parameter s to be copied
             */
            Solver( const Solver &s );

            /* Destructor */
            ~Solver( );

            void AssignFreq( );

            /* Update time step */
            void UpdateTimeStep( const RealDouble T );

            /* Update diffusion field */
            void UpdateDiff( const RealDouble dH, const RealDouble dV );

            /* Update advection field */
            void UpdateAdv( const RealDouble vH, const RealDouble vV );
           
            void Solve( Real_2DVector &V, const bool realInput );

            /* Fill value between threshold with val */
            void Fill( Real_2DVector &V, const RealDouble val, const RealDouble threshold = 0.0 );

            void Wisdom( Real_2DVector &V );

            Real_2DVector getDiffFactor( ) const;

            Complex_2DVector getAdvFactor( ) const;

            unsigned int getNx() const;
            unsigned int getNy() const;
            RealDouble getXlim() const;
            RealDouble getYlim() const;
            RealDouble getDt() const;

        protected:

            unsigned int n_x, n_y;
            RealDouble xlim, ylim;
            const bool doFill;
            const RealDouble fillVal;
            unsigned FFTW_flag;
            double dt;

        private:

            static const char *wisdomFile;
            Real_2DVector DiffFactor;
            Complex_2DVector AdvFactor;
            Real_1DVector kx, ky, kxx, kyy;
            

    };

#endif /* SOLVER_H_INCLUDED */

