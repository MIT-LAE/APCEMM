/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*              Spectral Advection aNd Diffusion Solver             */
/*                             (SANDS)                              */
/*                                                                  */
/* Solver Program File                                              */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 8/12/2018                                 */
/* File                 : Solver.cpp                                */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "SANDS/Solver.hpp"

namespace SANDS
{

    Solver::Solver( ):
        n_x( NX ),
        n_y( NY ),
        xlim_right( XLIM_RIGHT ),
        xlim_left( XLIM_LEFT ),
        ylim_up( YLIM_UP ),
        ylim_down( YLIM_DOWN ),
        doFill( 1 ),
        fillOpt( 1 ),
        fillVal( 0.0E+00 ),
        FFT_1D( NULL ),
        FFT_2D( NULL )
    {

        /* Constructor */

        shear = 0.0E+00;
        dH    = 0.0E+00;
        dV    = 0.0E+00;
        vH    = 0.0E+00;
        vV    = 0.0E+00;

    } /* End of Solver::Solver */
    Solver::Solver(int nx, int ny,
            double xlim_left, double xlim_right,
            double ylim_up, double ylim_down,
            bool doFill,
            int fillOpt,
            double fillVal): 
    n_x(nx),
    n_y(ny),
    xlim_left(xlim_left),
    xlim_right(xlim_right),
    ylim_up(ylim_up),
    ylim_down(ylim_down),
    doFill(doFill),
    fillOpt(fillOpt),
    fillVal(fillVal),
    shear(0),
    dH(0),
    dV(0),
    vH(0),
    vV(0)
    {

    }
    void Solver::Initialize( const bool MULTITHREADED_FFT, \
                             const bool USE_FFTW_WISDOM,   \
                             const char* FFTW_DIR,         \
                             const bool fill_,             \
                             const double fillVal_,    \
                             const UInt fillOpt_ )
    {
        //Have to do this badness to save the FFTs in the right directory since
        //the code originally didn't use a filesystem library
        std::string fftw_dir_string(FFTW_DIR);
        fftw_dir_string += '/';
        FFT_1D = new FourierTransform_1D<double>( MULTITHREADED_FFT, \
                                                      USE_FFTW_WISDOM,   \
                                                      fftw_dir_string.c_str(),          \
                                                      n_x );
        FFT_2D = new FourierTransform_2D<double>( MULTITHREADED_FFT, \
                                                      USE_FFTW_WISDOM,   \
                                                      fftw_dir_string.c_str(),          \
                                                      n_x,               \
                                                      n_y );

        doFill  = fill_;
        fillOpt = fillOpt_;
        fillVal = fillVal_;

        /* Initialize frequencies */
        AssignFreq();

        /* Initialize shear, diffusion and advection parameters */
        shear = -1.234E+56;
        dH    = -1.234E+56;
        dV    = -1.234E+56;
        vH    = -1.234E+56;
        vV    = -1.234E+56;

        /* Initialize diffusion and advection fields */
        Vector_1D  tempRow( n_x, 0.0E+00 );
        Vector_1Dc tempRowc( n_x, 0.0E+00 );

        for ( UInt i = 0; i < n_y; i++ ) {
            DiffFactor .push_back( tempRow );
            AdvFactor  .push_back( tempRowc );
            ShearFactor.push_back( tempRowc );
        }
    
    } /* End of Solver::Initialize */

    Solver::~Solver( )
    {

        /* Destructor */

        delete FFT_1D;
        /* ^ Calls ~FFT_1D */

        delete FFT_2D;
        /* ^ Calls ~FFT_2D */

    } /* End of Solver::~Solver */

    void Solver::AssignFreq( )
    {

        UInt i0;
        int k;

        /* The domain extends from -xlim to xlim
         * The length of the interval is thus 2*xlim */
    
        /* The frequencies are defined as:
         * kx = 2.0 * PI / (length of interval) * [0:Nx/2-1 -Nx/2:-1]
         */

        /* TODO: Parallelize these blocks */

        i0 = n_x/2;
        for ( UInt i = 0; i < n_x; i++ ) {
            kx.push_back( 0.0 );
            kxx.push_back( 0.0 );
            k = (i0%n_x) - n_x/2;
            kx[i] = 2.0 * physConst::PI / ( xlim_left + xlim_right ) * k;
            kxx[i] = - kx[i] * kx[i];
            i0++;
        }

        i0 = n_y/2;
        for ( UInt j = 0; j < n_y; j++ ) {
            ky.push_back( 0.0 );
            kyy.push_back( 0.0 );
            k = (i0%n_y) - n_y/2;
            ky[j] = 2.0 * physConst::PI / ( ylim_down + ylim_up ) * k;
            kyy[j] = - ky[j] * ky[j];
            i0++;
        }

    } /* End of Solver::AssignFreq */

    void Solver::UpdateTimeStep( const double T )
    {

        if ( T <= 0.0E+00 ) {
            std::cout << " In Solver::UpdateTimeStep: Non positive time step!\n";
            exit(-1);
        }

        dt = T;

    } /* End of Solver::UpdateTimeStep */

    void Solver::UpdateDiff( const double dH_, const double dV_ )
    {

        UInt iNx = 0;
        UInt jNy = 0;

        // dH horizontal diffusion factor >= 0.0
        // dV vertical diffusion factor >= 0.0

        if ( dH_ < 0.0 ) {
            std::cout << "SANDS: Horizontal diffusion coefficient, dH, is negative: dH = " << dH_ << std::endl;
            exit(-1);
        }

        if ( dV_ < 0.0 ) {
            std::cout << "SANDS: Vertical diffusion coefficient, dV, is negative: dV = " << dV_ << std::endl;
            exit(-1);
        }

        if ( ( dH_ != dH ) || ( dV_ != dV ) ) {
            dH = dH_;
            dV = dV_;

            #pragma omp parallel for                \
            if      ( !PARALLEL_CASES ) \
            default ( shared          ) \
            private ( iNx, jNy        ) \
            schedule( dynamic, 1      )
            for ( iNx = 0; iNx < n_x; iNx++ ) {
                for ( jNy = 0; jNy < n_y; jNy++ )
                    DiffFactor[jNy][iNx] = \
                            exp( dt * ( dH * kxx[iNx] + dV * kyy[jNy] ) );
            }
        }

    } /* End of Solver::UpdateDiff */

    void Solver::UpdateAdv( const double vH_, const double vV_, const bool artFilt )
    {

        UInt iNx = 0;
        UInt jNy = 0;
	    double alpha = 0.0;
	    double p = 1.0;
	    Vector_1D etay( n_y, 0.0E+00 );
        double etamaxy = 1.0;
	    Vector_1D etax( n_x, 0.0E+00 );
        double etamaxx = 1.0;

        // vH > 0 means left, < 0 means right
        // vV > 0 means downwards, < 0 means upwards

        /* Set up artificial filtering terms */
	if ( artFilt ) {
	    alpha = 10.0;
	    double hy = ( ylim_up + ylim_down ) / n_y;
	    for ( jNy = 0; jNy < n_y; jNy++ ) {
            if ( jNy < n_y/2 ) {
                etay[jNy] = jNy;
            }
            else {
                etay[jNy] = (double) jNy - (double) n_y;
            }
            // std::cout << "jNy=" << jNy << "etay=" << etay[jNy] << std::endl;
            if ( etay[jNy] > etamaxy ) {
                etamaxy = etay[jNy];
            }
	    }
	    // std::cout << "etamaxy=" << etamaxy << std::endl;
	    double hx = ( xlim_left + xlim_right ) / n_x;
	    for ( iNx = 0; iNx < n_x; iNx++ ) {
            if ( iNx < n_x/2 ) {
                etax[iNx] = iNx;
            }
            else {
                etax[iNx] = (double) iNx - (double) n_x;
            }
		    // std::cout << "iNx=" << iNx << "etax=" << etax[iNx] << std::endl;
            if ( etax[iNx] > etamaxx ) {
                etamaxx = etax[iNx];
            }
	    }
	    // std::cout << "etamaxx=" << etamaxx << std::endl;
	}

    if ( ( vH_ != vH ) || ( vV_ != vV ) ) {
        vH = vH_;
        vV = vV_;

        #pragma omp parallel for                \
        if      ( !PARALLEL_CASES ) \
        default ( shared          ) \
        private ( iNx, jNy        ) \
        schedule( dynamic, 1      )
        for ( iNx = 0; iNx < n_x; iNx++ ) {
            for ( jNy = 0; jNy < n_y; jNy++ ) {
                AdvFactor[jNy][iNx] =                             \
                        exp( physConst::_1j * dt * ( vH * kx[iNx] \
                                                    + vV * ky[jNy] ) \
                    - alpha * pow( etay[jNy]/etamaxy, 2*p ) \
                    - alpha * pow( etax[iNx]/etamaxx, 2*p ) );
            }
        }
    }

} /* End of Solver::UpdateAdv */

    void Solver::UpdateShear( const double shear_, const Vector_1D &y, const bool artFilt )
    {

        UInt iNx = 0;
        UInt jNy = 0;
	    double alpha = 0.0;
	    double p = 1.0;
	    Vector_1D eta( n_x, 0.0E+00 );
        double etamaxx = 1.0;

        /* Declare and initialize horizontal velocity corresponding to shear.
         * This value is dependent on the layer considered. */
        double V = 0.0E+00;

        /* Set up artificial filtering terms */
	    if ( artFilt ) {
	        alpha = 1.0;
	        double hx = ( xlim_left + xlim_right ) / n_x;
            for ( iNx = 0; iNx < n_x; iNx++ ) {
                if ( iNx < n_x/2 ) {
                    eta[iNx] = iNx;
                }
                else {
                    eta[iNx] = (double) iNx - (double) n_x;
                }
                if ( eta[iNx] > etamaxx ) {
                    etamaxx = eta[iNx];
                }
            }
	    }

        if ( shear != shear_ ) {
            /* Update shear value [1/s] */
            shear = shear_;

#pragma omp parallel for               \
            if     ( !PARALLEL_CASES ) \
            default( shared          ) \
            private( iNx, jNy, V     ) \
            schedule( dynamic, 1     )
            for ( jNy = 0; jNy < n_y; jNy++ ) {

                /* Compute horizonal velocity. V > 0 means that layer is going
                 * left. Since positive shear means that the plume is rotating 
                 * clockwise when looking from behind the aircraft, add a
                 * negative sign.
                 * If meteorological data is symmetric, that should not change
                 * the outcome */

                V = - shear * y[jNy];

                /* Computing frequencies */
                for ( iNx = 0; iNx < n_x; iNx++ )
                    ShearFactor[jNy][iNx] = exp( physConst::_1j * dt * V * kx[iNx] - alpha * pow( eta[iNx]/etamaxx, 2*p ) ); 
            }
        }

    } /* End of Solver::UpdateShear */

    void Solver::UpdateShear( const Vector_1D shear_, const Vector_1D &y, const bool artFilt )
    {

        UInt iNx = 0;
        UInt jNy = 0;
	    double alpha = 0.0;
	    double p = 1.0;
	    Vector_1D eta( n_x, 0.0E+00 );
        double etamaxx = 1.0;

        /* Declare and initialize horizontal velocity corresponding to shear.
         * This value is dependent on the layer considered. */
        double V = 0.0E+00;

        /* Set up artificial filtering terms */
	if ( artFilt ) {
	    alpha = 1.0;
	    double hx = ( xlim_left + xlim_right ) / n_x;
	    for ( iNx = 0; iNx < n_x; iNx++ ) {
            if ( iNx < n_x/2 ) {
                eta[iNx] = iNx;
            }
            else {
                eta[iNx] = (double) iNx - (double) n_x;
            }
            if ( eta[iNx] > etamaxx ) {
                etamaxx = eta[iNx];
            }
	    }
	}

#pragma omp parallel for               \
        if     ( !PARALLEL_CASES ) \
        default( shared          ) \
        private( iNx, jNy, V     ) \
        schedule( dynamic, 1     )
        for ( jNy = 0; jNy < n_y; jNy++ ) {

                /* Compute horizonal velocity. V > 0 means that layer is going
                 * left. Since positive shear means that the plume is rotating 
                 * clockwise when looking from behind the aircraft, add a
                 * negative sign.
                 * If meteorological data is symmetric, that should not change
                 * the outcome */

                shear = shear_[jNy];
                V = - shear * y[jNy];
                /* Computing frequencies */
                for ( iNx = 0; iNx < n_x; iNx++ )
                    ShearFactor[jNy][iNx] = exp( physConst::_1j * dt * V * kx[iNx] - alpha * pow( eta[iNx]/etamaxx, 2*p ) ); 
            
        }

    } /* End of Solver::UpdateShear */

    void Solver::Run( Vector_2D &V, const Vector_2D &cellAreas, const int fillOpt_ )
    {

        UInt iNx = 0;
        UInt jNy = 0;

        double mass0   = 0.0E+00;

        /* For diagnostic or enforce mass exact conservation, compute mass */
        if ( doFill && fillOpt_ == 1 ) {
#pragma omp parallel for                 \
            if       ( !PARALLEL_CASES ) \
            default  ( shared          ) \
            private  ( iNx, jNy        ) \
            reduction( +:mass0         ) \
            schedule ( dynamic, 1      )
            for ( jNy = 0; jNy < n_y; jNy++ ) {
                for ( iNx = 0; iNx < n_x; iNx++ )
                    mass0 += V[jNy][iNx] * cellAreas[jNy][iNx];
            }
        }

        /* Operator splitting approach:
         * 1) Solve outward diffusion and advection/settling
         * 2) Apply shear forces
         * 3) Apply any correction (positivity, Gibbs oscillation removal) 
         */

        /* 1) Apply diffusion and settling */
        FFT_2D->SANDS( DiffFactor, AdvFactor, V );

        /* 2) Apply shear forces */
        if ( shear != 0 ) {
            FFT_1D->ApplyShear( ShearFactor, V );
        }

        /* 3) Apply corrections */
        /* Fill negative values with fillVal */
        if ( doFill && ( fillOpt_ == 0 ) ) {
            Fill( V, fillVal );
        }

        /* Apply correction scheme to get rid of Gibbs oscillations */
        if ( doFill && ( fillOpt_ == 1 ) ) {
            ScinoccaCorr( V, mass0, cellAreas );
        }

    } /* End of Solver::Run */

    void Solver::Fill( Vector_2D &V, const double val, \
                       const double threshold )
    {

        UInt iNx = 0;
        UInt jNy = 0;

#pragma omp parallel for                \
            if      ( !PARALLEL_CASES ) \
            default ( shared          ) \
            private ( iNx, jNy        ) \
            schedule( dynamic, 1      )
        for ( iNx = 0; iNx < n_x; iNx++ ) {
            for ( jNy = 0; jNy < n_y; jNy++ ) {
                if ( V[jNy][iNx] <= threshold ) 
                    V[jNy][iNx] = val;
            }
        }

    } /* End of Solver::Fill */

    void Solver::ScinoccaCorr( Vector_2D &V, const double mass0, \
                               const Vector_2D &cellAreas )
    {

        /* The correction scheme is based on:
         * Scinocca, J. F., et al. "The CCCma third generation AGCM and its 
         * extension into the middle atmosphere." 
         * Atmospheric Chemistry and Physics 8.23 (2008): 7055-7074.*/

        /* 
         * Let V denote the data before advection and W the result of V 
         * after advection, solved using a spectral solver.
         *
         * Ensure that data is positive
         * W = | V0 * exp( V / V0 - 1 ),  if V <= V0
         *     | V                     ,  otherwise 
         *
         * Ensure that mass is conserved by applying correction
         * V_new = | W + C * ( W - V_low ) * ( V0 - W ),  if V_low <= W <= V0
         *         | W                                 ,  otherwise 
         *
         * In order for mass to be conserved, we need to have:
         * mass0 = \iint V dA = \iint V_new dA
         *       = \iint W dA + \iint_{Vlow <= W <= V0} C * ( W - V_low ) * ( V0 - W ) dA
         *       = mass    + C \iint_{Vlow <= W <= V0} ( W - V_low ) * ( V0 - W ) dA
         * and thus,
         * C = (mass0 - mass) / \iint_{Vlow <= W <= V0} ( W - V_low ) * ( V0 - W ) dA
         * 
         */

        UInt iNx = 0;
        UInt jNy = 0;

        double mass    = 0.0E+00;
        double C       = 0.0E+00;
        double Vlow    = 0.0E+00;
        double V0      = 0.0E+00;
        double negMass = 0.0E+00;

        if ( V[n_y-1][n_x-1] > 0 )
            V0 = V[n_y-1][n_x-1];
        else
            V0 = -V[n_y-1][n_x-1];

        if ( V0 <= 1.00E-70 )
            V0 = 1.00E-70;

        Vlow = V0/1.0E+06;

#pragma omp parallel default(shared) if ( !PARALLEL_CASES )
        {

#pragma omp for                          \
            reduction( +:C             ) \
            reduction( +:mass          ) \
            private  ( iNx, jNy        ) \
            schedule ( dynamic, 1      )
        for ( jNy = 0; jNy < n_y; jNy++ ) {
            for ( iNx = 0; iNx < n_x; iNx++ ) {
                if ( V[jNy][iNx] <= V0 ) {
                    /* */
                    V[jNy][iNx] = V0 * exp( V[jNy][iNx] / V0 - 1.0 );

                    if ( V[jNy][iNx] >= Vlow ) {
                        /* Mass of element between Vlow and V0 in [V]*[m^2] */
                        C += ( V[jNy][iNx] - Vlow ) *\
                             ( V0 - V[jNy][iNx] )   *\
                             cellAreas[jNy][iNx];
                    }
                }

                /* Mass of element after removal of low values ( <= V0 )
                 * in [V]*[m^2] */
                mass += V[jNy][iNx] * cellAreas[jNy][iNx];
            }
        }

#pragma omp single 
        {
        if ( C != 0.0E+00 ) {
            /* Correction factor. C should always be <= 0
             * If mass0 == mass, then correction factor is 0. */
            C = ( mass0 - mass ) / C;
        }

        /* If correction factor is not finite, set it to 0
         * This is unlikely to make a difference. */
        if ( !std::isfinite(C) || std::isnan(C) )
            C = 0.0E+00;
        }

#pragma omp barrier

#pragma omp for                          \
            reduction( +:negMass       ) \
            private  ( iNx, jNy        ) \
            schedule ( dynamic, 1      )
        for ( jNy = 0; jNy < n_y; jNy++ ) {
            for ( iNx = 0; iNx < n_x; iNx++ ) {
                if ( ( V[jNy][iNx] <= V0 ) && ( V[jNy][iNx] >= Vlow ) ) {

                    /* Apply correction, this can yield negative values,
                     * especially if Vlow + C * V0^2 < 0 */
                    V[jNy][iNx] += C * ( V[jNy][iNx] - Vlow ) \
                                     * ( V0 - V[jNy][iNx] );

                    /* We still want positive values, compute negative mass */
                    if ( V[jNy][iNx] < 0.0E+00 ) {
                        negMass     -= V[jNy][iNx] * cellAreas[jNy][iNx];
                        V[jNy][iNx]  = 0.0E+00;
                    }

                }
            }
        }

        } /* End of pragma omp parallel */

        if ( negMass > 0.0E+00 ) {

            /* The following lines aim to correct for the mass that has been
             * removed by clipping negative values.
             * If this correction were not applied, exact mass correction would not
             * be ensured. */
            bool success     = 0;
            UInt counter     = 0;
            int guess       = n_y*n_x;
            double tArea = 0.0E+00;
            
            while ( !success ) {
                counter = 0;
                tArea   = 0.0E+00;

                for ( UInt jNy = 0; jNy < n_y; jNy++ ) {
                    for ( UInt iNx = 0; iNx < n_x; iNx++ ) {
                        if ( V[jNy][iNx] > negMass / ( guess * cellAreas[jNy][iNx] ) ) {
                            counter += 1;
                            tArea   += cellAreas[jNy][iNx];
                        }
                    }
                }
                
                if ( counter >= guess )
                    success = 1;
                else
                    guess  -= (int) n_y*n_x/16;


                if ( guess <= 0 ) {
                    tArea = 0.0E+00;
                    break;
                }

            }
            

            if ( tArea > 0.0E+00 ) {

                for ( jNy = 0; jNy < n_y; jNy++ ) {
                    for ( iNx = 0; iNx < n_x; iNx++ ) {
                        if ( V[jNy][iNx] > negMass / ( guess * cellAreas[jNy][iNx] ) )
                            V[jNy][iNx] -= negMass / tArea;
                    }
                }

            } else {

                /* For debugging purposes, uncomment the following line */
                // std::cout << "tArea = " << tArea << " (dil: " << negMass / cellAreas[1][1] << ")" << std::endl;

            }

        }

    } /* End of Solver::ScinoccaCorr */


} /* SANDS */

/* End of Solver.cpp */
