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
        xlim( XLIM ),
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

    void Solver::Initialize( const bool fill_, const RealDouble fillVal_, \
                             const UInt fillOpt_ )
    {
    
        FFT_1D = new FourierTransform_1D<RealDouble>( n_x );
        FFT_2D = new FourierTransform_2D<RealDouble>( n_x, n_y );

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

        for ( unsigned int i = 0; i < n_y; i++ ) {
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

        unsigned int i0;
        int k;
     
        /* The domain extends from -xlim to xlim
         * The length of the interval is thus 2*xlim */
    
        /* The frequencies are defined as:
         * kx = 2.0 * PI / (length of interval) * [0:Nx/2-1 -Nx/2:-1]
         */

        i0 = n_x/2; 
        for ( unsigned int i = 0; i < n_x; i++ ) {
            kx.push_back( 0.0 );
            kxx.push_back( 0.0 );
            k = (i0%n_x) - n_x/2;
            kx[i] = physConst::PI / xlim * k;
            kxx[i] = - kx[i] * kx[i];
            i0++;
        }

        i0 = n_y/2;
        for ( unsigned int j = 0; j < n_y; j++ ) {
            ky.push_back( 0.0 );
            kyy.push_back( 0.0 );
            k = (i0%n_y) - n_y/2;
            ky[j] = 2.0 * physConst::PI / ( ylim_down + ylim_up ) * k;
            kyy[j] = - ky[j] * ky[j];
            i0++;
        }

    } /* End of Solver::AssignFreq */
    
    Vector_2D Solver::getDiffFactor( ) const
    {
    
        return DiffFactor;

    } /* End of Solver::getDiffFactor */

    Vector_2Dc Solver::getAdvFactor( ) const
    {
    
        return AdvFactor;

    } /* End of Solver::getAdvFactor */

    void Solver::UpdateTimeStep( const RealDouble T )
    {

        if ( T <= 0.0E+00 ) {
            std::cout << " In Solver::UpdateTimeStep: Non positive time step!\n";
            exit(-1);
        }

        dt = T;

    } /* End of Solver::UpdateTimeStep */

    void Solver::UpdateDiff( const RealDouble dH_, const RealDouble dV_ )
    {

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
            for ( unsigned int iNx = 0; iNx < n_x; iNx++ ) {
                for ( unsigned int jNy = 0; jNy < n_y; jNy++ )
                    DiffFactor[jNy][iNx] = exp( dt * ( dH * kxx[iNx] + dV * kyy[jNy] ) );
            }
        }

    } /* End of Solver::UpdateDiff */

    void Solver::UpdateAdv( const RealDouble vH_, const RealDouble vV_ )
    {

        // vH > 0 means left, < 0 means right
        // vV > 0 means downwards, < 0 means upwards

        if ( ( vH_ != vH ) || ( vV_ != vV ) ) {
            vH = vH_;
            vV = vV_;
            for ( unsigned int iNx = 0; iNx < n_x; iNx++ ) {
                for ( unsigned int jNy = 0; jNy < n_y; jNy++ )
                    AdvFactor[jNy][iNx] = exp( physConst::_1j * dt * ( vH * kx[iNx] + vV * ky[jNy] ) );
            }
        }

    } /* End of Solver::UpdateAdv */

    void Solver::UpdateShear( const RealDouble shear_, const Vector_1D &y )
    {

        /* Declare and initialize horizontal velocity corresponding to shear.
         * This value is dependent on the layer considered. */
        RealDouble V = 0.0E+00;

        if ( shear != shear_ ) {
            /* Update shear value [1/s] */
            shear = shear_;
            for ( unsigned int jNy = 0; jNy < n_y; jNy++ ) {
                /* Compute horizonal velocity. V > 0 means that layer is going left.
                 * Since positive shear means that the plume is rotating 
                 * clockwise when looking from behind the aircraft, add a negative sign.
                 * If meteorological data is symmetric, that should not change the 
                 * outcome */
                V = - shear * y[jNy];
                /* Computing frequencies */
                for ( unsigned int iFreq = 0; iFreq < n_x; iFreq++ )
                    ShearFactor[jNy][iFreq] = exp( physConst::_1j * dt * V * kx[iFreq] ); 
            }
        }

    } /* End of Solver::UpdateShear */

    void Solver::Run( Vector_2D &V, const Vector_2D &cellAreas, const int fillOpt_ )
    {

        RealDouble mass0 = 0.0E+00;

        /* For diagnostic or enforce mass exact conservation, compute mass */
        if ( doFill && fillOpt_ == 1 ) {
            for ( UInt jNy = 0; jNy < n_y; jNy++ ) {
                for ( UInt iNx = 0; iNx < n_x; iNx++ )
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
        if ( shear != 0 )
            FFT_1D->ApplyShear( ShearFactor, V );

        /* 3) Apply corrections */
        /* Fill negative values with fillVal */
        if ( doFill && ( fillOpt_ == 0 ) )
            Fill( V, fillVal );
 
        /* Apply correction scheme to get rid of Gibbs oscillations */
        if ( doFill && ( fillOpt_ == 1 ) )
            ScinoccaCorr( V, mass0, cellAreas );

    } /* End of Solver::Run */

    void Solver::Fill( Vector_2D &V, const RealDouble val, const RealDouble threshold )
    {

        for ( unsigned int iNx = 0; iNx < n_x; iNx++ ) {
            for ( unsigned int jNy = 0; jNy < n_y; jNy++ ) {
                if ( V[jNy][iNx] <= threshold ) 
                    V[jNy][iNx] = val;
            }
        }

    } /* End of Solver::Fill */

    void Solver::ScinoccaCorr( Vector_2D &V, const RealDouble mass0, const Vector_2D &cellAreas )
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
         * V_new = | W + C * ( W - V_low ) * ( V0 - W ),  if W <= V0
         *         | W                                 ,  otherwise 
         *
         * In order for mass to be conserved, we need to have:
         * mass0 = \iint V dA = \iint V_new dA
         *       = \iint W dA + \iint_{W <= V0} C * ( W - V_low ) * ( V0 - W ) dA
         *       = mass    + C \iint_{W <= V0} ( W - V_low ) * ( V0 - W ) dA
         * and thus,
         * C = (mass0 - mass) / \iint_{W <= V0} ( W - V_low ) * ( V0 - W ) dA
         * 
         */

        RealDouble mass = 0.0E+00;
        RealDouble C    = 0.0E+00;
        RealDouble Vlow = 0.0E+00;
        RealDouble V0   = 0.0E+00;

        if ( V[n_y-1][n_x-1] > 0 )
            V0 = V[n_y-1][n_x-1];
        else
            V0 = -V[n_y-1][n_x-1];

        if ( V0 <= 1.00E-100 )
            V0 = 1.00E-100;

        Vlow = V0/1.0E+06;

        for ( UInt jNy = 0; jNy < n_y; jNy++ ) {
            for ( UInt iNx = 0; iNx < n_x; iNx++ ) {
                if ( V[jNy][iNx] <= V0 ) {
                    /* */
                    V[jNy][iNx] = V0 * exp( V[jNy][iNx] / V0 - 1.0 );
                    /* Mass of element lower than V0 in [V]*[m^2] */
                    C += ( V[jNy][iNx] - Vlow ) *\
                         ( V0 - V[jNy][iNx] )   *\
                         cellAreas[jNy][iNx];
                }
                /* Mass of element after removal of negative values 
                 * in [V]*[m^2] */
                mass += V[jNy][iNx] * cellAreas[jNy][iNx];
            }
        }

        if ( C != 0.0E+00 ) {
            /* Correction factor.
             * If mass0 == mass, then correction factor is 0. */
            C  = ( mass0 - mass ) / C;
        }

        /* If correction factor is not finite, set it to 0
         * This is unlikely to make a difference. */
        if ( !std::isfinite(C) || std::isnan(C) )
            C = 0.0E+00;

        for ( UInt jNy = 0; jNy < n_y; jNy++ ) {
            for ( UInt iNx = 0; iNx < n_x; iNx++ ) {
                if ( V[jNy][iNx] <= V0 ) {
                    /* Apply correction */
                    V[jNy][iNx] += C * ( V[jNy][iNx] - Vlow ) \
                                     * ( V0 - V[jNy][iNx] );
                    /* We still want positive values */
                    if ( V[jNy][iNx] < 0.0E+00 )
                        V[jNy][iNx] = 0.0E+00;
                }
            }
        }

    } /* End of Solver::ScinoccaCorr */


} /* SANDS */

/* End of Solver.cpp */
