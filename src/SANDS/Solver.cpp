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
        ylim( YLIM ),
        doFill( 1 ), 
        fillVal( 0.0E+00 ),
        FFT( NULL )
    {

        /* Constructor */

    } /* End of Solver::Solver */

    void Solver::Initialize( const bool fill, const RealDouble fillValue )
    {
    
        FFT = new FourierTransform<double>( n_x, n_y );

        doFill = fill;
        fillVal = fillValue;

        /* Initialize frequencies */
        AssignFreq();

        /* Initialize diffusion and advection fields */
        for ( unsigned int i = 0; i < n_y; i++ ) {
            DiffFactor.push_back( Vector_1D( n_x ) );
            AdvFactor .push_back( Vector_1Dc( n_x ) );
        }
    
    } /* End of Solver::Initialize */

    Solver::~Solver( )
    {

        /* Destructor */

        delete FFT;
        /* ^ Calls ~FFT */

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
            ky[j] = physConst::PI / ylim * k;
            kyy[j] = - ky[j] * ky[j];
            i0++;
        }

    }

    
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

    void Solver::UpdateDiff( const RealDouble dH, const RealDouble dV )
    {

        // dH horizontal diffusion factor >= 0.0
        // dV vertical diffusion factor >= 0.0

        if ( dH < 0.0 ) {
            std::cout << "SANDS: Horizontal diffusion coefficient, dH, is negative: dH = " << dH << std::endl;
            exit(-1);
        }

        if ( dV < 0.0 ) {
            std::cout << "SANDS: Vertical diffusion coefficient, dV, is negative: dV = " << dV << std::endl;
            exit(-1);
        }

        for ( unsigned int i = 0; i < n_x; i++ ) {
            for ( unsigned int j = 0; j < n_y; j++ )
                DiffFactor[j][i] = exp( dt * ( dH * kxx[i] + dV * kyy[j] ) );
        }    

    } /* End of Solver::UpdateDiff */

    void Solver::UpdateAdv( const double vH, const double vV )
    {

        // vH > 0 means left, < 0 means right
        // vV > 0 means upwards, < 0 means downwards
        
        for ( unsigned int i = 0; i < n_x; i++ ) {
            for ( unsigned int j = 0; j < n_y; j++ )
                AdvFactor[j][i] = exp( physConst::_1j * dt * ( vH * kx[i] + vV * ky[j] ) );
        }

    } /* End of Solver::UpdateAdv */

    void Solver::Run( Real_2DVector &V )
    {

        FFT->SANDS( DiffFactor, AdvFactor, V );

        if ( doFill )
            Fill( V, fillVal );


    } /* End of Solver::Run */

    void Solver::Fill( Real_2DVector &V, const RealDouble val, const RealDouble threshold )
    {

        for ( unsigned int iNx = 0; iNx < n_x; iNx++ ) {
            for ( unsigned int jNy = 0; jNy < n_y; jNy++ ) {
                if ( V[jNy][iNx] <= threshold ) 
                    V[jNy][iNx] = val;
            }
        }
        

    } /* End of Solver::Fill */


} /* SANDS */

/* End of Solver.cpp */
