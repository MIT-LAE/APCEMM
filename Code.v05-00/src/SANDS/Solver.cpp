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
        FFT_1D( NULL ),
        FFT_2D( NULL )
    {


        /* Constructor */

    } /* End of Solver::Solver */

    void Solver::Initialize( const bool fill, const RealDouble fillValue )
    {
    
        FFT_1D = new FourierTransform_1D<double>( n_x );
        FFT_2D = new FourierTransform_2D<double>( n_x, n_y );

        doFill = fill;
        fillVal = fillValue;

        /* Initialize frequencies */
        AssignFreq();

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

    void Solver::UpdateShear( const double shear_, const Mesh &m )
    {

        /* Declare and initialize horizontal velocity corresponding to shear.
         * This value is dependent on the layer considered. */
        double V = 0.0E+00;

        /* Update shear value [1/s] */
        shear = shear_;

        /* Get coordinates of each horizontal line */
        Vector_1D y = m.y();

        for ( unsigned int jNy = 0; jNy < m.Ny(); jNy++ ) {
            /* Compute horizonal velocity. V > 0 means that layer is going left */
            V = shear * y[jNy];
            /* Computing frequencies */
            for ( unsigned int iFreq = 0; iFreq < m.Nx(); iFreq++ )
                ShearFactor[jNy][iFreq] = exp( physConst::_1j * dt * V * kx[iFreq] ); 
        }

    } /* End of Solver::UpdateShear */

    void Solver::Run( Real_2DVector &V )
    {

        /* Operator splitting approach:
         * 1) solve outward diffusion and advection/settling
         * 2) Apply shear forces
         */

        /* 1) Apply diffusion and settling */
        FFT_2D->SANDS( DiffFactor, AdvFactor, V );

        /* 2) Apply shear forces */
        FFT_1D->ApplyShear( ShearFactor, V );

        /* Fill negative values with fillVal */
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
