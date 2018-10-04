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

    const char *Solver::wisdomFile = WISDOMFILE;

    Solver::Solver( )
    {
        /* Base Constructor */

        bool fill = 0;
        RealDouble fillValue = 0.0;
        unsigned flag = FFTW_ESTIMATE;

        Solver( fill, fillValue, flag );

    } /* End of Solver::Solver */

    Solver::Solver( bool fill, RealDouble fillValue, unsigned flag )
    {
        /* Constructor */

        n_x = NX;
        n_y = NY;

        xlim = XLIM;
        ylim = YLIM;

        doFill = fill;
        fillVal = fillValue;
        FFTW_flag = flag; 

        AssignFreq();

        for ( unsigned int i = 0; i < n_y; i++ ) {
            DiffFactor.push_back( Real_1DVector( n_x ) );
            AdvFactor .push_back( Complex_1DVector( n_x ) );
        }

    } /* End of Solver::Solver */

    Solver::Solver( const Solver &s )
    {

        n_x = s.n_x;
        n_y = s.n_y;

        xlim = s.xlim;
        ylim = s.ylim;

        doFill = s.doFill;
        fillVal = s.fillVal;
        
        FFTW_flag = s.FFTW_flag;
        dt = s.dt;

        DiffFactor = s.DiffFactor;
        AdvFactor = s.AdvFactor;

        kx = s.kx;
        ky = s.ky;
        kxx = s.kxx;
        kyy = s.kyy;

    } /* End of Solver::Solver */

    Solver& Solver::operator=( const Solver &s )
    {

        if ( &s == this )
            return *this;

        n_x = s.n_x;
        n_y = s.n_y;

        xlim = s.xlim;
        ylim = s.ylim;

        doFill = s.doFill;
        fillVal = s.fillVal;
        
        FFTW_flag = s.FFTW_flag;
        dt = s.dt;

        DiffFactor = s.DiffFactor;
        AdvFactor = s.AdvFactor;

        kx = s.kx;
        ky = s.ky;
        kxx = s.kxx;
        kyy = s.kyy;
        return *this;

    } /* End of Solver::operator= */

    Solver::~Solver( )
    {
        /* Destructor */

    } /* End of Solver::~Solver */

    void Solver::AssignFreq( )
    {

        unsigned int i0;
        int k;
     
        /* The domain extends from -xlim to xlim
         * so, length of the interval is 2*xlim */
        /* Formula is:
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

        /* Allocate n_y elements */
        
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

    
    Real_2DVector Solver::getDiffFactor( ) const
    {
    
        return DiffFactor;

    } /* End of Solver::getDiffFactor */

    Complex_2DVector Solver::getAdvFactor( ) const
    {
    
        return AdvFactor;

    } /* End of Solver::getAdvFactor */

    void Solver::UpdateTimeStep( RealDouble T )
    {

        dt = T;

    } /* End of Solver::UpdateStep */

    void Solver::UpdateDiff( RealDouble dH, RealDouble dV )
    {

        // dH horizontal diffusion factor >= 0.0
        // dV vertical diffusion factor >= 0.0

        if ( dH < 0.0 ) {
            std::cout << "SANDS: Horizontal diffusion coefficient, dH, is negative: dH = " << dH << std::endl;
            return;
        }

        if ( dV < 0.0 ) {
            std::cout << "SANDS: Vertical diffusion coefficient, dV, is negative: dV = " << dV << std::endl;
            return;
        }

        for ( unsigned int i = 0; i < n_x; i++ ) {
            for ( unsigned int j = 0; j < n_y; j++ )
                DiffFactor[j][i] = exp( dt * ( dH * kxx[i] + dV * kyy[j] ) );
        }    

    } /* End of Solver::UpdateDiff */

    void Solver::UpdateAdv( double vH, double vV )
    {

        // vH > 0 means left, < 0 means right
        // vV > 0 means upwards, < 0 means downwards
        
        for ( unsigned int i = 0; i < n_x; i++ ) {
            for ( unsigned int j = 0; j < n_y; j++ )
                AdvFactor[j][i] = exp( _1j * dt * ( vH * kx[i] + vV * ky[j] ) );
        }

    } /* End of Solver::UpdateAdv */

    void Solver::Solve( Real_2DVector &V, const bool realInput )
    {

       SANDS( V, DiffFactor, AdvFactor, wisdomFile, realInput ); 

       if ( doFill )
           Fill( V, fillVal );

    } /* End of Solver::SANDS */

    void Solver::Fill( Real_2DVector &V, RealDouble val, RealDouble threshold )
    {

        for ( unsigned int iNx = 0; iNx < n_x; iNx++ ) {
            for ( unsigned int jNy = 0; jNy < n_y; jNy++ ) {
                if ( V[jNy][iNx] <= threshold ) 
                    V[jNy][iNx] = val;
            }
        }
        

    } /* End of Solver::Fill */


    void Solver::Wisdom( Real_2DVector &V )
    {

        std::ifstream ifile( wisdomFile );

        if ( (bool)ifile ) 
            std::cout << "Wisdom file: " << wisdomFile << " will be overwritten! " << std::endl;

        SaveWisdom( V, wisdomFile );

    } /* End of Solver::Wisdom */

    unsigned int Solver::getNx() const
    {

        return n_x;

    } /* End of Solver::getNx */

    unsigned int Solver::getNy() const
    {

        return n_y;

    } /* End of Solver::getNy */

    RealDouble Solver::getXlim() const
    {

        return xlim;

    } /* End of Solver::getXlim */

    RealDouble Solver::getYlim() const
    {

        return ylim;

    } /* End of Solver::getYlim */

    RealDouble Solver::getDt() const
    {

        return dt;

    } /* End of Solver::getDt */

/* End of Solver.cpp */
