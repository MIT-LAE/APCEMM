/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Mesh Program File                                                */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Mesh.cpp                                  */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Mesh.hpp"

Mesh::Mesh( )
{

    /* Default Constructor */

    nx = NX;
    ny = NY;

    xlim = XLIM;
    ylim = YLIM;

    hx = 2 * xlim / nx;
    hy = 2 * ylim / ny;

    /* Cell center x-coordinates */
    for ( unsigned int i = 0; i < nx; i++ ) {
        x.push_back( (double) 0.0 );
        x[i] = i * hx - xlim + hx / 2;
        x_e.push_back( (double) 0.0 );
        x_e[i] = x[i] - hx/2;
    }
    x_e.push_back( (double) 0.0 );
    x_e[nx] = x_e[nx-1] + hx;
    
    /* Cell center y-coordinates */
    for ( unsigned int j = 0; j < ny; j++ ) {
        y.push_back( (double) 0.0 );
        y[j] = j * hy - ylim + hy / 2;
        y_e.push_back( (double) 0.0 );
        y_e[j] = y[j] - hy/2;
    }
    y_e.push_back( (double) 0.0 );
    y_e[ny] = y_e[ny-1] + hy;

    totArea = 0;
    for ( unsigned int jNy = 0; jNy < ny; jNy++ ) {
        areas.push_back( Real_1DVector( nx, 0.0 ) );
        for ( unsigned int iNx = 0; iNx < nx; iNx++ ) {
            areas[jNy][iNx] = ( y_e[jNy+1] - y_e[jNy] ) * \
                              ( x_e[iNx+1] - x_e[iNx] );
            /* Comes down to hx * hy, if mesh is cartesian uniform */
            totArea += areas[jNy][iNx];
        }
    }

} /* End of Mesh::Mesh */

Mesh::Mesh( const Mesh &m )
{

    x = m.x;
    y = m.y;
    xlim = m.xlim;
    ylim = m.ylim;
    hx = m.hx;
    hy = m.hy;
    nx = m.nx;
    ny = m.ny;

} /* End of Mesh::Mesh */

Mesh& Mesh::operator=( const Mesh &m )
{

    if ( &m == this )
        return *this;

    x = m.x;
    y = m.y;
    xlim = m.xlim;
    ylim = m.ylim;
    hx = m.hx;
    hy = m.hy;
    nx = m.nx;
    ny = m.ny;
    return *this;

} /* End of Mesh::operator= */

Mesh::~Mesh( )
{

    /* Destructor */

} /* End of Mesh::~Mesh */

void Mesh::Ring2Mesh( Cluster &c )
{
    
    unsigned int nRing = c.getnRing();

    std::vector<std::vector<bool >> v2d;
    std::vector<bool> v1d;
    std::vector<std::pair<unsigned int, unsigned int>> v1d_int;
    v1d = std::vector<bool>( NX, 0 );

    unsigned int nx_max, ny_max;

    /* Assert that NX and NY are multiples of 2! */
#if Y_SYMMETRY
    nx_max = std::ceil(NX/2);
#else
    nx_max = NX;
#endif /* Y_SYMMETRY */

#if X_SYMMETRY
    ny_max = std::ceil(NY/2);
#else
    ny_max = NY;
#endif /* X_SYMMETRY */

    for ( unsigned int iRing = 0; iRing < nRing; iRing++ ) {
        nCellMap.push_back( 0 );
        RingMeshMap.push_back( v2d );
        indList.push_back( v1d_int );
        for ( unsigned int iNy = 0; iNy < NY; iNy++ ) {
            RingMeshMap[iRing].push_back( v1d );
        }
    }

    std::vector<Ring> RingV;
    RingV = c.getRings();
    double hAxis, vAxis, hAxis_in, vAxis_in;
    double xRatio, xRatio_in;
    unsigned int val;
    for ( unsigned int iRing = 0; iRing < nRing; iRing++ ) {
        val = 0;
        if ( iRing == 0 ) {
            hAxis = RingV[iRing].getHAxis();
            vAxis = RingV[iRing].getVAxis();
            /* Use symmetry */
            for ( unsigned int iNx = 0; iNx < nx_max; iNx++ ) {
                xRatio = ( x[iNx]  /  hAxis ) * ( x[iNx]  /  hAxis );
                for ( unsigned int jNy = 0; jNy < ny_max; jNy++ ) {
                    /* If point is within the region and out of the smaller region */
                    if ( xRatio + ( y[jNy] / vAxis ) * ( y[jNy] / vAxis ) <= 1 ) {
                        RingMeshMap[iRing][jNy][iNx] = 1;
                        indList[iRing].push_back( std::make_pair(iNx, jNy) );
                        val++;
                    }
                }
            }
        }
        else {
            hAxis = RingV[iRing].getHAxis();
            vAxis = RingV[iRing].getVAxis();
            hAxis_in = RingV[iRing-1].getHAxis();
            vAxis_in = RingV[iRing-1].getVAxis();
            for ( unsigned int iNx = 0; iNx < nx_max; iNx++ ) {
                xRatio    = ( x[iNx]  /  hAxis    ) * ( x[iNx]  /  hAxis    );
                xRatio_in = ( x[iNx]  /  hAxis_in ) * ( x[iNx]  /  hAxis_in );
                for ( unsigned int jNy = 0; jNy < ny_max; jNy++ ) {
                    /* If point is within the region and out of the smaller region */
                    if ( ( xRatio + ( y[jNy] / vAxis ) * ( y[jNy] / vAxis ) <= 1 ) && ( xRatio_in + ( y[jNy] / vAxis_in ) * ( y[jNy] / vAxis_in ) > 1 ) ) {
                        RingMeshMap[iRing][jNy][iNx] = 1;
                        indList[iRing].push_back( std::make_pair(iNx, jNy) );
                        val++;
                    }
                }
            }
        }
        if ( val != 0 )
            nCellMap[iRing] = val;
        else
            std::cout << "Ring " << iRing << " has no cell in it (nMap = 0)!!" << std::endl;

#if ( Y_SYMMETRY && X_SYMMETRY )
        /* If both symmetries are valid, loop over 3 regions */
        /*
         *                         |
         *              (2)        |       (4)
         *                         |
         * NY/2  __________________|__________________
         *                         |
         *                         |
         *              (1)        |       (3)
         *                         |
         * (0,0)                 NX/2
         */
        /* Region 1 has been done previously */

        /* Do region 2 */
        for ( unsigned int iNx = 0; iNx < nx_max; iNx++ ) {
            for ( unsigned int jNy = ny_max; jNy < NY; jNy++ ) {
                RingMeshMap[iRing][jNy][iNx] = RingMeshMap[iRing][(NY - 1) - jNy][iNx];
            }
        }
        
        for ( unsigned int iList = 0; iList < nCellMap[iRing]; iList++ ) {
            indList[iRing].push_back( std::make_pair(indList[iRing][iList].first, (NY - 1) - indList[iRing][iList].second) );
        }
        nCellMap[iRing] *= 2;
        
        /* Do regions 3 and 4 */
        /* 3 */
        for ( unsigned int iNx = nx_max; iNx < NX; iNx++ ) {
            for ( unsigned int jNy = 0; jNy < ny_max; jNy++ ) {
                RingMeshMap[iRing][jNy][iNx] = RingMeshMap[iRing][jNy][(NX - 1) - iNx];
            }
        }
       
        /* 4 */
        for ( unsigned int iNx = nx_max; iNx < NX; iNx++ ) {
            for ( unsigned int jNy = ny_max; jNy < NY; jNy++ ) {
                RingMeshMap[iRing][jNy][iNx] = RingMeshMap[iRing][(NY - 1) - jNy][(NX - 1) - iNx];
            }
        }
        
        for ( unsigned int iList = 0; iList < nCellMap[iRing]; iList++ ) {
            indList[iRing].push_back( std::make_pair((NX - 1) - indList[iRing][iList].first, indList[iRing][iList].second) );
        }
        nCellMap[iRing] *= 2;

#endif

#if ( X_SYMMETRY && !Y_SYMMETRY )
        /* Do regions 2 and 4 */
        /* 2 and 4 */
        for ( unsigned int iNx = 0; iNx < NX; iNx++ ) {
            for ( unsigned int jNy = ny_max; jNy < NY; jNy++ ) {
                RingMeshMap[iRing][jNy][iNx] = RingMeshMap[iRing][(NY - 1) - jNy][iNx];
            }
        }
        
        for ( unsigned int iList = 0; iList < nCellMap[iRing]; iList++ ) {
            indList[iRing].push_back( std::make_pair(indList[iRing][iList].first, (NY - 1) - indList[iRing][iList].second) );
        }
        nCellMap[iRing] *= 2;

#endif

#if ( !X_SYMMETRY && Y_SYMMETRY )
        /* Do regions 3 and 4 */
        /* 3 and 4 */
        for ( unsigned int iNx = nx_max; iNx < NX; iNx++ ) {
            for ( unsigned int jNy = 0; jNy < NY; jNy++ ) {
                RingMeshMap[iRing][jNy][iNx] = RingMeshMap[iRing][jNy][(NX - 1) - iNx];
            }
        }

        for ( unsigned int iList = 0; iList < nCellMap[iRing]; iList++ ) {
            indList[iRing].push_back( std::make_pair((NX - 1) - indList[iRing][iList].first, indList[iRing][iList].second) );
        }
        nCellMap[iRing] *= 2;

#endif

    }


} /* End of Mesh::Ring2Mesh */

Real_1DVector Mesh::getX( ) const
{

    return x;

} /* End of Mesh::getX */

Real_1DVector Mesh::getY( ) const
{

    return y;

} /* End of Mesh::getY */
        
Real_1DVector Mesh::getX_e( ) const
{

    return x_e;

} /* End of Mesh::getX_e */

Real_1DVector Mesh::getY_e( ) const
{

    return y_e;

} /* End of Mesh::getY_e */

Real_2DVector Mesh::getAreas( ) const
{

    return areas;

} /* End of Mesh::getAreas */

RealDouble Mesh::getTotArea( ) const
{

    return totArea;

} /* End of Mesh::getTotArea */

RealDouble Mesh::gethx( ) const
{

    return hx;

} /* End of Mesh::gethx */

RealDouble Mesh::gethy( ) const
{

    return hy;

} /* End of Mesh::gethy */

unsigned int Mesh::getNx( ) const
{

    return nx;

} /* End of Mesh::getNx */

unsigned int Mesh::getNy( ) const
{

    return ny;

} /* End of Mesh::getNy */

std::vector<std::vector<std::vector<bool> > > Mesh::getMap( ) const
{

    return RingMeshMap;

} /* End of Mesh::getMap */

std::vector<std::vector<std::pair<unsigned int, unsigned int> > > Mesh::getList( ) const
{

    return indList;

} /* End of Mesh::getList */

std::vector<unsigned int> Mesh::getnMap( ) const
{

    return nCellMap;

} /* End of Mesh::getnMap */

void Mesh::Debug( ) const
{

    std::cout << std::endl;
    std::cout << "**** Input Debugger ****" << std::endl;
    std::cout << "Mesh-Ring structure: " << std::endl;
    std::cout << std::endl;

    std::cout << " Number of cells in each ring: " << std::endl;
    for ( unsigned int iRing = 0; iRing < nCellMap.size(); iRing++ ) {
        std::cout << " ";
        std::cout << std::setw(4);
        std::cout << nCellMap[iRing];
        std::cout << " cells are in ring ";
        std::cout << std::setw(2);
        std::cout << iRing;
        std::cout << std::endl;
    }
    double cell = 0;
    for ( unsigned int i = 0; i < nCellMap.size(); i++ ) {
        cell += nCellMap[i];
    }
    std::cout << std::endl;
    std::cout << " Total: " << std::endl;
    std::cout << " ";
    std::cout << cell;
    std::cout << " cells over ";
    std::cout << NCELL;
    std::cout << " ( ";
    std::cout << ( cell / ((double)NCELL) );
    std::cout << " % )";
    std::cout << std::endl;
    std::cout << std::endl;

} /* End of Mesh::Debug */


/* End of Mesh.cpp */
