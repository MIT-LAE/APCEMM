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

#include "Core/Mesh.hpp"

Mesh::Mesh( )
{

    /* Default Constructor */

    nx = NX;
    ny = NY;

    xlim      = XLIM;
    ylim_up   = YLIM_UP;
    ylim_down = YLIM_DOWN;

    hx_ = 2 * xlim / nx;
    hy_ = ( ylim_up + ylim_down ) / ny;

    /* Cell center x-coordinates */
    for ( unsigned int i = 0; i < nx; i++ ) {
        x_.push_back( (double) 0.0 );
        x_[i] = i * hx_ - xlim + hx_ / 2.0;
        x_e_.push_back( (double) 0.0 );
        x_e_[i] = x_[i] - hx_ / 2.0;
    }
    x_e_.push_back( (double) 0.0 );
    x_e_[nx] = x_e_[nx-1] + hx_;
    
    /* Cell center y-coordinates */
    for ( unsigned int j = 0; j < ny; j++ ) {
        y_.push_back( (double) 0.0 );
        y_[j] = j * hy_ - ylim_down + hy_ / 2.0;
        y_e_.push_back( (double) 0.0 );
        y_e_[j] = y_[j] - hy_ / 2.0;
    }
    y_e_.push_back( (double) 0.0 );
    y_e_[ny] = y_e_[ny-1] + hy_;

    totArea_ = 0;
    for ( unsigned int jNy = 0; jNy < ny; jNy++ ) {
        areas_.push_back( Real_1DVector( nx, 0.0 ) );
        for ( unsigned int iNx = 0; iNx < nx; iNx++ ) {
            areas_[jNy][iNx] = ( y_e_[jNy+1] - y_e_[jNy] ) * \
                              ( x_e_[iNx+1] - x_e_[iNx] );
            /* Comes down to hx_ * hy_, if mesh is cartesian uniform */
            totArea_ += areas_[jNy][iNx];
        }
    }

} /* End of Mesh::Mesh */

Mesh::Mesh( const Mesh &m )
{

    x_ = m.x_;
    y_ = m.y_;
    x_e_ = m.x_e_;
    y_e_ = m.y_e_;
    areas_ = m.areas_;
    xlim = m.xlim;
    ylim_up = m.ylim_up;
    ylim_down = m.ylim_down;
    hx_ = m.hx_;
    hy_ = m.hy_;
    nx = m.nx;
    ny = m.ny;
    nCellMap = m.nCellMap;
    indList = m.indList;
    RingMeshMap = m.RingMeshMap;

} /* End of Mesh::Mesh */

Mesh& Mesh::operator=( const Mesh &m )
{

    if ( &m == this )
        return *this;

    x_ = m.x_;
    y_ = m.y_;
    x_e_ = m.x_e_;
    y_e_ = m.y_e_;
    areas_ = m.areas_;
    xlim = m.xlim;
    ylim_up = m.ylim_up;
    ylim_down = m.ylim_down;
    hx_ = m.hx_;
    hy_ = m.hy_;
    nx = m.nx;
    ny = m.ny;
    nCellMap = m.nCellMap;
    indList = m.indList;
    RingMeshMap = m.RingMeshMap;
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

    /* If we use half-rings, break X_SYMMETRY */
    if ( c.halfRing() ) {

#undef X_SYMMETRY
#define X_SYMMETRY              0

    }

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

    /* For rings and ambient */
    for ( unsigned int iRing = 0; iRing < nRing + 1; iRing++ ) {
        nCellMap.push_back( 0 );
        RingMeshMap.push_back( v2d );
        indList.push_back( v1d_int );
        for ( unsigned int iNy = 0; iNy < NY; iNy++ ) {
            RingMeshMap[iRing].push_back( v1d );
        }
    }

    std::vector<Ring> RingV;
    RingV = c.getRings();
    if ( RingV[nRing - 1].getHAxis() > x_[NX - 1] ) {
        std::cout << "The largest ring's horizontal axis is larger than the grid's dimensions!\n";
        std::cout << "Horizontal axis: " << RingV[nRing-1].getHAxis() << " >= " << x_[NX - 1] << "\n";
    }
    if ( RingV[nRing - 1].getVAxis() > y_[NY - 1] ) {
        std::cout << "The largest ring's vertical axis is larger than the grid's dimensions!\n";
        std::cout << "Vertical axis: " << RingV[nRing-1].getVAxis() << " >= " << y_[NY - 1] << "\n";
    }
    double hAxis, vAxis, hAxis_in, vAxis_in;
    double xRatio, xRatio_in;
    unsigned int val;
    for ( unsigned int iRing = 0; iRing < nRing + 1; iRing++ ) {
        val = 0;
        if ( !c.halfRing() ) {
            if ( iRing == 0 ) {
                hAxis = RingV[iRing].getHAxis();
                vAxis = RingV[iRing].getVAxis();
                /* Use symmetry */
                for ( unsigned int iNx = 0; iNx < nx_max; iNx++ ) {
                    xRatio = ( x_[iNx]  /  hAxis ) * ( x_[iNx]  /  hAxis );
                    for ( unsigned int jNy = 0; jNy < ny_max; jNy++ ) {
                        /* If point is within the region and out of the smaller region */
                        if ( xRatio + ( y_[jNy] / vAxis ) * ( y_[jNy] / vAxis ) <= 1 ) {
                            RingMeshMap[iRing][jNy][iNx] = 1;
                            indList[iRing].push_back( std::make_pair(iNx, jNy) );
                            val++;
                        }
                    }
                }
            }
            else if ( ( iRing > 0 ) && ( iRing < nRing ) ) {
                hAxis = RingV[iRing].getHAxis();
                vAxis = RingV[iRing].getVAxis();
                hAxis_in = RingV[iRing-1].getHAxis();
                vAxis_in = RingV[iRing-1].getVAxis();
                for ( unsigned int iNx = 0; iNx < nx_max; iNx++ ) {
                    xRatio    = ( x_[iNx]  /  hAxis    ) * ( x_[iNx]  /  hAxis    );
                    xRatio_in = ( x_[iNx]  /  hAxis_in ) * ( x_[iNx]  /  hAxis_in );
                    for ( unsigned int jNy = 0; jNy < ny_max; jNy++ ) {
                        /* If point is within the region */
                        if ( ( xRatio + ( y_[jNy] / vAxis ) * ( y_[jNy] / vAxis ) <= 1 ) && ( xRatio_in + ( y_[jNy] / vAxis_in ) * ( y_[jNy] / vAxis_in ) > 1 ) ) {
                            RingMeshMap[iRing][jNy][iNx] = 1;
                            indList[iRing].push_back( std::make_pair(iNx, jNy) );
                            val++;
                        }
                    }
                }
            }
            else {
                hAxis = RingV[nRing-1].getHAxis();
                vAxis = RingV[nRing-1].getVAxis();
                for ( unsigned int iNx = 0; iNx < nx_max; iNx++ ) {
                    xRatio    = ( x_[iNx]  /  hAxis    ) * ( x_[iNx]  /  hAxis    );
                    for ( unsigned int jNy = 0; jNy < ny_max; jNy++ ) {
                        /* If point is within the region and out of the smaller region */
                        if ( ( xRatio + ( y_[jNy] / vAxis ) * ( y_[jNy] / vAxis ) > 1 ) ) {
                            RingMeshMap[iRing][jNy][iNx] = 1;
                            indList[iRing].push_back( std::make_pair(iNx, jNy) );
                            val++;
                        }
                    }
                }
            }
        } else {
            if ( iRing == 0 ) {
                hAxis = RingV[iRing].getHAxis();
                vAxis = RingV[iRing].getVAxis();
                /* Use symmetry */
                for ( unsigned int iNx = 0; iNx < nx_max; iNx++ ) {
                    xRatio = ( x_[iNx]  /  hAxis ) * ( x_[iNx]  /  hAxis );
                    for ( unsigned int jNy = 0; jNy < ny_max; jNy++ ) {
                        /* If point is within the region */
                        if ( ( xRatio + ( y_[jNy] / vAxis ) * ( y_[jNy] / vAxis ) <= 1 ) && ( y_[jNy] >= 0 ) ) {
                            RingMeshMap[iRing][jNy][iNx] = 1;
                            indList[iRing].push_back( std::make_pair(iNx, jNy) );
                            RingMeshMap[iRing+1][(NY - 1) - jNy][iNx] = 1;
                            indList[iRing+1].push_back( std::make_pair(iNx, (NY - 1) - jNy) );
                            val++;
                        }
                    }
                }
                iRing++;
            }
            else if ( ( iRing > 0 ) && ( iRing < nRing ) ) {
                hAxis = RingV[iRing].getHAxis();
                vAxis = RingV[iRing].getVAxis();
                hAxis_in = RingV[iRing-2].getHAxis();
                vAxis_in = RingV[iRing-2].getVAxis();
                for ( unsigned int iNx = 0; iNx < nx_max; iNx++ ) {
                    xRatio    = ( x_[iNx]  /  hAxis    ) * ( x_[iNx]  /  hAxis    );
                    xRatio_in = ( x_[iNx]  /  hAxis_in ) * ( x_[iNx]  /  hAxis_in );
                    for ( unsigned int jNy = 0; jNy < ny_max; jNy++ ) {
                        /* If point is within the region and out of the smaller region */
                        if ( ( ( xRatio + ( y_[jNy] / vAxis ) * ( y_[jNy] / vAxis ) <= 1 ) && ( xRatio_in + ( y_[jNy] / vAxis_in ) * ( y_[jNy] / vAxis_in ) > 1 ) ) && ( y_[jNy] >= 0 ) ) {
                            RingMeshMap[iRing][jNy][iNx] = 1;
                            indList[iRing].push_back( std::make_pair(iNx, jNy) );
                            RingMeshMap[iRing+1][(NY - 1) - jNy][iNx] = 1;
                            indList[iRing+1].push_back( std::make_pair(iNx, (NY - 1) - jNy) );
                            val++;
                        }
                    }
                }
                iRing++;
            }
            else {
                hAxis = RingV[nRing-1].getHAxis();
                vAxis = RingV[nRing-1].getVAxis();
                for ( unsigned int iNx = 0; iNx < nx_max; iNx++ ) {
                    xRatio    = ( x_[iNx]  /  hAxis    ) * ( x_[iNx]  /  hAxis    );
                    for ( unsigned int jNy = 0; jNy < ny_max; jNy++ ) {
                        /* If point is out of the smaller region */
                        if ( ( xRatio + ( y_[jNy] / vAxis ) * ( y_[jNy] / vAxis ) > 1 ) ) {
                            RingMeshMap[iRing][jNy][iNx] = 1;
                            indList[iRing].push_back( std::make_pair(iNx, jNy) );
                            val++;
                        }
                    }
                }
            }
        } 

        if ( val != 0 ) {
            if ( !c.halfRing() ) {
                nCellMap[iRing] = val;
            } else {
                if ( iRing < nRing ) {
                    nCellMap[iRing-1] = val;
                    nCellMap[iRing]   = val;
                } 
                else {
                    nCellMap[iRing] = val;
                }
            }
        }
        else {
            std::cout << "Ring " << iRing << " has no cell in it (nMap = 0)!!" << std::endl;
        }

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
        if ( !c.halfRing() ) {
            for ( unsigned int iNx = nx_max; iNx < NX; iNx++ ) {
                for ( unsigned int jNy = 0; jNy < NY; jNy++ ) {
                    RingMeshMap[iRing][jNy][iNx] = RingMeshMap[iRing][jNy][(NX - 1) - iNx];
                }
            }

            for ( unsigned int iList = 0; iList < nCellMap[iRing]; iList++ ) {
                indList[iRing].push_back( std::make_pair((NX - 1) - indList[iRing][iList].first, indList[iRing][iList].second) );
            }
            nCellMap[iRing] *= 2;
        } else {
            if ( iRing < nRing ) {
                for ( unsigned int iNx = nx_max; iNx < NX; iNx++ ) {
                    for ( unsigned int jNy = 0; jNy < NY; jNy++ ) {
                        RingMeshMap[iRing-1][jNy][iNx] = RingMeshMap[iRing-1][jNy][(NX - 1) - iNx];
                        RingMeshMap[iRing][jNy][iNx] = RingMeshMap[iRing][jNy][(NX - 1) - iNx];
                    }
                }

                for ( unsigned int iList = 0; iList < nCellMap[iRing-1]; iList++ ) {
                    indList[iRing-1].push_back( std::make_pair((NX - 1) - indList[iRing-1][iList].first, indList[iRing-1][iList].second) );
                }
                for ( unsigned int iList = 0; iList < nCellMap[iRing]; iList++ ) {
                    indList[iRing].push_back( std::make_pair((NX - 1) - indList[iRing][iList].first, indList[iRing][iList].second) );
                }
                nCellMap[iRing-1] *= 2;
                nCellMap[iRing] *= 2;
            }
            else {
                for ( unsigned int iNx = nx_max; iNx < NX; iNx++ ) {
                    for ( unsigned int jNy = 0; jNy < NY; jNy++ ) {
                        RingMeshMap[iRing][jNy][iNx] = RingMeshMap[iRing][jNy][(NX - 1) - iNx];
                    }
                }

                for ( unsigned int iList = 0; iList < nCellMap[iRing]; iList++ ) {
                    indList[iRing].push_back( std::make_pair((NX - 1) - indList[iRing][iList].first, indList[iRing][iList].second) );
                }
                nCellMap[iRing] *= 2;
            }
        }

#endif

    }


} /* End of Mesh::Ring2Mesh */

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
    for ( unsigned int i = 0; i < nCellMap.size() - 1; i++ ) {
        cell += nCellMap[i];
    }
    std::cout << std::endl;
    std::cout << " Rings: " << std::endl;
    std::cout << " ";
    std::cout << cell;
    std::cout << " cells over ";
    std::cout << NCELL;
    std::cout << " ( ";
    std::cout << 100 * ( cell / ((double)NCELL) );
    std::cout << " % )";
    std::cout << std::endl;
    std::cout << std::endl;

} /* End of Mesh::Debug */


/* End of Mesh.cpp */
