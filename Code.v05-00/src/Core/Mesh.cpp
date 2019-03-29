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
    for ( UInt i = 0; i < nx; i++ ) {
        x_.push_back( (RealDouble) 0.0 );
        x_[i] = i * hx_ - xlim + hx_ / 2.0;
        x_e_.push_back( (RealDouble) 0.0 );
        x_e_[i] = x_[i] - hx_ / 2.0;
    }
    x_e_.push_back( (RealDouble) 0.0 );
    x_e_[nx] = x_e_[nx-1] + hx_;
    
    /* Cell center y-coordinates */
    for ( UInt j = 0; j < ny; j++ ) {
        y_.push_back( (RealDouble) 0.0 );
        y_[j] = j * hy_ - ylim_down + hy_ / 2.0;
        y_e_.push_back( (RealDouble) 0.0 );
        y_e_[j] = y_[j] - hy_ / 2.0;
    }
    y_e_.push_back( (RealDouble) 0.0 );
    y_e_[ny] = y_e_[ny-1] + hy_;

    totArea_ = 0;
    /* TO CHANGE if mesh is non-uniform */
    cellArea_ = ( y_e_[1] - y_e_[0] ) * ( x_e_[1] - x_e_[0] );
    totArea_  = ny * nx * cellArea_;
    for ( UInt jNy = 0; jNy < ny; jNy++ ) {
        areas_.push_back( Vector_1D( nx, 0.0 ) );
        for ( UInt iNx = 0; iNx < nx; iNx++ ) {
            areas_[jNy][iNx] = cellArea_;
            //areas_[jNy][iNx] = ( y_e_[jNy+1] - y_e_[jNy] ) * \
            //                  ( x_e_[iNx+1] - x_e_[iNx] );
            ///* Comes down to hx_ * hy_, if mesh is cartesian uniform */
            //totArea_ += areas_[jNy][iNx];
        }
    }

} /* End of Mesh::Mesh */

Mesh::Mesh( const Mesh &m )
{

    x_          = m.x_;
    y_          = m.y_;
    x_e_        = m.x_e_;
    y_e_        = m.y_e_;
    areas_      = m.areas_;
    totArea_    = m.totArea_;
    cellArea_   = m.cellArea_;
    xlim        = m.xlim;
    ylim_up     = m.ylim_up;
    ylim_down   = m.ylim_down;
    hx_         = m.hx_;
    hy_         = m.hy_;
    nx          = m.nx;
    ny          = m.ny;
    nCellMap    = m.nCellMap;
    weights     = m.weights;
    mapIndex_   = m.mapIndex_;

} /* End of Mesh::Mesh */

Mesh& Mesh::operator=( const Mesh &m )
{

    if ( &m == this )
        return *this;

    x_          = m.x_;
    y_          = m.y_;
    x_e_        = m.x_e_;
    y_e_        = m.y_e_;
    areas_      = m.areas_;
    totArea_    = m.totArea_;
    cellArea_   = m.cellArea_;
    xlim        = m.xlim;
    ylim_up     = m.ylim_up;
    ylim_down   = m.ylim_down;
    hx_         = m.hx_;
    hy_         = m.hy_;
    nx          = m.nx;
    ny          = m.ny;
    nCellMap    = m.nCellMap;
    weights     = m.weights;
    mapIndex_   = m.mapIndex_;

    return *this;

} /* End of Mesh::operator= */

Mesh::~Mesh( )
{

    /* Destructor */

} /* End of Mesh::~Mesh */

void Mesh::Ring2Mesh( Cluster &c )
{
    
    UInt nRing = c.getnRing();

    Vector_2D v2d;
    Vector_1D v1d( NX, 0.0E+00 );

    UInt nx_max, ny_max;

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
    for ( UInt iRing = 0; iRing < nRing + 1; iRing++ ) {
        nCellMap.push_back( 0 );
        weights.push_back( v2d );
        for ( UInt iNy = 0; iNy < NY; iNy++ ) {
            weights[iRing].push_back( v1d );
        }
    }
    for ( UInt iNy = 0; iNy < NY; iNy++ )
        mapIndex_.push_back( Vector_1Dui( NX, 0 ) );

    std::vector<Ring> RingV;
    RingV = c.getRings();
    if ( RingV[nRing - 1].getHAxis() > x_[NX - 1] ) {
        std::cout << "The largest ring's horizontal axis is larger than the grid's dimensions!\n";
        std::cout << "Horizontal axis: " << RingV[nRing-1].getHAxis() << " >= " << x_[NX - 1] << std::endl;
        exit(-1);
    }
    if ( RingV[nRing - 1].getVAxis() > y_[NY - 1] ) {
        std::cout << "The largest ring's vertical axis is larger than the grid's dimensions!\n";
        std::cout << "Vertical axis: " << RingV[nRing-1].getVAxis() << " >= " << y_[NY - 1] << std::endl;
        exit(-1);
    }

    RealDouble hAxis, vAxis, hAxis_in, vAxis_in;
    RealDouble xRatio, xRatio_in;
    UInt val;

    for ( UInt iRing = 0; iRing < nRing + 1; iRing++ ) {
        val = 0;
        if ( !c.halfRing() ) {
            /* Full rings */
            if ( iRing == 0 ) {
                hAxis = RingV[iRing].getHAxis();
                vAxis = RingV[iRing].getVAxis();
                /* Use symmetry */
                for ( UInt iNx = 0; iNx < nx_max; iNx++ ) {
                    xRatio = ( x_[iNx]  /  hAxis ) * ( x_[iNx]  /  hAxis );
                    for ( UInt jNy = 0; jNy < ny_max; jNy++ ) {
                        /* If point is within the region and out of the smaller region */
                        if ( xRatio + ( y_[jNy] / vAxis ) * ( y_[jNy] / vAxis ) <= 1 ) {
                            weights[iRing][jNy][iNx] = 1.0E+00;
                            mapIndex_[jNy][iNx] = iRing;
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
                for ( UInt iNx = 0; iNx < nx_max; iNx++ ) {
                    xRatio    = ( x_[iNx]  /  hAxis    ) * ( x_[iNx]  /  hAxis    );
                    xRatio_in = ( x_[iNx]  /  hAxis_in ) * ( x_[iNx]  /  hAxis_in );
                    for ( UInt jNy = 0; jNy < ny_max; jNy++ ) {
                        /* If point is within the region */
                        if ( ( xRatio + ( y_[jNy] / vAxis ) * ( y_[jNy] / vAxis ) <= 1 ) && ( xRatio_in + ( y_[jNy] / vAxis_in ) * ( y_[jNy] / vAxis_in ) > 1 ) ) {
                            weights[iRing][jNy][iNx] = 1.0E+00;
                            mapIndex_[jNy][iNx] = iRing;
                            val++;
                        }
                    }
                }
            }
            else {
                hAxis = RingV[nRing-1].getHAxis();
                vAxis = RingV[nRing-1].getVAxis();
                for ( UInt iNx = 0; iNx < nx_max; iNx++ ) {
                    xRatio    = ( x_[iNx]  /  hAxis    ) * ( x_[iNx]  /  hAxis    );
                    for ( UInt jNy = 0; jNy < ny_max; jNy++ ) {
                        /* If point is within the region and out of the smaller region */
                        if ( ( xRatio + ( y_[jNy] / vAxis ) * ( y_[jNy] / vAxis ) > 1 ) ) {
                            weights[iRing][jNy][iNx] = 1.0E+00;
                            mapIndex_[jNy][iNx] = iRing;
                            val++;
                        }
                    }
                }
            }
        } else {
            /* Half rings */
            if ( iRing == 0 ) {
                hAxis = RingV[iRing].getHAxis();
                vAxis = RingV[iRing].getVAxis();
                /* Use symmetry */
                for ( UInt iNx = 0; iNx < nx_max; iNx++ ) {
                    xRatio = ( x_[iNx]  /  hAxis ) * ( x_[iNx]  /  hAxis );
                    for ( UInt jNy = 0; jNy < ny_max; jNy++ ) {
                        /* If point is within the region */
                        if ( ( xRatio + ( y_[jNy] / vAxis ) * ( y_[jNy] / vAxis ) <= 1 ) && ( y_[jNy] >= 0 ) ) {
                            weights[iRing][jNy][iNx] = 1.0E+00;
                            mapIndex_[jNy][iNx] = iRing;
                            val++;
                        }
                    }
                }
            } else if ( iRing == 1 ) {
                hAxis = RingV[iRing].getHAxis();
                vAxis = RingV[iRing].getVAxis();
                /* Use symmetry */
                for ( UInt iNx = 0; iNx < nx_max; iNx++ ) {
                    xRatio = ( x_[iNx]  /  hAxis ) * ( x_[iNx]  /  hAxis );
                    for ( UInt jNy = 0; jNy < ny_max; jNy++ ) {
                        /* If point is within the region */
                        if ( ( xRatio + ( y_[jNy] / vAxis ) * ( y_[jNy] / vAxis ) <= 1 ) && ( y_[jNy] < 0 ) ) {
                            weights[iRing][jNy][iNx] = 1.0E+00;
                            mapIndex_[jNy][iNx] = iRing;
                            val++;
                        }
                    }
                }
            } else if ( ( iRing > 1 ) && ( iRing < nRing ) ) {
                if ( iRing % 2 ) {
                    /* Upper rings */
                    hAxis = RingV[iRing].getHAxis();
                    vAxis = RingV[iRing].getVAxis();
                    hAxis_in = RingV[iRing-2].getHAxis();
                    vAxis_in = RingV[iRing-2].getVAxis();
                    for ( UInt iNx = 0; iNx < nx_max; iNx++ ) {
                        xRatio    = ( x_[iNx]  /  hAxis    ) * ( x_[iNx]  /  hAxis    );
                        xRatio_in = ( x_[iNx]  /  hAxis_in ) * ( x_[iNx]  /  hAxis_in );
                        for ( UInt jNy = 0; jNy < ny_max; jNy++ ) {
                            /* If point is within the region and out of the smaller region */
                            if ( ( ( xRatio + ( y_[jNy] / vAxis ) * ( y_[jNy] / vAxis ) <= 1 ) && ( xRatio_in + ( y_[jNy] / vAxis_in ) * ( y_[jNy] / vAxis_in ) > 1 ) ) && ( y_[jNy] >= 0 ) ) {
                                weights[iRing][jNy][iNx] = 1.0E+00;
                                mapIndex_[jNy][iNx] = iRing;
                                val++;
                            }
                        }
                    }
                } else {
                    /* Lower rings */
                    hAxis = RingV[iRing].getHAxis();
                    vAxis = RingV[iRing].getVAxis();
                    hAxis_in = RingV[iRing-2].getHAxis();
                    vAxis_in = RingV[iRing-2].getVAxis();
                    for ( UInt iNx = 0; iNx < nx_max; iNx++ ) {
                        xRatio    = ( x_[iNx]  /  hAxis    ) * ( x_[iNx]  /  hAxis    );
                        xRatio_in = ( x_[iNx]  /  hAxis_in ) * ( x_[iNx]  /  hAxis_in );
                        for ( UInt jNy = 0; jNy < ny_max; jNy++ ) {
                            /* If point is within the region and out of the smaller region */
                            if ( ( ( xRatio + ( y_[jNy] / vAxis ) * ( y_[jNy] / vAxis ) <= 1 ) && ( xRatio_in + ( y_[jNy] / vAxis_in ) * ( y_[jNy] / vAxis_in ) > 1 ) ) && ( y_[jNy] < 0 ) ) {
                                weights[iRing][jNy][iNx] = 1;
                                mapIndex_[jNy][iNx] = iRing;
                                val++;
                            }
                        }
                    }
                }
            }
            else {
                hAxis = RingV[nRing-1].getHAxis();
                vAxis = RingV[nRing-1].getVAxis();
                for ( UInt iNx = 0; iNx < nx_max; iNx++ ) {
                    xRatio    = ( x_[iNx]  /  hAxis    ) * ( x_[iNx]  /  hAxis    );
                    for ( UInt jNy = 0; jNy < ny_max; jNy++ ) {
                        /* If point is out of the smaller region */
                        if ( ( xRatio + ( y_[jNy] / vAxis ) * ( y_[jNy] / vAxis ) > 1 ) ) {
                            weights[iRing][jNy][iNx] = 1.0E+00;
                            mapIndex_[jNy][iNx] = iRing;
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
            /* Map should be such that \iint w dA = A_ring
             * <=> \sum w_ij A_ij = A_ring 
             * For a uniform mesh, this comes down to:
             * \sum w_ij = A_ring / A_ij = # of grid cell in ring */
        }
        else {
            std::cout << "Ring " << iRing << " has no cell in it (nMap = 0)!!" << std::endl;
            exit(-1);
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
        for ( UInt iNx = 0; iNx < nx_max; iNx++ ) {
            for ( UInt jNy = ny_max; jNy < NY; jNy++ ) {
                weights[iRing][jNy][iNx] = weights[iRing][(NY - 1) - jNy][iNx];
                mapIndex_[jNy][iNx] = mapIndex_[(NY - 1) - jNy][iNx];
            }
        }
        
        nCellMap[iRing] *= 2;
        
        /* Do regions 3 and 4 */
        /* 3 */
        for ( UInt iNx = nx_max; iNx < NX; iNx++ ) {
            for ( UInt jNy = 0; jNy < ny_max; jNy++ ) {
                weights[iRing][jNy][iNx] = weights[iRing][jNy][(NX - 1) - iNx];
                mapIndex_[jNy][iNx] = mapIndex_[jNy][(NX - 1) - iNx];
            }
        }
       
        /* 4 */
        for ( UInt iNx = nx_max; iNx < NX; iNx++ ) {
            for ( UInt jNy = ny_max; jNy < NY; jNy++ ) {
                weights[iRing][jNy][iNx] = weights[iRing][(NY - 1) - jNy][(NX - 1) - iNx];
                mapIndex_[jNy][iNx] = mapIndex_[(NY - 1) - jNy][(NX - 1 ) - iNx];
            }
        }
        
        nCellMap[iRing] *= 2;

#endif

#if ( X_SYMMETRY && !Y_SYMMETRY )
        /* Do regions 2 and 4 */
        /* 2 and 4 */
        for ( UInt iNx = 0; iNx < NX; iNx++ ) {
            for ( UInt jNy = ny_max; jNy < NY; jNy++ ) {
                weights[iRing][jNy][iNx] = weights[iRing][(NY - 1) - jNy][iNx];
                mapIndex_[jNy][iNx] = mapIndex_[(NY - 1) - jNy][iNx];
            }
        }
        
        nCellMap[iRing] *= 2;

#endif

#if ( !X_SYMMETRY && Y_SYMMETRY )
        /* Do regions 3 and 4 */
        /* 3 and 4 */
        if ( !c.halfRing() ) {
            for ( UInt iNx = nx_max; iNx < NX; iNx++ ) {
                for ( UInt jNy = 0; jNy < NY; jNy++ ) {
                    weights[iRing][jNy][iNx] = weights[iRing][jNy][(NX - 1) - iNx];
                    mapIndex_[jNy][iNx] = mapIndex_[jNy][(NX - 1) - iNx];
                }
            }

            nCellMap[iRing] *= 2;
        } else {
            if ( iRing < nRing ) {
                for ( UInt iNx = nx_max; iNx < NX; iNx++ ) {
                    for ( UInt jNy = 0; jNy < NY; jNy++ ) {
                        weights[iRing-1][jNy][iNx] = weights[iRing-1][jNy][(NX - 1) - iNx];
                        weights[iRing][jNy][iNx]   = weights[iRing][jNy][(NX - 1) - iNx];
                        mapIndex_[jNy][iNx] = mapIndex_[jNy][(NX - 1) - iNx];
                    }
                }

                nCellMap[iRing-1] *= 2;
                nCellMap[iRing] *= 2;
            }
            else {
                for ( UInt iNx = nx_max; iNx < NX; iNx++ ) {
                    for ( UInt jNy = 0; jNy < NY; jNy++ ) {
                        weights[iRing][jNy][iNx] = weights[iRing][jNy][(NX - 1) - iNx];
                        mapIndex_[jNy][iNx] = mapIndex_[jNy][(NX - 1) - iNx];
                    }
                }

                nCellMap[iRing] *= 2;
            }
        }

#endif

    }


} /* End of Mesh::Ring2Mesh */

void Mesh::MapWeights( )
{

    UInt jNy, iNx, iRing;
    UInt nRing = weights.size(); // This is actually equal to nRing+1

    RealDouble max = 0.0E+00;
    for ( jNy = 0; jNy < NY; jNy++ ) {
        for ( iNx = 0; iNx < NX; iNx++ ) {
            max = weights[0][jNy][iNx];
            mapIndex_[jNy][iNx] = 0;
            for ( iRing = 1; iRing < nRing; iRing++ ) {
                if ( weights[iRing][jNy][iNx] > max ) {
                    max = weights[iRing][jNy][iNx];
                    mapIndex_[jNy][iNx] = iRing;
                }
            }
        }
    }

} /* End of Mesh::MapWeights */

void Mesh::Debug( ) const
{

    std::cout << std::endl;
    std::cout << "**** Input Debugger ****" << std::endl;
    std::cout << "Mesh-Ring structure: " << std::endl;
    std::cout << std::endl;

    std::cout << " Number of cells in each ring: " << std::endl;
    for ( UInt iRing = 0; iRing < nCellMap.size(); iRing++ ) {
        std::cout << " ";
        std::cout << std::setw(4);
        std::cout << nCellMap[iRing];
        std::cout << " cells are in ring ";
        std::cout << std::setw(2);
        std::cout << iRing;
        std::cout << std::endl;
    }
    RealDouble cell = 0;
    for ( UInt i = 0; i < nCellMap.size() - 1; i++ ) {
        cell += nCellMap[i];
    }
    std::cout << std::endl;
    std::cout << " Rings: " << std::endl;
    std::cout << " ";
    std::cout << cell;
    std::cout << " cells over ";
    std::cout << NCELL;
    std::cout << " ( ";
    std::cout << 100 * ( cell / ((RealDouble)NCELL) );
    std::cout << " % )";
    std::cout << std::endl;
    std::cout << std::endl;

} /* End of Mesh::Debug */


/* End of Mesh.cpp */
