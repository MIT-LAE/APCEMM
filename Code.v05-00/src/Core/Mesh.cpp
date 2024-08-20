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

#include <iostream>
#include "Util/PhysConstant.hpp"
#include "Core/Mesh.hpp"

Mesh::Mesh(const OptInput& optInput):
    xlim_right( optInput.ADV_GRID_XLIM_RIGHT ),
    xlim_left( optInput.ADV_GRID_XLIM_LEFT ),
    ylim_up( optInput.ADV_GRID_YLIM_UP ),
    ylim_down( optInput.ADV_GRID_YLIM_DOWN ),
    nx( optInput.ADV_GRID_NX ),
    ny( optInput.ADV_GRID_NY )
{
    initCoordVectors(MeshDomainLimitsSpec::CENTERED_LIMITS);
}
Mesh::Mesh(int nx, int ny, double xlim_left, double xlim_right, double ylim_up, double ylim_down, MeshDomainLimitsSpec limitsSpec):
    xlim_right( xlim_right ),
    xlim_left( xlim_left ),
    ylim_up( ylim_up ),
    ylim_down( ylim_down ),
    nx( nx ),
    ny( ny )
{
    initCoordVectors(limitsSpec);

} /* End of Mesh::Mesh */

void Mesh::initCoordVectors(MeshDomainLimitsSpec limitsSpec){

    hx_ = limitsSpec == MeshDomainLimitsSpec::CENTERED_LIMITS
                    ? ( xlim_right + xlim_left ) / nx
                    : ( xlim_right - xlim_left ) / nx;
    hy_ = limitsSpec == MeshDomainLimitsSpec::CENTERED_LIMITS
                    ? ( ylim_up + ylim_down ) / ny 
                    : ( ylim_up - ylim_down ) / ny;

    /* Cell center x-coordinates */
    for ( UInt i = 0; i < nx; i++ ) {
        x_.push_back( i * hx_ - xlim_left + hx_ / 2.0 );
        x_e_.push_back( x_[i] - hx_ / 2.0 );
        if ( i > 0 )
            dx_.push_back( x_e_[i] - x_e_[i-1] );
    }
    x_e_.push_back( x_e_[nx-1] + hx_ );
    dx_.push_back( x_e_[nx] - x_e_[nx-1] );
    
    /* Cell center y-coordinates */
    for ( UInt j = 0; j < ny; j++ ) {
        y_.push_back( j * hy_ - ylim_down + hy_ / 2.0 );
        y_e_.push_back( y_[j] - hy_ / 2.0 );
        if ( j > 0 )
            dy_.push_back( y_e_[j] - y_e_[j-1] );
    }
    y_e_.push_back( y_e_[ny-1] + hy_ );
    dy_.push_back( y_e_[ny] - y_e_[ny-1] );
    
    areas_ = Vector_2D(ny, Vector_1D(nx, 0));
    calcAreas();
    
}
void Mesh::Ring2Mesh( Cluster &c )
{

    UInt maxRing = 0;
    UInt nRing   = c.getnRing();

    UInt nx_max, ny_max;

    nx_max = nx;
    ny_max = ny;

    std::vector<Ring> RingV;
    RingV = c.getRings();

    double hAxis, vAxis, hAxis_in, vAxis_in;
    double xRatio, xRatio_in;
    UInt nCell;

    /* In this scenario, chemistry and microphysics are performed at the
     * grid-scale level.
     * We thus only care about the most inner ring to release emissions */

    if ( c.halfRing() )
        maxRing = 2;
    else
        maxRing = 1;

    Vector_2D v2d;
    Vector_1D v1d( nx, 0.0E+00 );

    /* For rings and ambient */
    for ( UInt iRing = 0; iRing < maxRing; iRing++ ) {
        nCellMap.push_back( 0 );
        weights.push_back( v2d );
        for ( UInt iNy = 0; iNy < ny; iNy++ ) {
            weights[iRing].push_back( v1d );
        }
    }
    for ( UInt iNy = 0; iNy < ny; iNy++ )
        mapIndex_.push_back( Vector_1Dui( nx, 0 ) );


    for ( UInt iRing = 0; iRing < maxRing; iRing++ ) {
        nCell = 0;
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
                            nCell++;
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
                            nCell++;
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
                            nCell++;
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
                            nCell++;
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
                            nCell++;
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
                                nCell++;
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
                                nCell++;
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
                            nCell++;
                        }
                    }
                }
            }
        } 

        if ( nCell != 0 ) {
            nCellMap[iRing] = nCell;
            /* Map should be such that \iint w dA = A_ring
             * <=> \sum w_ij A_ij = A_ring 
             * For a uniform mesh, this comes down to:
             * \sum w_ij = A_ring / A_ij = # of grid cell in ring */
        }
        else {
            std::cout << "Ring " << iRing << " has no cell in it (nCell = 0)!!" << std::endl;
            if ( iRing == 0 ) {
                /* Just add to central ring */
                UInt jNy = 0;
                UInt iNx = 0;
                double ringArea;

                /* Calculate ring Area */
                hAxis = RingV[iRing].getHAxis();
                vAxis = RingV[iRing].getVAxis();
                ringArea = physConst::PI * hAxis * vAxis;
                std::cout << ringArea << std::endl;

                /* Split over 4 cells */
                /* First cell */
                jNy = std::floor( y_[0]/(y_[0]-y_[1]) ); //std::ceil(ny/2);
                iNx = std::floor( x_[0]/(x_[0]-x_[1]) ); //std::ceil(nx/2);
                if ( y_[jNy] < 0 ) {
                    jNy+=1;
                }
                if ( x_[iNx] < 0 ) {
                    iNx+=1;
                }
                //First cell
                weights[iRing][jNy][iNx] = 0.25 * ringArea / cellArea_; //1.0E+00;
                mapIndex_[jNy][iNx] = iRing;

                /* Second cell */
                weights[iRing][jNy-1][iNx] = 0.25 * ringArea / cellArea_; //1.0E+00;
		        mapIndex_[jNy-1][iNx] = iRing;

                /* Third cell */
                weights[iRing][jNy][iNx-1] = 0.25 * ringArea / cellArea_; //1.0E+00;
                mapIndex_[jNy][iNx-1] = iRing;

                /* Fourth cell */
                weights[iRing][jNy-1][iNx-1] = 0.25 * ringArea / cellArea_; //1.0E+00;
                mapIndex_[jNy-1][iNx-1] = iRing;

                /* nCellMap */
                nCellMap[iRing] = 4;

            }
            else {
                exit(-1);
            }
        }
    }


} /* End of Mesh::Ring2Mesh */

void Mesh::MapWeights( )
{

    UInt jNy, iNx, iRing;
    UInt nRing = weights.size(); // This is actually equal to nRing+1

    double max = 0.0E+00;
    for ( jNy = 0; jNy < ny; jNy++ ) {
        for ( iNx = 0; iNx < nx; iNx++ ) {
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

void Mesh::updateVertGrid( const Vector_1D& yE_new ) {
    for(std::size_t i = 0; i < yE_new.size(); i++) {
        y_[i] = (yE_new[i] + yE_new[i+1]) / 2;
        dy_[i] = yE_new[i+1] - yE_new[i];
    }
    calcAreas();
    y_e_ = yE_new;
}

/* End of Mesh.cpp */
