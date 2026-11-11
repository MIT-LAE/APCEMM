/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*                        AIrcraft Microphysics                     */
/*                              (AIM)                               */
/*                                                                  */
/* Coagulation Program File                                         */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 9/27/2018                                 */
/* File                 : Coagulation.cpp                           */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "AIM/buildKernel.hpp"
#include "Util/PhysFunction.hpp"
#include "AIM/Coagulation.hpp"


namespace AIM
{

    Coagulation::Coagulation( CoagulationPair pair, Vector_1D const &bin_Centers_1, Vector_1D const &bin_VCenters_1, double rho_1, double temperature_K_, double pressure_Pa_ ): Kernel( bin_Centers_1.size() )
    {

        /* Constructor */

        double temperature_K = temperature_K_;
        double pressure_Pa = pressure_Pa_;

        /* Set shape */
        for ( unsigned int iBin_1 = 0; iBin_1 < bin_Centers_1.size(); iBin_1++ ) {
            Kernel[iBin_1].resize( bin_Centers_1.size() );
        }

        if ( ( pair == CoagulationPair::LiquidLiquid ) || ( pair == CoagulationPair::LiquidSoot ) ) {
            /* Liquid-liquid and liquid-soot coagulation */

            /* Declare and assign coagulation kernels for different physical processes */
            Vector_2D K_Brow = buildBrownianKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_1, rho_1 );
            Vector_2D K_DE   = buildDEKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_1, rho_1, K_Brow );

            /* Fill in total coagulation kernel */
            for ( unsigned int iBin_1 = 0; iBin_1 < bin_Centers_1.size(); iBin_1++ ) {
                for ( unsigned int iBin_2 = 0; iBin_2 < bin_Centers_1.size(); iBin_2++ ) {
                    Kernel[iBin_1][iBin_2] = ( K_Brow[iBin_1][iBin_2] + K_DE[iBin_1][iBin_2] ) * 1.00E+06; /* Conversion from m^3/s to cm^3/s */
                }
            }

        }
        else if ( pair == CoagulationPair::IceIce ) {
            /* Ice-Ice coagulation */

            /* Declare and assign coagulation kernels for different physical processes */
            Vector_2D K_Brow = buildBrownianKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_1, rho_1 );
            Vector_2D K_DE   = buildDEKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_1, rho_1, K_Brow );
            Vector_2D K_GC   = buildGCKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_1, rho_1 );
            Vector_2D K_TI   = buildTIKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_1, rho_1 );
            Vector_2D K_TS   = buildTSKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_1, rho_1 );
            
            /* Fill in total coagulation kernel */
            for ( unsigned int iBin_1 = 0; iBin_1 < bin_Centers_1.size(); iBin_1++ ) {
                for ( unsigned int iBin_2 = 0; iBin_2 < bin_Centers_1.size(); iBin_2++ ) {
                    Kernel[iBin_1][iBin_2] = ( K_Brow[iBin_1][iBin_2] + K_DE[iBin_1][iBin_2] + K_GC[iBin_1][iBin_2] + K_TI[iBin_1][iBin_2] + K_TS[iBin_1][iBin_2] ) * 1.00E+06; /* Conversion from m^3/s to cm^3/s */
                }
            }

        }
        else {
            throw std::invalid_argument( "In AIM::Coagulation::Coagulation: invalid CoagulationPair value." );
        }

        for ( unsigned int iBin_1 = 0; iBin_1 < bin_Centers_1.size(); iBin_1++ ) {
            beta.push_back( Vector_1D( bin_Centers_1.size() ) );
            for ( unsigned int iBin_2 = 0; iBin_2 < bin_Centers_1.size(); iBin_2++ ) {
                /* Assuming an aggregation efficiency of 1 */
                beta[iBin_1][iBin_2] = Kernel[iBin_1][iBin_2];
            }
        }
        buildF   ( bin_VCenters_1 );

    } /* End of Coagulation::Coagulation */
            
    Coagulation::Coagulation( CoagulationPair pair, Vector_1D const &bin_Centers_1, double rho_1, double radius_2, double rho_2, double temperature_K_, double pressure_Pa_ ):
        Kernel_1D( bin_Centers_1.size() )
    {

        /* Constructor */

        double temperature_K = temperature_K_;
        double pressure_Pa = pressure_Pa_;

        /* Declare coagulation kernel and set shape */

        if ( ( pair == CoagulationPair::LiquidLiquid ) || ( pair == CoagulationPair::LiquidSoot ) ) {
            /* Liquid-liquid and liquid-soot coagulation */

            /* Declare and assign coagulation kernels for different physical processes */
            Vector_1D K_Brow = buildBrownianKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, radius_2, rho_2 );
            Vector_1D K_DE   = buildDEKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, radius_2, rho_2, K_Brow );

            /* Fill in total coagulation kernel */
            for ( unsigned int iBin_1 = 0; iBin_1 < bin_Centers_1.size(); iBin_1++ ) {
                Kernel_1D[iBin_1] = ( K_Brow[iBin_1] + K_DE[iBin_1] ) * 1.00E+06; /* Conversion from m^3/s to cm^3/s */
            }

        }
        else if ( pair == CoagulationPair::IceIce ) {
            /* Ice-Ice coagulation */

            /* Declare and assign coagulation kernels for different physical processes */
            Vector_1D K_Brow = buildBrownianKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, radius_2, rho_2 );
            Vector_1D K_DE   = buildDEKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, radius_2, rho_2, K_Brow );
            Vector_1D K_GC   = buildGCKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, radius_2, rho_2 );
            Vector_1D K_TI   = buildTIKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, radius_2, rho_2 );
            Vector_1D K_TS   = buildTSKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, radius_2, rho_2 );
            
            /* Fill in total coagulation kernel */
            for ( unsigned int iBin_1 = 0; iBin_1 < bin_Centers_1.size(); iBin_1++ ) {
                Kernel_1D[iBin_1] = ( K_Brow[iBin_1] + K_DE[iBin_1] + K_GC[iBin_1] + K_TI[iBin_1] + K_TS[iBin_1] ) * 1.00E+06; /* Conversion from m^3/s to cm^3/s */
            }

        }
        else {
            throw std::invalid_argument( "In AIM::Coagulation::Coagulation: invalid CoagulationPair value." );
        }

    } /* End of Coagulation::Coagulation */

    void Coagulation::buildF(const Vector_1D &bin_VCenters )
    {

        double vij;
        UInt iBin, jBin, index;
        UInt size = bin_VCenters.size();

        Vector_2D v2d;
        std::vector<std::vector<UInt> > v2d_index;
        std::vector<UInt> v1d_index;

        for ( iBin = 0; iBin < size; iBin++ ) {
            indices.push_back( v2d_index );
            for ( jBin = 0; jBin < size; jBin++ )
                indices[iBin].push_back( v1d_index );
        }

        for ( iBin = 0; iBin < size; iBin++ ) {
            f.push_back( v2d );
            for ( jBin = 0; jBin < size; jBin++ ) {
                f[iBin].push_back( Vector_1D( size, 0.0 ) );
                vij = bin_VCenters[iBin] + bin_VCenters[jBin];
                /* Using std functions: */
                index = std::distance( bin_VCenters.begin(), std::upper_bound( bin_VCenters.begin(), bin_VCenters.end(), vij ) ) - 1;

                if ( index < size-1 ) {
                    f[iBin][jBin][index  ] = ( bin_VCenters[index+1] - vij ) / ( bin_VCenters[index+1] - bin_VCenters[index] ) * bin_VCenters[index] / vij;
                    f[iBin][jBin][index+1] = 1.0 - f[iBin][jBin][index];
                    indices[jBin][index  ].push_back( iBin );
                    indices[jBin][index+1].push_back( iBin );
                } else if ( index >= size-1 ) {
                    index = size-1;
                    f[iBin][jBin][index  ] = 1.0;
                    indices[jBin][index  ].push_back( iBin );
                }

                /* For debug purposes */
                if ( 0 ) {
                    std::cout << "(iBin, jBin) = (" << iBin << ", " << jBin << ") , Index: " << index << ", ";
                    std::cout << bin_VCenters[index] << " <= " << vij << " < " << bin_VCenters[index+1] << "\n";

                    std::cout << "(iBin, jBin) = (" << iBin << ", " << jBin << ") , Index    : " << index << ", f = " << f[iBin][jBin][index] << "\n";
                    std::cout << "(iBin, jBin) = (" << iBin << ", " << jBin << ") , Index + 1: " << index+1 << ", f = " << f[iBin][jBin][index+1] << "\n";
                }

            }
        }

    } /* End of Coagulation::buildF */

    void Coagulation::buildF( const Vector_3D &bin_VCenters, const UInt jNy, const UInt iNx )
    {

        double vij;
        UInt iBin, jBin, kBin, index;
        UInt size = bin_VCenters.size();
        Vector_1D bin_VCentersCopy( size, 0.0E+00 );

        Vector_2D v2d;
        std::vector<std::vector<UInt> > v2d_index;
        std::vector<UInt> v1d_index;

        for ( iBin = 0; iBin < size; iBin++ ) {
            bin_VCentersCopy[iBin] = bin_VCenters[iBin][jNy][iNx];
            for ( jBin = 0; jBin < size; jBin++ )
                indices[iBin][jBin].clear();
        }

        for ( iBin = 0; iBin < size; iBin++ ) {
            for ( jBin = 0; jBin < size; jBin++ ) {
                for ( kBin = 0; kBin < size; kBin++ )
                    f[iBin][jBin][kBin] = 0.0E+00;
                vij = bin_VCentersCopy[iBin] + bin_VCentersCopy[jBin];
                /* Using std functions: */
                index = std::distance( bin_VCentersCopy.begin(), std::upper_bound( bin_VCentersCopy.begin(), bin_VCentersCopy.end(), vij ) ) - 1;
                if ( index < size-1 ) {
                    f[iBin][jBin][index  ] = ( bin_VCentersCopy[index+1] - vij ) / ( bin_VCentersCopy[index+1] - bin_VCentersCopy[index] ) * bin_VCentersCopy[index] / vij;
                    f[iBin][jBin][index+1] = 1.0 - f[iBin][jBin][index];
                    indices[jBin][index  ].push_back( iBin );
                    indices[jBin][index+1].push_back( iBin );
                } else if ( index >= size-1 ) {
                    index = size-1;
                    f[iBin][jBin][index  ] = 1.0;
                    indices[jBin][index  ].push_back( iBin );
                }

                /* For debug purposes */
                if ( 0 ) {
                    std::cout << "(iBin, jBin) = (" << iBin << ", " << jBin << ") , Index: " << index << ", ";
                    std::cout << bin_VCentersCopy[index] << " <= " << vij << " < " << bin_VCentersCopy[index+1] << "\n";

                    std::cout << "(iBin, jBin) = (" << iBin << ", " << jBin << ") , Index    : " << index << ", f = " << f[iBin][jBin][index] << "\n";
                    std::cout << "(iBin, jBin) = (" << iBin << ", " << jBin << ") , Index + 1: " << index+1 << ", f = " << f[iBin][jBin][index+1] << "\n";
                }

            }
        }

    } /* End of Coagulation::buildF */

    Vector_2D Coagulation::getKernel() const
    {

        return Kernel;

    } /* End of Coagulation::getKernel */
            
    Vector_1D Coagulation::getKernel_1D() const
    {

        return Kernel_1D;

    } /* End of Coagulation::getKernel_1D */

    Vector_2D Coagulation::getBeta() const
    {

        return beta;

    } /* End of Coagulation::getBeta */

    Vector_3D Coagulation::getF() const
    {

        return f;

    } /* End of Coagulation::getF */

}

/* End of Coagulation.cpp */
