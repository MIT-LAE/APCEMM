/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*                        AIrcraft Microphysics                     */
/*                              (AIM)                               */
/*                                                                  */
/* buildKernel Program File                                         */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 9/27/2018                                 */
/* File                 : buildKernel.cpp                           */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <vector>

#include "../../Headers/PhysConstant.hpp"
#include "../../Headers/PhysFunction.hpp"

typedef double RealDouble;
typedef std::vector<RealDouble> Vector_1D;
typedef std::vector<Vector_1D> Vector_2D;

namespace AIM
{

    Vector_1D buildBrownianKernel( RealDouble temperature_K, RealDouble pressure_Pa, Vector_1D bin_Centers, RealDouble rho_1 , RealDouble bin_R, RealDouble rho_2 )
    {

        /* Returns the 1D coagulation kernel associated with Brownian coagulation */

        /* Declare coagulation kernel */
        Vector_1D K_Brownian( bin_Centers.size() );

        /* Declare particle diffusion coefficients */
        Vector_1D  diffCoef_1( bin_Centers.size() );
        RealDouble diffCoef_2;
       
        /* Declare mean distance from the center of a sphere reach by particles leaving */
        Vector_1D delta_1( bin_Centers.size() );
        RealDouble delta_2;

        /* Declare particle velocities */
        Vector_1D velp_1( bin_Centers.size() );
        RealDouble velp_2;

        /* Initialize particle diffusion coefficients */
        for ( unsigned int iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
            diffCoef_1[iBin] = physFunc::partDiffCoef( bin_Centers[iBin], temperature_K, pressure_Pa );
        }
        diffCoef_2 = physFunc::partDiffCoef( bin_R, temperature_K, pressure_Pa );

        /* Initialize mean distance */
        for ( unsigned int iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
            delta_1[iBin] = physFunc::delta_p( bin_Centers[iBin], physFunc::mass_sphere( bin_Centers[iBin], rho_1 ), temperature_K, pressure_Pa );
        }
        delta_2 = physFunc::delta_p( bin_R, physFunc::mass_sphere( bin_R, rho_2 ), temperature_K, pressure_Pa );

        /* Initialize particle velocities */
        for ( unsigned int iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
            velp_1[iBin] = physFunc::thermalSpeed( temperature_K, physFunc::mass_sphere( bin_Centers[iBin], rho_1 ) );
        }
        velp_2 = physFunc::thermalSpeed( temperature_K, physFunc::mass_sphere( bin_R, rho_2 ) );

        /* Initialize Brownian kernel */
        for ( unsigned int iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
            K_Brownian[iBin] = 4.0 * physConst::PI * ( bin_Centers[iBin] + bin_R ) * ( diffCoef_1[iBin] + diffCoef_2 ) / \
                               ( ( bin_Centers[iBin] + bin_R ) / ( bin_Centers[iBin] + bin_R + sqrt( delta_1[iBin] * delta_1[iBin] + delta_2 * delta_2 ) ) + \
                                 4.0 * ( diffCoef_1[iBin] + diffCoef_2 ) / ( sqrt( velp_1[iBin] * velp_1[iBin] + velp_2 + velp_2 ) * ( bin_Centers[iBin] + bin_R ) ) );
        }

        return K_Brownian;

    } /* End of buildBrownianKernel */

    Vector_2D buildBrownianKernel( RealDouble temperature_K, RealDouble pressure_Pa, Vector_1D bin_Centers_1, RealDouble rho_1 , Vector_1D bin_Centers_2, RealDouble rho_2 )
    {

        /* Returns the 2D coagulation kernel associated with Brownian coagulation */

        /* Declare coagulation kernel */
        Vector_2D K_Brownian;

        /* Declare particle diffusion coefficients */
        Vector_1D diffCoef_1( bin_Centers_1.size() );
        Vector_1D diffCoef_2( bin_Centers_2.size() );
       
        /* Declare mean distance from the center of a sphere reach by particles leaving */
        Vector_1D delta_1( bin_Centers_1.size() );
        Vector_1D delta_2( bin_Centers_2.size() );

        /* Declare particle velocities */
        Vector_1D velp_1( bin_Centers_1.size() );
        Vector_1D velp_2( bin_Centers_2.size() );

        /* Initialize particle diffusion coefficients */
        for ( unsigned int iBin = 0; iBin < bin_Centers_1.size(); iBin++ ) {
            diffCoef_1[iBin] = physFunc::partDiffCoef( bin_Centers_1[iBin], temperature_K, pressure_Pa );
        }
        for ( unsigned int iBin = 0; iBin < bin_Centers_2.size(); iBin++ ) {
            diffCoef_2[iBin] = physFunc::partDiffCoef( bin_Centers_2[iBin], temperature_K, pressure_Pa );
        }

        /* Initialize mean distance */
        for ( unsigned int iBin = 0; iBin < bin_Centers_1.size(); iBin++ ) {
            delta_1[iBin] = physFunc::delta_p( bin_Centers_1[iBin], physFunc::mass_sphere( bin_Centers_1[iBin], rho_1 ), temperature_K, pressure_Pa );
        }
        for ( unsigned int iBin = 0; iBin < bin_Centers_2.size(); iBin++ ) {
            delta_2[iBin] = physFunc::delta_p( bin_Centers_2[iBin], physFunc::mass_sphere( bin_Centers_2[iBin], rho_2 ), temperature_K, pressure_Pa );
        }

        /* Initialize particle velocities */
        for ( unsigned int iBin = 0; iBin < bin_Centers_1.size(); iBin++ ) {
            velp_1[iBin] = physFunc::thermalSpeed( temperature_K, physFunc::mass_sphere( bin_Centers_1[iBin], rho_1 ) );
        }
        for ( unsigned int iBin = 0; iBin < bin_Centers_2.size(); iBin++ ) {
            velp_2[iBin] = physFunc::thermalSpeed( temperature_K, physFunc::mass_sphere( bin_Centers_2[iBin], rho_2 ) );
        }

        /* Initialize Brownian kernel */
        for ( unsigned int iBin_1 = 0; iBin_1 < bin_Centers_1.size(); iBin_1++ ) {
            K_Brownian.push_back( Vector_1D( bin_Centers_2.size() ) );
            for ( unsigned int iBin_2 = 0; iBin_2 < bin_Centers_2.size(); iBin_2++ ) {
                K_Brownian[iBin_1][iBin_2] = 4.0 * physConst::PI * ( bin_Centers_1[iBin_1] + bin_Centers_2[iBin_2] ) * ( diffCoef_1[iBin_1] + diffCoef_2[iBin_2] ) / \
                                             ( ( bin_Centers_1[iBin_1] + bin_Centers_2[iBin_2] ) / ( bin_Centers_1[iBin_1] + bin_Centers_2[iBin_2] + sqrt( delta_1[iBin_1] * delta_1[iBin_1] + delta_2[iBin_2] * delta_2[iBin_2] ) ) \
                                               + 4.0 * ( diffCoef_1[iBin_1] + diffCoef_2[iBin_2] ) / ( sqrt( velp_1[iBin_1] * velp_1[iBin_1] + velp_2[iBin_2] + velp_2[iBin_2] ) * ( bin_Centers_1[iBin_1] + bin_Centers_2[iBin_2] ) ) );
            }
        }

        return K_Brownian;

    } /* End of buildBrownianKernel */

    Vector_1D buildTIKernel( RealDouble temperature_K, RealDouble pressure_Pa, Vector_1D bin_Centers, RealDouble rho_1, RealDouble bin_R, RealDouble rho_2 )
    {

        /* Returns the 1D coagulation kernel linked to turbulent initial motion.
         * Coagulation of particles occurs through this physical process because 
         * particles of different size accelerate differently */

        /* Declare coagulation kernel */
        Vector_1D K_TI( bin_Centers.size() );

        /* Initialize TI kernel */
        for ( unsigned int iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
            K_TI[iBin] = physConst::PI * pow ( physConst::EPSILON , 0.75 ) / ( physConst::g * pow( physFunc::kinVisc( temperature_K, pressure_Pa ), 0.25 ) ) * ( bin_Centers[iBin] + bin_R ) * ( bin_Centers[iBin] + bin_R ) * std::abs( physFunc::vFall( bin_Centers[iBin], rho_1, temperature_K, pressure_Pa ) - physFunc::vFall( bin_R, rho_2, temperature_K, pressure_Pa ) );
        }

        return K_TI;

    } /* End of buildTIKernel */
    
    Vector_2D buildTIKernel( RealDouble temperature_K, RealDouble pressure_Pa, Vector_1D bin_Centers_1, RealDouble rho_1, Vector_1D bin_Centers_2, RealDouble rho_2 )
    {

        /* Returns the 2D coagulation kernel linked to turbulent initial motion.
         * Coagulation of particles occurs through this physical process because 
         * particles of different size accelerate differently */

        /* Declare coagulation kernel */
        Vector_2D K_TI;

        /* Initialize TI kernel */
        for ( unsigned int iBin_1 = 0; iBin_1 < bin_Centers_1.size(); iBin_1++ ) {
            K_TI.push_back( Vector_1D( bin_Centers_2.size() ) );
            for ( unsigned int iBin_2 = 0; iBin_2 < bin_Centers_2.size(); iBin_2++ ) {
                K_TI[iBin_1][iBin_2] = physConst::PI * pow ( physConst::EPSILON , 0.75 ) / ( physConst::g * pow( physFunc::kinVisc( temperature_K, pressure_Pa ), 0.25 ) ) * ( bin_Centers_1[iBin_1] + bin_Centers_2[iBin_2] ) * ( bin_Centers_1[iBin_1] + bin_Centers_2[iBin_2] ) * std::abs( physFunc::vFall( bin_Centers_1[iBin_1], rho_1, temperature_K, pressure_Pa ) - physFunc::vFall( bin_Centers_2[iBin_2], rho_2, temperature_K, pressure_Pa ) );
            }
        }

        return K_TI;

    } /* End of buildTIKernel */

}

