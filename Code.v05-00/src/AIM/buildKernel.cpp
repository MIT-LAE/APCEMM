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

#include "AIM/buildKernel.hpp"
#include "Util/PhysConstant.hpp"
#include "Util/PhysFunction.hpp"

namespace AIM
{

    Vector_1D buildBrownianKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers, double rho_1 , double bin_R, double rho_2 )
    {

        /* Returns the 1D coagulation kernel associated with Brownian coagulation */

        /* Declare coagulation kernel */
        Vector_1D K_Brownian( bin_Centers.size() );

        /* Declare particle diffusion coefficients */
        Vector_1D  diffCoef_1( bin_Centers.size() );
        double diffCoef_2;
       
        /* Declare mean distance from the center of a sphere reach by particles leaving */
        Vector_1D delta_1( bin_Centers.size() );
        double delta_2;

        /* Declare particle velocities */
        Vector_1D velp_1( bin_Centers.size() );
        double velp_2;

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

    Vector_2D buildBrownianKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers_1, double rho_1 , Vector_1D const &bin_Centers_2, double rho_2 )
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
    
    Vector_1D buildDEKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers, double rho_1, double bin_R, double rho_2, Vector_1D const &K_Brow )
    {

        /* Returns the 1D coagulation kernel linked to convective 
         * Brownian diffusion enhancement.
         * When a large particle falls through the air, eddies 
         * created in its wake enhance diffusion of other particles
         * to its surface */

        /* Declare DE kernel */
        Vector_1D K_DE( bin_Centers.size() );

        /* Initialize DE kernel */
        for ( unsigned int iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
            if ( bin_Centers[iBin] <= bin_R ) {
                if ( physFunc::Reynolds_p( bin_R, rho_2, temperature_K, pressure_Pa ) <= 1 ) {
                    K_DE[iBin] = K_Brow[iBin] * 0.45 * pow( physFunc::Reynolds_p( bin_R, rho_2, temperature_K, pressure_Pa ), 1.0 / double(3.0) ) * pow( physFunc::Schmidt_p( bin_Centers[iBin], temperature_K, pressure_Pa ), 1.0 / double(3.0) );
                } else {
                    K_DE[iBin] = K_Brow[iBin] * 0.45 * pow( physFunc::Reynolds_p( bin_R, rho_2, temperature_K, pressure_Pa ), 0.5               ) * pow( physFunc::Schmidt_p( bin_Centers[iBin], temperature_K, pressure_Pa ), 1.0 / double(3.0) );
                }
            } else {
                if ( physFunc::Reynolds_p( bin_Centers[iBin], rho_1, temperature_K, pressure_Pa ) <= 1 ) {
                    K_DE[iBin] = K_Brow[iBin] * 0.45 * pow( physFunc::Reynolds_p( bin_Centers[iBin], rho_1, temperature_K, pressure_Pa ), 1.0 / double(3.0) ) * pow( physFunc::Schmidt_p( bin_R, temperature_K, pressure_Pa ), 1.0 / double(3.0) );
                } else {
                    K_DE[iBin] = K_Brow[iBin] * 0.45 * pow( physFunc::Reynolds_p( bin_Centers[iBin], rho_1, temperature_K, pressure_Pa ), 0.5               ) * pow( physFunc::Schmidt_p( bin_R, temperature_K, pressure_Pa ), 1.0 / double(3.0) );
                }
            }
        }

        return K_DE;

    } /* End of buildDEKernel */

    Vector_2D buildDEKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers_1, double rho_1, Vector_1D const &bin_Centers_2, double rho_2, Vector_2D const &K_Brow )
    {

        /* Returns the 2D coagulation kernel linked to convective 
         * Brownian diffusion enhancement.
         * When a large particle falls through the air, eddies 
         * created in its wake enhance diffusion of other particles
         * to its surface */

        /* Declare DE kernel */
        Vector_2D K_DE;

        /* Initialize DE kernel */
        for ( unsigned int iBin_1 = 0; iBin_1 < bin_Centers_1.size(); iBin_1++ ) {
            K_DE.push_back( Vector_1D( bin_Centers_2.size() ) );
            for ( unsigned int iBin_2 = 0; iBin_2 < bin_Centers_2.size(); iBin_2++ ) {
                if ( bin_Centers_1[iBin_1] <= bin_Centers_2[iBin_2]) {
                    if ( physFunc::Reynolds_p( bin_Centers_2[iBin_2], rho_2, temperature_K, pressure_Pa ) <= 1 ) {
                        K_DE[iBin_1][iBin_2] = K_Brow[iBin_1][iBin_2] * 0.45 * pow( physFunc::Reynolds_p( bin_Centers_2[iBin_2], rho_2, temperature_K, pressure_Pa ), 1.0 / double(3.0) ) * pow( physFunc::Schmidt_p( bin_Centers_1[iBin_1], temperature_K, pressure_Pa ), 1.0 / double(3.0) );
                    } else {
                        K_DE[iBin_1][iBin_2] = K_Brow[iBin_1][iBin_2] * 0.45 * pow( physFunc::Reynolds_p( bin_Centers_2[iBin_2], rho_2, temperature_K, pressure_Pa ), 0.5               ) * pow( physFunc::Schmidt_p( bin_Centers_1[iBin_1], temperature_K, pressure_Pa ), 1.0 / double(3.0) );
                    }
                } else {
                    if ( physFunc::Reynolds_p( bin_Centers_1[iBin_1], rho_1, temperature_K, pressure_Pa ) <= 1 ) {
                        K_DE[iBin_1][iBin_2] = K_Brow[iBin_1][iBin_2] * 0.45 * pow( physFunc::Reynolds_p( bin_Centers_1[iBin_1], rho_1, temperature_K, pressure_Pa ), 1.0 / double(3.0) ) * pow( physFunc::Schmidt_p( bin_Centers_2[iBin_2], temperature_K, pressure_Pa ), 1.0 / double(3.0) );
                    } else {
                        K_DE[iBin_1][iBin_2] = K_Brow[iBin_1][iBin_2] * 0.45 * pow( physFunc::Reynolds_p( bin_Centers_1[iBin_1], rho_1, temperature_K, pressure_Pa ), 0.5               ) * pow( physFunc::Schmidt_p( bin_Centers_2[iBin_2], temperature_K, pressure_Pa ), 1.0 / double(3.0) );
                    }

                }
            }
        }

        return K_DE;

    } /* End of buildDEKernel */
    
    Vector_1D buildGCKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers, double rho_1, double bin_R, double rho_2 )
    {

        /* Returns the 1D coagulation kernel linked to gravitational 
         * collection. 
         * When two particles of different size fall, one may catch 
         * up with and collide with the other */ 

        /* Declare GC kernel */
        Vector_1D K_GC( bin_Centers.size() );

        /* Initialize GC kernel */
        for ( unsigned int iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
            K_GC[iBin] = physFunc::E_agg( bin_Centers[iBin], rho_1, bin_R, rho_2, temperature_K, pressure_Pa ) * physConst::PI * ( bin_Centers[iBin] + bin_R ) * ( bin_Centers[iBin] + bin_R ) * std::abs( physFunc::vFall( bin_Centers[iBin], rho_1, temperature_K, pressure_Pa ) - physFunc::vFall( bin_R, rho_2, temperature_K, pressure_Pa ) );
        }

        return K_GC;

    } /* End of buildGCKernel */

    Vector_2D buildGCKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers_1, double rho_1, Vector_1D const &bin_Centers_2, double rho_2 )
    {
        
        /* Returns the 2D coagulation kernel linked to gravitational 
         * collection. 
         * When two particles of different size fall, one may catch 
         * up with and collide with the other */ 

        /* Declare GC kernel */
        Vector_2D K_GC;

        /* Initialize GC kernel */
        for ( unsigned int iBin_1 = 0; iBin_1 < bin_Centers_1.size(); iBin_1++ ) {
            K_GC.push_back( Vector_1D( bin_Centers_2.size() ) );
            for ( unsigned int iBin_2 = 0; iBin_2 < bin_Centers_2.size(); iBin_2++ ) {
                K_GC[iBin_1][iBin_2] = physFunc::E_agg( bin_Centers_1[iBin_1], rho_1, bin_Centers_2[iBin_2], rho_2, temperature_K, pressure_Pa ) * physConst::PI * ( bin_Centers_1[iBin_1] + bin_Centers_2[iBin_2] ) * ( bin_Centers_1[iBin_1] + bin_Centers_2[iBin_2]) * std::abs( physFunc::vFall( bin_Centers_1[iBin_1], rho_1, temperature_K, pressure_Pa ) - physFunc::vFall( bin_Centers_2[iBin_2], rho_2, temperature_K, pressure_Pa ) );
            }
        }

        return K_GC;

    } /* End of buildGCKernel */


    Vector_1D buildTIKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers, double rho_1, double bin_R, double rho_2 )
    {

        /* Returns the 1D coagulation kernel linked to turbulent inertial motion.
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
    
    Vector_2D buildTIKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers_1, double rho_1, Vector_1D const &bin_Centers_2, double rho_2 )
    {

        /* Returns the 2D coagulation kernel linked to turbulent inertial motion.
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

    Vector_1D buildTSKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers, double rho_1, double bin_R, double rho_2 )
    {

        /* Returns the 1D coagulation kernel linked to turbulent shear. 
         * Wind shear in turbulent air causes particles moving with the 
         * air to collide */

        /* Declare TS kernel */
        Vector_1D K_TS( bin_Centers.size() );

        /* Initialize TS kernel */
        for ( unsigned int iBin = 0; iBin < bin_Centers.size(); iBin++ ) {
            K_TS[iBin] = sqrt( 8.0 * physConst::PI * physConst::EPSILON / ( 15.0 * physFunc::kinVisc( temperature_K, pressure_Pa ) ) ) * ( bin_Centers[iBin] + bin_R ) * ( bin_Centers[iBin] + bin_R ) * ( bin_Centers[iBin] + bin_R );
        }

        return K_TS;

    } /* End of buildTSKernel */
    
    Vector_2D buildTSKernel( double temperature_K, double pressure_Pa, Vector_1D const &bin_Centers_1, double rho_1, Vector_1D const &bin_Centers_2, double rho_2 )
    {

        /* Returns the 2D coagulation kernel linked to turbulent shear. 
         * Wind shear in turbulent air causes particles moving with the 
         * air to collide */

        /* Declare TS kernel */
        Vector_2D K_TS;

        /* Initialize TS kernel */
        for ( unsigned int iBin_1 = 0; iBin_1 < bin_Centers_1.size(); iBin_1++ ) {
            K_TS.push_back( Vector_1D( bin_Centers_2.size() ) );
            for ( unsigned int iBin_2 = 0; iBin_2 < bin_Centers_2.size(); iBin_2++ ) {
                K_TS[iBin_1][iBin_2] = sqrt( 8.0 * physConst::PI * physConst::EPSILON / ( 15.0 * physFunc::kinVisc( temperature_K, pressure_Pa ) ) ) * ( bin_Centers_1[iBin_1] + bin_Centers_2[iBin_2] ) * ( bin_Centers_1[iBin_1] + bin_Centers_2[iBin_2] ) * ( bin_Centers_1[iBin_1] + bin_Centers_2[iBin_2] );
            }
        }

        return K_TS;

    } /* End of buildTSKernel */
}

