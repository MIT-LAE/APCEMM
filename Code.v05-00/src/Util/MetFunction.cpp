/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* MetFunction Program File                                         */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 11/2/2018                                 */
/* File                 : MetFunction.cpp                           */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Util/MetFunction.hpp"
 

namespace met
{

    void ISA( const RealDouble z, RealDouble &pressure, \
              RealDouble &temperature )
    {
        
        /* DESCRIPTION: Implements the mathematical representation of the
         * lapse rate atmospheric equations for ambient temperature and 
         * pressure for the input geopotential altitude. */

        /* INPUTS:
         * RealDouble z : Geopotential altitude in m */

        RealDouble expon = 0.0E+00;
        RealDouble theta = 0.0E+00;

        if ( z > ISA_HTS ) {
            temperature = ISA_T0 - ISA_LAPSERATE * ISA_HTS;
            expon = exp( physConst::g / ( physConst::R_Air * temperature ) * ( ISA_HTS - z ) );
        }
        else {
            temperature = ISA_T0 - ISA_LAPSERATE * z      ;
            expon = 1.0E+00;
        }

        theta = temperature / ISA_T0;

        pressure = ISA_P0 \
                   * pow( theta, physConst::g / ( ISA_LAPSERATE * physConst::R_Air ) ) \
                   * expon;

    } /* End of ISA */
    
    void ISA( const Vector_1D z, Vector_1D &pressure, \
              Vector_1D &temperature )
    {

        /* DESCRIPTION: Implements the mathematical representation of the
         * lapse rate atmospheric equations for ambient temperature and 
         * pressure for the input geopotential altitude. */

        /* INPUTS:
         * Vector_1D z : Geopotential altitude in m */

        Vector_1D expon( z.size(), 0.0E+00 );
        Vector_1D theta( z.size(), 0.0E+00 );

        for ( unsigned int i_z = 0; i_z < z.size(); i_z++ ) {
            if ( z[i_z] > ISA_HTS ) {
                temperature[i_z] = ISA_T0 - ISA_LAPSERATE * ISA_HTS;
                expon[i_z] = exp( physConst::g / ( physConst::R_Air * temperature[i_z] ) \
                             * ( ISA_HTS - z[i_z] ) );
            }
            else {
                temperature[i_z] = ISA_T0 - ISA_LAPSERATE * z[i_z] ;
                expon[i_z] = 1.0E+00;
            }

            theta[i_z] = temperature[i_z] / ISA_T0;

            pressure[i_z] = ISA_P0 \
                            * pow( theta[i_z], physConst::g / ( ISA_LAPSERATE * physConst::R_Air ) ) \
                            * expon[i_z];
        }

    } /* End of ISA */

    void ISA( const RealDouble z, RealDouble &pressure )
    {
        
        /* DESCRIPTION: Implements the mathematical representation of the
         * lapse rate atmospheric equations for ambient temperature and 
         * pressure for the input geopotential altitude. */

        /* INPUTS:
         * RealDouble z : Geopotential altitude in m */

        RealDouble expon = 0.0E+00;
        RealDouble theta = 0.0E+00;
        RealDouble temperature = 0.0E+00;

        if ( z > ISA_HTS ) {
            temperature = ISA_T0 - ISA_LAPSERATE * ISA_HTS;
            expon = exp( physConst::g / ( physConst::R_Air * temperature ) * ( ISA_HTS - z ));
        }
        else {
            temperature = ISA_T0 - ISA_LAPSERATE * z      ;
            expon = 1.0E+00;
        }

        theta = temperature / ISA_T0;

        pressure = ISA_P0 \
                   * pow( theta, physConst::g / ( ISA_LAPSERATE * physConst::R_Air ) ) \
                   * expon;

    } /* End of ISA */
    
    void ISA( const Vector_1D z, Vector_1D &pressure )
    {

        /* DESCRIPTION: Implements the mathematical representation of the
         * lapse rate atmospheric equations for ambient temperature and 
         * pressure for the input geopotential altitude. */

        /* INPUTS:
         * Vector_1D z : Geopotential altitude in m */

        Vector_1D expon( z.size(), 0.0E+00 );
        Vector_1D theta( z.size(), 0.0E+00 );
        Vector_1D temperature( z.size(), 0.0E+00 );

        for ( unsigned int i_z = 0; i_z < z.size(); i_z++ ) {
            if ( z[i_z] > ISA_HTS ) {
                temperature[i_z] = ISA_T0 - ISA_LAPSERATE * ISA_HTS;
                expon[i_z] = exp( physConst::g / ( physConst::R_Air * temperature[i_z] ) \
                             * ( ISA_HTS - z[i_z] ) );
            }
            else {
                temperature[i_z] = ISA_T0 - ISA_LAPSERATE * z[i_z] ;
                expon[i_z] = 1.0E+00;
            }

            theta[i_z] = temperature[i_z] / ISA_T0;

            pressure[i_z] = ISA_P0 \
                            * pow( theta[i_z], physConst::g / ( ISA_LAPSERATE * physConst::R_Air ) ) \
                            * expon[i_z];
        }

    } /* End of ISA */

    RealDouble ComputeLapseRate( const RealDouble TEMP, const RealDouble RHi, \
                                 const RealDouble DEPTH )
    {

        /* DESCRIPTION: Computes the temperature lapse rate from the 
         * temperature and RH at flight level and the depth of the 
         * supersaturated layer, assuming a constant background H2O 
         * concentration */

        /* INPUTS:
         * RealDouble TEMP  : Temperature at flight level in [K]
         * RealDouble RHi   : RH w.r.t. ice at flight level in [%]
         * RealDouble DEPTH : Depth of the moist layer in [m] */

        /* Solve pSat_H2Os(T)/T = x_star using a Newton-Raphson iteration 
         * scheme */
        UInt counter    = 0;
        RealDouble T    = TEMP;
        RealDouble Tpre = 0.0E+00;
        RealDouble pSat = 0.0E+00;
        RealDouble pSat_= 0.0E+00;

        /* H2O concentration, assumed constant, computed from flight level
         * conditions in molec/cm^3 */
        const RealDouble H2O = RHi / RealDouble(100.0) * \
            physFunc::pSat_H2Os( TEMP ) / physConst::kB / TEMP * 1.00E-06;

        const RealDouble x_star = H2O * physConst::kB * 1.00E+06;

        while ( counter < 100 ) {

            Tpre = T;
            counter += 1;

            pSat = physFunc::pSat_H2Os( T ); 
            pSat_= physFunc::dpSat_H2Os( T );
            T = T - ( pSat / T - x_star ) / ( T * pSat_ - pSat ) * T * T;

            if ( ABS( T - Tpre ) < 1.00E-03 )
                break;

        }

        return ( TEMP - T ) / DEPTH;

    } /* End of ComputeLapseRate */

}

/* End of MetFunction.cpp */
