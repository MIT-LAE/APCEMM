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

    void ISA( const RealDouble &z, RealDouble &pressure, \
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
    
    void ISA( const Vector_1D &z, Vector_1D &pressure, \
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

    void ISA( const RealDouble &z, RealDouble &pressure )
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
    
    void ISA( const Vector_1D &z, Vector_1D &pressure )
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
    
    void ISA_pAlt( RealDouble &z, const RealDouble &pressure )
    {

        /* DESCRIPTION: Find altitude corresponding to the prescribed 
         * pressure using dichotomy search algorithm */

        /* INPUTS:
         * RealDouble p : Pressure in Pa */

        const UInt nLUT = 500;
        const RealDouble zMAX = 50.0E+03;
        const RealDouble DZ = zMAX / RealDouble( nLUT );
        Vector_1D zLUT( nLUT, 0.0E+00 );
        Vector_1D pLUT( nLUT, 0.0E+00 );

        UInt i_z = 0;

        for ( i_z = 1; i_z < nLUT; i_z++ )
            zLUT[i_z] = zLUT[i_z-1] + DZ;

        ISA( zLUT, pLUT );

        i_z = 0;
        while ( i_z < nLUT ) {
            if ( pressure > pLUT[i_z] ) {
                /* Then: 
                 * pLUT[i_z-1] > pressure > pLUT[i_z]
                 * */
                break;
            }
            i_z += 1;
        }
        z = zLUT[i_z-1] + DZ * log( pressure / pLUT[i_z-1] ) / log( pLUT[i_z] / pLUT[i_z-1] );


    } /* End of ISA_pAlt */

    void ISA_pAlt( Vector_1D &z, const Vector_1D &pressure )
    {

        /* DESCRIPTION: Find altitude corresponding to the prescribed 
         * pressure using dichotomy search algorithm */

        /* INPUTS:
         * Vector_1D p : Pressure in Pa */

        const UInt nLUT = 500;
        const RealDouble zMAX = 50.0E+03;
        const RealDouble DZ = zMAX / RealDouble( nLUT );
        Vector_1D zLUT( nLUT, 0.0E+00 );
        Vector_1D pLUT( nLUT, 0.0E+00 );

        UInt i_z = 0;

        for ( i_z = 1; i_z < nLUT; i_z++ )
            zLUT[i_z] = zLUT[i_z-1] + DZ;

        ISA( zLUT, pLUT );

        for ( UInt iN = 0; iN < z.size(); iN++ ) {
            i_z = 0;
            while ( i_z < nLUT ) {
                if ( pressure[iN] > pLUT[i_z] ) {
                    /* Then: 
                     * pLUT[i_z-1] > pressure[iN] > pLUT[i_z]
                     * */
                    break;
                }
                i_z += 1;
            }
            z[iN] = zLUT[i_z-1] + DZ * log( pressure[iN] / pLUT[i_z-1] ) / log( pLUT[i_z] / pLUT[i_z-1] );
        }


    } /* End of ISA_pAlt */

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

    UInt nearestNeighbor( float xq[], const float &x , size_t xq_size) {

        /* DESCRIPTION: Finds the closest of x in xq, returning the index */

        /* INPUTS:
         * (1-D Vector) xq: query values
         * float x:  desired value */

        /*This function passes in a float array which means the size of the array is unknown. Therefore, we actually 
        have no idea when to stop if the query value is past the end of the array.
        */
        UInt i_Z = 0;

        /* Identify increasing direction */
        if ( xq[0]-xq[1]<0 ) {
            while ( xq[i_Z] <= x ) {

                if(i_Z == xq_size-1) return i_Z;
                i_Z += 1;
            }
            /* Check if out of bounds*/
            if(x < xq[0]){
                return 0;
            }
            /* Check if previous altitude closer */
            else if ( xq[i_Z]-x >= x-xq[i_Z-1] ) {
                i_Z -= 1;
            }
        }
        else {
            while ( xq[i_Z] >= x ) {
                if(i_Z == xq_size-1) return i_Z;
                i_Z += 1;
            }
            /* Check if out of bounds */
            if(x>xq[0]){
                return 0;
            }
            else if ( x-xq[i_Z] >= xq[i_Z-1]-x ) {
                i_Z -= 1;
            }
        }
        return i_Z;

    } /* End of nearestNeighbor */

    UInt nearestNeighbor( Vector_1D xq, const float &x ) {
        //Just a temporary solution to get rid of this redundant high maint code. Ideally want to template it later or something - Michael
        float* floatArray = new float[xq.size()];
        for (int i = 0 ; i < xq.size(); i++)
        {
            floatArray[i] = (float) xq[i];
        }
        UInt i = nearestNeighbor(floatArray, x, xq.size());
        delete[] floatArray;
        return i;
    } /* End of nearestNeighbor */ 

    float linearInterp( float xq[], float yq[], const float &x , size_t xq_size) {

        /* DESCRIPTION: Linearly interpolated around the desired x value */

        /* INPUTS:
         * float xq: x query values
         * float yq: y query values
         * float x: desired value to interpolate around*/

        /* Initialize variables */
        UInt i_X;
        UInt i_X2;
        float y;

        /* Find closest point */
        i_X = nearestNeighbor( xq, x , xq_size);
        i_X2 = i_X;

        /* Disallow extrapolation for now*/
        bool out_of_bounds = xq[0]-xq[1]<0 ? x < xq[0] || x > xq[xq_size-1] :
                                             x > xq[0] || x < xq[xq_size-1];
        if(out_of_bounds) {
            throw std::range_error("met::linearInterp : x out of range of xq. Extrapolation not supported.");
        }
        /* Check direction xq increasing */
        while ( xq[i_X2]==xq[i_X] ) {
        if ( xq[0]-xq[1]<0 ) {
            /* Find the next closest point */
            if ( xq[i_X] > x ) {
                i_X2 = i_X2 - 1;
            }
            else {
                i_X2 = i_X2 + 1;
            }
        }
        else {
            /* Find the next closest point */
            if ( xq[i_X] >= x ) {
                i_X2 = i_X2 + 1;
            }
            else {
                i_X2 = i_X2 - 1;
            }
        }
        }

        /* Linear interpolation */
        y = yq[i_X] + ( x-xq[i_X] ) * ( yq[i_X2]-yq[i_X] ) / ( xq[i_X2]-xq[i_X] );
        return y;

    }

    RealDouble satdepth_calc( float RHw[], float T[], float alt[], UInt iFlight, UInt var_length ) {

        /* DESCRIPTION: Finds the saturation depth RHw and T profiles */

        /* INPUTS:
         * float RHw[]: Relative humidity wrt water [%] profile
         * float T[]: Temperature [K] profile
         * float alt[]: Altitude [m] associated with profile
         * UInt iFlight: index at flight altitude
         * UInt var_length: length of RHw and T profiles */

        RealDouble satdepth = 0.00E+00;
        RealDouble curdepth = 0.00E+00;
        int iCur = iFlight;
        RealDouble RHi_cur, RHi_prev;

        /* Loop over altitudes till sat depth found or end of profile reached */
        while ( satdepth==0.00E+00 && iCur >= 0 ) {

            /* Calculate RHi at current altitude */
            RHi_cur = RHw[iCur] * physFunc::pSat_H2Ol( T[iCur] ) / physFunc::pSat_H2Os( T[iCur] );
            /* Check if current altitude is subsaturated */
            if ( iCur != iFlight && RHi_cur < 100 ) {

                /* Get distance from previous altitude for saturation depth */
                //Because met input automatically sets RH to nearest point,
                //We need to "interpolate" on the border point between the saturated and non-sat layer.

                //satdepth = alt[iFlight] - alt[iCur];

                RHi_prev = RHw[iCur + 1] * physFunc::pSat_H2Ol( T[iCur + 1] ) / physFunc::pSat_H2Os( T[iCur + 1] );
                satdepth = alt[iFlight] - (alt[iCur] + (100 - RHi_cur)/(RHi_prev - RHi_cur) * (alt[iCur + 1] - alt[iCur]));
                break;
            }

            /* Check first point is ISS */
            if ( RHi_cur < 100 && iCur==iFlight ) {
                satdepth = 1.0; /* Set to some arbitrary value, contrail should not survive VS anyway */
                break;
            }

            /* Iterate and check satdepth */
            iCur = iCur - 1;
            
            // std::cout << alt[iCur] << "m, " << curdepth << "m" << std::endl;
        }
        /* Check a genuine satdepth found */
        if ( iCur < 0) {
            throw std::out_of_range("In met::satdepth_calc: No end of ice supersaturated layer found.");
        }
        else if (satdepth > YLIM_DOWN){
            throw std::out_of_range("In met::satdepth_calc: Ice supersatured layer depth exceeds domain limits.");
        }

        return satdepth;

    } /* End of satdepth_calc */ 

}

/* End of MetFunction.cpp */
