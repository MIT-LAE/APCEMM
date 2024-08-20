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

#include "Util/PhysConstant.hpp"
#include "Util/PhysFunction.hpp"
#include "Util/MetFunction.hpp"

namespace met
{

    void ISA( const double &z, double &pressure, \
              double &temperature )
    {
        
        /* DESCRIPTION: Implements the mathematical representation of the
         * lapse rate atmospheric equations for ambient temperature and 
         * pressure for the input geopotential altitude. */

        /* INPUTS:
         * double z : Geopotential altitude in m */

        double expon = 0.0E+00;
        double theta = 0.0E+00;

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

    void ISA( const double &z, double &pressure )
    {
        
        /* DESCRIPTION: Implements the mathematical representation of the
         * lapse rate atmospheric equations for ambient temperature and 
         * pressure for the input geopotential altitude. */

        /* INPUTS:
         * double z : Geopotential altitude in m */

        double expon = 0.0E+00;
        double theta = 0.0E+00;
        double temperature = 0.0E+00;

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
    
    void ISA_pAlt( double &z, const double &pressure )
    {

        /* DESCRIPTION: Find altitude corresponding to the prescribed 
         * pressure using dichotomy search algorithm */

        /* INPUTS:
         * double p : Pressure in Pa */

        const UInt nLUT = 500;
        const double zMAX = 50.0E+03;
        const double DZ = zMAX / double( nLUT );
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
        const double zMAX = 50.0E+03;
        const double DZ = zMAX / double( nLUT );
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

    double ComputeLapseRate( const double TEMP, const double RHi, \
                                 const double DEPTH )
    {

        /* DESCRIPTION: Computes the temperature lapse rate from the 
         * temperature and RH at flight level and the depth of the 
         * supersaturated layer, assuming a constant background H2O 
         * concentration */

        /* INPUTS:
         * double TEMP  : Temperature at flight level in [K]
         * double RHi   : RH w.r.t. ice at flight level in [%]
         * double DEPTH : Depth of the moist layer in [m] */

        /* Solve pSat_H2Os(T)/T = x_star using a Newton-Raphson iteration 
         * scheme */
        UInt counter    = 0;
        double T    = TEMP;
        double Tpre = 0.0E+00;
        double pSat = 0.0E+00;
        double pSat_= 0.0E+00;

        /* H2O concentration, assumed constant, computed from flight level
         * conditions in molec/cm^3 */
        const double H2O = RHi / double(100.0) * \
            physFunc::pSat_H2Os( TEMP ) / physConst::kB / TEMP * 1.00E-06;

        const double x_star = H2O * physConst::kB * 1.00E+06;

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

    std::size_t nearestNeighbor( const Vector_1D& xq, double x ) {

        /* DESCRIPTION: Finds the closest of x in xq, returning the index */

        /* INPUTS:
         * (1-D Vector) xq: query values, must be sorted
         * double x:  desired value */

        double diff = std::numeric_limits<double>::max();
        for(std::size_t i = 0; i < xq.size(); i++) {
            double newDiff = std::abs(xq[i] - x);
            if(newDiff >= diff) {
                return i - 1;
            }
            diff = newDiff;
        }
        return xq.size() - 1;
    } /* End of nearestNeighbor */

    double linearInterp( const Vector_1D& xq, const Vector_1D& yq, double x ) {

        /* DESCRIPTION: Linearly interpolated around the desired x value */

        /* INPUTS:
         * float xq: x query values
         * float yq: y query values
         * float x: desired value to interpolate around*/

        /* Initialize variables */
        int i_X, i_X2;
        auto xq_size = xq.size();

        /* Find closest point */
        i_X = nearestNeighbor( xq, x );
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
        return yq[i_X] + ( x-xq[i_X] ) * ( yq[i_X2]-yq[i_X] ) / ( xq[i_X2]-xq[i_X] );
    }

    double linearInterp(double x1, double y1, double x2, double y2, double xq, bool support_extrapolation){
        //Check for out of range inputs
        if(!support_extrapolation){
            if( (x1 < x2) && (xq < x1 || xq > x2)){
                throw std::range_error("In met::linearInterp: Query point is out of range! Extrapolation support set to false.");
            }
            else if( (x1 > x2) && (xq > x1 || xq < x2)){
                throw std::range_error("In met::linearInterp: Query point is out of range! Extrapolation support set to false.");
            }
        }

        return y1 + (xq - x1) / (x2 - x1) * (y2 - y1);
    }
    double linInterpMetData(const Vector_1D& altitude_init, const Vector_1D& metVar_init, double altitude_query){
        // Alt input is increasing
        std::size_t i_Z = nearestNeighbor( altitude_init, altitude_query );
        std::size_t idx_x1;

        //Edge cases
        const double epsilon = 1e-3;
        if((i_Z == 0 || i_Z == altitude_init.size() - 1) && std::abs(altitude_init[i_Z] - altitude_query) < epsilon) {
            return metVar_init[i_Z];
        }

        if(altitude_init[0] < altitude_init[1]){
            idx_x1 = altitude_query > altitude_init[i_Z] ? i_Z : i_Z - 1;
        }
        // Alt input is decreasing
        else {
            idx_x1 = altitude_query <= altitude_init[i_Z] ? i_Z : i_Z - 1;
        }
        

        std::size_t idx_x2 = idx_x1 + 1;
        if(idx_x1 < 0) { 
            throw std::range_error("Input flight altitude out of range of met. data!"); 
        }
        if(idx_x2 >= altitude_init.size()) { 
            throw std::range_error("Input flight altitude out of range of met. data!"); 
        }
        /* Loop round horizontal coordinates to assign temperature */
        double metVar_local = met::linearInterp(altitude_init[idx_x1], metVar_init[idx_x1], \
                                                altitude_init[idx_x2], metVar_init[idx_x2], \
                                                altitude_query);
        return metVar_local;
    }
    
    double satdepth_calc( const Vector_1D& RHw, const Vector_1D& T, const Vector_1D& alt, int iFlight, double YLIM_DOWN ) {

        /* DESCRIPTION: Finds the saturation depth RHw and T profiles */

        /* INPUTS:
         * RHw: Relative humidity wrt water [%] profile
         * T: Temperature [K] profile
         * alt: Altitude [m] associated with profile
         * iFlight: index at flight altitude
         * YLIM_DOWN: bottom limit of domain */

        double satdepth = 0.00E+00;
        // double curdepth = 0.00E+00;
        int iCur = iFlight;
        double RHi_cur, RHi_prev;

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
            throw std::out_of_range(std::string("\nIn met::satdepth_calc: Ice supersatured layer depth exceeds domain limits. \n") + 
                                    "YLIM_DOWN: " + std::to_string(YLIM_DOWN) + "\n" + 
                                    "satdepth: " + std::to_string(satdepth));
        }

        return satdepth;

    } /* End of satdepth_calc */ 


    newXCoordsPair calcNewXCoords(const Vector_1D& dy_old, const Vector_1D& dy_new, const Vector_1D& x0_old, const Vector_1D& dx_old, int nx){
        int ny = dy_old.size();
        Vector_1D dx_new(ny);
        Vector_1D x0_new(ny);

        for(int j = 0; j < ny; j++) {
            dx_new[j] = dx_old[j] * dy_new[j] / dy_old[j];
            x0_new[j] = x0_old[j] - (dx_new[j] - dx_old[j]) * static_cast<double>(nx) / 2.0;
        }

        return newXCoordsPair {std::move(x0_new), std::move(dx_new)};
    }
}

/* End of MetFunction.cpp */
