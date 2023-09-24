/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*                        AIrcraft Microphysics                     */
/*                              (AIM)                               */
/*                                                                  */
/* Nucleation Program File                                          */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 10/1/2018                                 */
/* File                 : Nucleation.cpp                            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "AIM/Nucleation.hpp"

namespace AIM
{

    /* Set of functions that calculates the binary nucleation rate and radius 
     * of the critical nucleation cluster using the parameterization of...
     *
     * Vehkamaki, H., M. Kulmala, I. Napari, K. E. J. Lehtinen, C. Timmreck, 
     * M. Noppel, and A. Laaksonen. "An Improved Parameterization for Sulfuric 
     * Acid-Water Nucleation Rates for Tropospheric and Stratospheric Conditions." 
     * Journal of Geophysical Research-Atmospheres 107, no. D22 (2002). */

    /* Fit is valid in the temperature range 230.15 - 320.15 K but can be extended
     * down to 190K (Vehkamaki, 2002) */

    /* If f(x) = a_0 + a_1 * x + a_2 * x^2 + a_3 * x^3 + a_4 * x_4 + ...
     * Data is computed as:
     * f(x) = a_0 + x * ( a_1 + x * ( a_2 + x * ( a_3 + x * ( a_4 + ... ) */

    double sigma( double x_m, double T )
    {

        /* DESCRIPTION:
         * Returns the surface tension \sigma in J/m^2 for all sulfuric acid mass
         * fractions x_m where the solution is liquid */

        /* INPUTS:
         * - double x_m :: sulfuric acid mass fraction
         * - double T   :: temperature in K
         *
         * OUTPUTS:
         * - double :: surface tension in J/m^2 */
        
        /* Apply limitations */
        if ( T > 305.15 )
            T = 305.15;

        if ( T < 190.0 )
            T = 190.0;

        const double a = + 1.1864E-01 + x_m * ( - 1.1651E-01 + x_m * ( + 7.6852E-01 \
                             + x_m * ( - 2.40909E-00 + x_m * ( + 2.95434E-00 + x_m * ( - 1.25852E-00 ) ) ) ) );
        const double b = - 1.5709E-04 + x_m * ( + 4.0102E-04 + x_m * ( - 2.3995E-03 \
                             + x_m * ( + 7.611235E-03 + x_m * ( - 9.37386E-03 + x_m * ( + 3.89722E-03 ) ) ) ) );

        return a + T * b;

    } /* End of sigma */

    double rho( double x_m, double T )
    {

        /* DESCRIPTION:
         * Returns the density for the sulfuric acid solution solution in kg/m^3 */

        /* INPUTS:
         * - double x_m :: sulfuric acid mass fraction
         * - double T   :: temperature in K
         *
         * OUTPUTS:
         * - double :: density in kg/m^3 */

        /* Note: Coefficients have been scaled by 1.00E+03 compared to the original paper
         * to account for the conversion from g/cm^3 to kg/m^3 */

        /* Apply limitations */
        if ( T > 305.15 )
            T = 305.15;

        if ( T < 190.0 )
            T = 190.0;

        const double a = + 7.681724E+02 + x_m * ( + 2.1847140E+03 + x_m * ( + 7.163002E+03 \
                             + x_m * ( - 4.431447E+04 + x_m * ( + 8.875606E+04 + x_m * ( - 7.573729E+04 \
                             + x_m * ( + 2.343228E+04 ) ) ) ) ) );
        const double b = + 1.808225E+00 + x_m * ( - 9.294656E+00 + x_m * ( - 3.742148E+01 \
                             + x_m * ( + 2.565321E+02 + x_m * ( - 5.362872E+02 + x_m * ( + 4.857736E+02 \
                             + x_m * ( - 1.629592E+02 ) ) ) ) ) );
        const double c = - 3.478524E-03 + x_m * ( + 1.335867E-02 + x_m * ( + 5.195706E-02 \
                             + x_m * ( - 3.717636E-01 + x_m * ( + 7.990811E-01 + x_m * ( - 7.458060E-01 \
                             + x_m * ( + 2.58139E-01 ) ) ) ) ) );

        return a + T * ( b + T * c );

    } /* End of rho */

    double x_star( double T, double RH, double nSulf )
    {

        /* DESCRIPTION:
         * Returns the mole fraction of sulfuric acid in the critical cluster */

        /* INPUTS:
         * - double T     :: temperature in K
         * - double RH    :: relative humidity on a ( 0, 1 ) scale
         * - double nSulf :: total gas phase concentration of sulfuric acid in molecules/cm^3
         *
         * OUTPUTS:
         * - double :: mole fraction of sulfuric acid */

        /* Apply limitations */
        if ( nSulf >= 1.0E+11 )
            nSulf = 1.0E+11;

        if ( RH > 1.0E+00 )
            RH = 1.0E+00;

        if ( RH < 1.0E-04 )
            RH = 1.0E-04;

        if ( T > 305.15 )
            T = 305.15;

        if ( T < 190.0 )
            T = 190.0;

        const double logRH    = log(RH);
        const double lognSulf = log(nSulf);

        const double a = + 7.40997E-01 - 3.49998E-03 * lognSulf + logRH * ( + 2.01048E-03 \
                             + logRH * ( + 1.57407E-03 + logRH * ( + 1.84403E-04 ) ) );
        const double b = - 2.66379E-03 + 5.04022E-05 * lognSulf + logRH * ( - 1.83289E-04 \
                             + logRH * ( - 1.79059E-05 + logRH * ( - 1.50345E-06 ) ) );

        return std::min( std::max( a + T * b, 0.0 ), 1.0 );

    } /* End of x_star */

    double nuclRate( double T, double x_m, double RH, double nSulf )
    {

        /* DESCRIPTION:
         * Returns the nucleation rate in #/(cm^3 s) */

        /* INPUTS:
         * - double T     :: temperature in K
         * - double x_m   :: mole fraction of sulfuric acid
         * - double RH    :: relative humidity on a ( 0, 1 ) scale
         * - double nSulf :: total gas phase concentration of sulfuric acid in molecules/cm^3
         *
         * OUTPUTS:
         * - double :: nucleation rate */
        
        /* Apply limitations */
        if ( nSulf >= 1.0E+11 )
            nSulf = 1.0E+11;

        if ( RH > 1.0E+00 )
            RH = 1;

        if ( RH < 1.0E-04 )
            RH = 1.0E-04;

        if ( T > 305.15 )
            T = 305.15;

        if ( T < 190.0 )
            T = 190.0;

        const double a = + 1.43090E-01 + T * ( + 2.21956E-00 + T * ( - 2.73911E-02 + T * 7.22811E-05 ) ) + 5.91822E-00 / x_m;
        const double b = + 1.17489E-01 + T * ( + 4.62532E-01 + T * ( - 1.18059E-02 + T * 4.04196E-05 ) ) + 1.57963E+01 / x_m;
        const double c = - 2.15554E-01 + T * ( - 8.10269E-02 + T * ( + 1.43581E-03 - T * 4.77580E-06 ) ) - 2.91297E-00 / x_m;
        const double d = - 3.58856E-00 + T * ( + 4.95080E-02 + T * ( - 2.13820E-04 + T * 3.10801E-07 ) ) - 2.93333E-02 / x_m;
        const double e = + 1.14598E-00 + T * ( - 6.00796E-01 + T * ( + 8.64245E-03 - T * 2.28947E-05 ) ) - 8.44985E-00 / x_m;
        const double f = + 2.15855E-00 + T * ( + 8.08121E-02 + T * ( - 4.07382E-04 - T * 4.01957E-07 ) ) + 7.21326E-01 / x_m;
        const double g = + 1.62410E-00 + T * ( - 1.60106E-02 + T * ( + 3.77124E-05 + T * 3.21794E-08 ) ) - 1.13255E-02 / x_m;
        const double h = + 9.71682E-00 + T * ( - 1.15048E-01 + T * ( + 1.57098E-04 + T * 4.00914E-07 ) ) + 7.11860E-01 / x_m;
        const double i = - 1.05611E-00 + T * ( + 9.03378E-03 + T * ( - 1.98417E-05 + T * 2.46048E-08 ) ) - 5.79087E-02 / x_m;
        const double j = - 1.48712E-01 + T * ( + 2.83508E-03 + T * ( - 9.24619E-06 + T * 5.00427E-09 ) ) - 1.27081E-02 / x_m;

        const double logRH    = log(RH);
        const double lognSulf = log(nSulf);

        return exp( a + logRH * ( b + logRH * ( c + logRH * d ) ) \
                   + lognSulf * ( ( e + logRH * ( f + logRH * g ) ) \
                   + lognSulf * ( h + logRH * i + lognSulf * j ) ) );

    } /* End of nuclRate */

    double nTot( double T, double x_m, double RH, double nSulf )
    {
        
        /* DESCRIPTION:
         * Returns the total number of molecules in the critical cluster */

        /* INPUTS:
         * - double T     :: temperature in K
         * - double x_m   :: mole fraction of sulfuric acid
         * - double RH    :: relative humidity on a ( 0, 1 ) scale
         * - double nSulf :: total gas phase concentration of sulfuric acid in molecules/cm^3
         *
         * OUTPUTS:
         * - double :: total number of molecules */
        
        /* Apply limitations */
        if ( nSulf >= 1.0E+11 )
            nSulf = 1.0E+11;

        if ( RH > 1.0E+00 )
            RH = 1;

        if ( RH < 1.0E-04 )
            RH = 1.0E-04;

        if ( T > 305.15 )
            T = 305.15;

        if ( T < 190.0 )
            T = 190.0;

        const double a = - 2.95413E-03 + T * ( - 9.76834E-02 + T * ( + 1.02485E-03 - T * 2.18646E-06 ) ) - 1.01717E-01 / x_m;
        const double b = - 2.05064E-03 + T * ( - 7.58504E-03 + T * ( + 1.92654E-04 - T * 6.70430E-07 ) ) - 2.55774E-01 / x_m;
        const double c = + 3.22308E-03 + T * ( + 8.52637E-04 + T * ( - 1.54757E-05 + T * 5.66661E-08 ) ) + 3.38444E-02 / x_m;
        const double d = + 4.74323E-02 + T * ( - 6.25104E-04 + T * ( + 2.65066E-06 - T * 3.67471E-09 ) ) - 2.67251E-04 / x_m;
        const double e = - 1.25211E-02 + T * ( + 5.80655E-03 + T * ( - 1.01674E-04 + T * 2.88195E-07 ) ) + 9.42243E-02 / x_m;
        const double f = - 3.85460E-02 + T * ( - 6.72316E-04 + T * ( + 2.60288E-06 + T * 1.19416E-08 ) ) - 8.51515E-03 / x_m;
        const double g = - 1.83749E-02 + T * ( + 1.72072E-04 + T * ( - 3.71766E-07 - T * 5.14875E-10 ) ) + 2.68660E-04 / x_m;
        const double h = - 6.19974E-02 + T * ( + 9.06958E-04 + T * ( - 9.11728E-07 - T * 5.36796E-09 ) ) - 7.74234E-03 / x_m;
        const double i = + 1.21827E-02 + T * ( - 1.06650E-04 + T * ( + 2.53460E-07 - T * 3.63519E-10 ) ) + 6.10065E-04 / x_m;
        const double j = + 3.20184E-04 + T * ( - 1.74762E-05 + T * ( + 6.06504E-08 - T * 1.42177E-11 ) ) + 1.35751E-04 / x_m;

        const double logRH    = log(RH);
        const double lognSulf = log(nSulf);
        
        return exp( a + logRH * ( b + logRH * ( c + logRH * d ) ) \
                   + lognSulf * ( ( e + logRH * ( f + logRH * g ) ) \
                   + lognSulf * ( h + logRH * i  + lognSulf * j ) ) );

    } /* End of nTot */

    double radCluster( double x_m, double n )
    {
        
        /* DESCRIPTION:
         * Returns the radius of the cluster in m  */

        /* INPUTS:
         * - double x_m :: mole fraction of sulfuric acid
         * - double n   :: total number of molecules in the critical cluster
         *
         * OUTPUTS:
         * - double :: cluster radius  */

        return 1.0E-09 * exp( -1.6524245E+00 + 4.2316402E-01 * x_m + 3.346648E-01 * log(n) );

    } /* End of radCluster */

    double nThresh( double T, double RH )
    {

        /* DESCRIPTION:
         * Returns the threshold concentration of sulfuric acid in molecules/cm^3 that produces a nucleation rate of 1 #/(cm^3 s)  */

        /* INPUTS:
         * - double T  :: temperature in K
         * - double RH :: relative humidity on a ( 0 - 1 ) scale 
         *
         * OUTPUTS:
         * - double :: threshold concentration  */

        /* Refer to errate for corrected formula!! */

        if ( T > 305.15 )
            T = 305.15;

        if ( T < 190.0 )
            T = 190.0;

        double invT  = 1.0 / T;
        double logRH = log(RH);

        return exp( - 2.79243E+02 + 1.17344E+01 * RH + 2.27009E+04 * invT - 1.08864E+03 * invT * RH + 1.14436E+00 * T \
                    - 3.02331E-02 * RH * T - 1.30254E-03 * T * T - 6.38697E+00 * logRH + 8.5498E+02 * invT * logRH \
                    + 8.79662E-03 * T * logRH );
    

    } /* End of nThresh */

}

/* End of Nucleation.cpp */
