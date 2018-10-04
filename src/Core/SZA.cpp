/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* SZA Program File                                                 */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : SZA.cpp                                   */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <cmath>

void SZA( double latitude_deg, int dayGMT,\
          double &sunRise, double &sunSet,\
          double &SZASINLAT, double &SZACOSLAT,\
          double &SZASINDEC, double &SZACOSDEC )
{
    double const A0 = 0.006918;
    double const A1 = 0.399912;
    double const A2 = 0.006758;
    double const A3 = 0.002697;
    double const B1 = 0.070257;
    double const B2 = 0.000907;
    double const B3 = 0.000148;

    const double PI = 3.141592653589793238460; /* \pi */

    double r_SZA = 2*PI*(std::floor(dayGMT) - 1)/365.0;

    double DEC = A0 - A1*cos(1*r_SZA) + B1*sin(1*r_SZA)\
                    - A2*cos(2*r_SZA) + B2*sin(2*r_SZA)\
                    - A3*cos(3*r_SZA) + B3*sin(3*r_SZA);

    SZASINLAT = std::sin(latitude_deg*PI/180);
    SZACOSLAT = std::cos(latitude_deg*PI/180);
    SZASINDEC = std::sin(DEC);
    SZACOSDEC = std::cos(DEC);

    sunRise = std::max((12.0 - 180.0/(PI*15.0)*acos(-(SZASINLAT * SZASINDEC)\
                                              / (SZACOSLAT * SZACOSDEC))), 0.0);
    sunSet  = std::min((12.0 + 180.0/(PI*15.0)*acos(-(SZASINLAT * SZASINDEC)\
                                              / (SZACOSLAT * SZACOSDEC))), 24.0);

} /* End of SZA */


