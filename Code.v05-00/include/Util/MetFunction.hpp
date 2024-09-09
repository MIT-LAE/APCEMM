/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* MetFunction Header File                                          */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 11/2/2018                                 */
/* File                 : MetFunction.hpp                           */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef METFUNCTION_H_INCLUDED
#define METFUNCTION_H_INCLUDED

#include <cmath>
#include "ForwardDecl.hpp"

#define ABS(x)   ( ((x) >=  0 ) ?(x):(-x) ) 

namespace met
{

    /* Defined meteorological functions */

    /* International Standard Atmosphere (ISA) */

    const double ISA_LAPSERATE = 6.5E-03;
    const double ISA_HTS = 1.1E+04;
    const double ISA_HTP = 2.0E+04;
    const double ISA_RHO0 = 1.225E+00;
    const double ISA_P0 = 101325;
    const double ISA_T0 = 288.15;

    void ISA( const double &z, double &pressure, \
              double &temperature );
    void ISA( const Vector_1D &z, Vector_1D &pressure, \
              Vector_1D &temperature );
    void ISA( const double &z, double &pressure );
    void ISA( const Vector_1D &z, Vector_1D &pressure );
    void ISA_pAlt( double &z, const double &pressure );
    void ISA_pAlt( Vector_1D &z, const Vector_1D &pressure );


    double ComputeLapseRate( const double TEMP, const double RHi, \
                                 const double DEPTH );
                                 
    std::size_t nearestNeighbor( const Vector_1D& xq, double x);
    double linearInterp( const Vector_1D& xq, const Vector_1D& yq, double x );
    double linearInterp( double x1, double y1, double x2, double y2, double xq, bool support_extrapolation = false);

    double linInterpMetData(const Vector_1D& altitude_init, const Vector_1D& metVar_init, double altitude_query);

    double satdepth_calc( const Vector_1D& RHi, const Vector_1D& alt, int iFlight, double YLIM_DOWN);
    

    inline double rhiCorrection( double rhi, double a = 0.9779, double b = 1.635 ) {
        /* 
            @param: rhi - Relative humidity wrt ice in %
            @return: "Corrected" rhi to be biased in favor of ISSRs in %

            Source: Teoh et al (2022)
            "Aviation contrail climate effects in the North Atlantic from 2016 to 2021"
            Atmospheric Chemistry and Physics
            https://doi.org/10.5194/acp-22-10919-2022
         */

        double rhi_abs = rhi / 100.0;

        double rhi_corr = (rhi_abs / a) <= 1
                        ? rhi_abs / a
                        : std::min(std::pow(rhi_abs / a, b), 1.65);
        return rhi_corr * 100.0;
    }
    struct newXCoordsPair {
        Vector_1D x0_new;
        Vector_1D dx_new;
    };
    newXCoordsPair calcNewXCoords(const Vector_1D& dy_old, const Vector_1D& dy_new, const Vector_1D& x0_old, const Vector_1D& dx_old, int nx);
}

#endif /* METFUNCTION_H_INCLUDED */
