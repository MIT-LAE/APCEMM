/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* DESCRIPTION: KPP module for heterogeneous chemistry */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "KPP/KPP_Parameters.h"
#include "KPP/KPP_Global.h"
#include "KPP/KPP_Sparse.h"
#include "KPP/KPP.hpp"
#include "Core/Parameters.hpp"
#include "Util/Error.hpp"

#define MAX(a,b) ( ((a) >= (b)) ? (a):(b)  )
#define MIN(b,c) ( ((b) < (c))  ? (b):(c)  )
 
static const double PSCMINLIFE = 1.0E-03;
static const double GAMMA_HO2 = 0.2; 

double ARSL1K( const double AREA,  const double RADI,  const double AIRDENS, \
               const double STKCF, const double XTEMP, const double SQM );
void CHECK_NAT( bool &IS_NAT, bool &IS_PSC,   bool &IS_STRAT,    \
                const unsigned int STATE_PSC, const double PATM, double tropopausePressure );
double N2O5( unsigned int N, const double TEMP, const double RH );
double CLD1K_BrNO3( const double AIRDENS, const double TEMP, \
                    const double QL,      const double CLDF, \
                    bool IS_LAND,         bool IS_ICE );
double CLDICE_HBrHOBr( const double AIRDENS,  const double TEMP,        \
                       const double QI,       const double CLDF,        \
                       const double HBr,      const double HOBr,        \
                       double &K_HBr,         double &K_HOBr,           \
                       const double AREA_ICE, const double EFFRADI_ICE, \
                       const double IWC_ICE );
double HETNO3( const double A, const double B, const double AREA[NAERO], \
               const double RADI[NAERO],       const double TEMP,        \
               const double AIRDENS );
double HETNO2( const double A, const double B, const double AREA[NAERO], \
               const double RADI[NAERO],       const double TEMP,        \
               const double AIRDENS );
double HETHO2( const double A, const double B, const double AREA[NAERO], \
               const double RADI[NAERO],       const double TEMP,        \
               const double AIRDENS );
double HETHBr( const double A, const double B, const double KHETI_SLA[11], \
               const double AREA[NAERO],       const double RADI[NAERO],   \
               const double TEMP,              const double AIRDENS,       \
               bool IS_STRAT );
double HETN2O5( const double A, const double B,  const double KHETI_SLA[11], \
                const double AREA[NAERO],        const double RADI[NAERO],   \
                const double TEMP,               const double AIRDENS,       \
                const double RH, double SPC_SO4, double SPC_NIT,             \
                bool NATSURFACE );
double HETBrNO3( const double A, const double B, const double KHETI_SLA[11], \
                 const double AREA[NAERO],       const double RADI[NAERO],   \
                 const double TEMP,              const double AIRDENS,       \
                 bool IS_STRAT, bool IS_PSC,     double CLD_BrNO3_RC,        \
                 bool NATSURFACE );
double HETHOBr( const double A, const double B, const double KHETI_SLA[11], \
                const double AREA[NAERO],       const double RADI[NAERO],   \
                const double TEMP,              const double AIRDENS,       \
                bool IS_STRAT );
double HETHOBr_ice( );
double HETHBr_ice( );
double HETN2O5_PSC( const double A, const double B, const double KHETI_SLA[11], \
                    const double AREA[NAERO],       const double RADI[NAERO],   \
                    const double TEMP,              const double AIRDENS,       \
                    bool IS_STRAT,                  bool NATSURFACE );
double HETClNO3_PSC1( const double A, const double B, const double KHETI_SLA[11], \
                      const double AREA[NAERO],       const double RADI[NAERO],   \
                      const double TEMP,              const double AIRDENS,       \
                      bool IS_STRAT,                  bool NATSURFACE );
double HETClNO3_PSC2( const double A, const double B, const double KHETI_SLA[11], \
                      const double AREA[NAERO],       const double RADI[NAERO],   \
                      const double TEMP,              const double AIRDENS,       \
                      bool IS_STRAT,                  bool NATSURFACE );
double HETClNO3_PSC3( const double A, const double B, const double KHETI_SLA[11], \
                      const double AREA[NAERO],       const double RADI[NAERO],   \
                      const double TEMP,              const double AIRDENS,       \
                      bool IS_STRAT,                  bool NATSURFACE );
double HETBrNO3_PSC( const double A, const double B, const double KHETI_SLA[11], \
                     const double AREA[NAERO],       const double RADI[NAERO],   \
                     const double TEMP,              const double AIRDENS,       \
                     bool IS_STRAT,                  bool NATSURFACE );
double HETHOCl_PSC1( const double A, const double B, const double KHETI_SLA[11], \
                     const double AREA[NAERO],       const double RADI[NAERO],   \
                     const double TEMP,              const double AIRDENS,       \
                     bool IS_STRAT,                  bool NATSURFACE );
double HETHOCl_PSC2( const double A, const double B, const double KHETI_SLA[11], \
                     const double AREA[NAERO],       const double RADI[NAERO],   \
                     const double TEMP,              const double AIRDENS,       \
                     bool IS_STRAT,                  bool NATSURFACE );
double HETHOBr_PSC( const double A, const double B, const double KHETI_SLA[11], \
                    const double AREA[NAERO],       const double RADI[NAERO],   \
                    const double TEMP,              const double AIRDENS,       \
                    bool IS_STRAT,                  bool NATSURFACE );


void GC_SETHET( const double TEMP, const double PATM, const double AIRDENS, \
                const double RELHUM, const unsigned int STATE_PSC,          \
                const double SPC[], const double AREA[NAERO],               \
                const double RADI[NAERO], const double IWC,                 \
                const double KHETI_SLA[11], double tropopausePressure )
{

    /* Sets up the array of heterogeneous chemistry rates for the KPP chemistry solver */

    /* INPUT PARAMETERS:
     *
     * const double TEMP          : Temperature in K 
     * const double PATM          : Pressure in Pa 
     * const double AIRDENS       : Air density in molec/cm^3 
     * const double RELHUM        : Relative humidity [0-1] 
     * const double SPC[]         : Species array in molec/cm^3 
     * const double AREA[NAERO]   : Aerosol area in m^2/cm^3
     * const double RADI[NAERO]   : Aerosol radius in m 
     * const double IWC           : Ice water content in kg/cm^3
     * const double KHETI_SLA[11] : Sticking coefficients */
     
    /* Aerosol list:
     * 0 : NAT/ice 
     * 1 : liquid stratospheric aerosol
     * 2 : sulfate aerosols (near surface)
     * 3 : soot */

    /* Scalars */
    double ADJUSTEDRATE, HBr_RTEMP, HOBr_RTEMP, \
           QICE,         QLIQ,      CLDF,                  \
           SPC_BrNO3,    SPC_ClNO3, SPC_H2O,   SPC_HBr,    \
           SPC_HCl,      SPC_HOBr,  SPC_HOCl,  SPC_N2O5,   \
           SPC_SO4,      SPC_NIT,   CLD_BrNO3_RC,          \
           KI_HBr,       KI_HOBr;

    /* UCX-based mechanisms */
    unsigned int PSCIDX;
    double EDUCTCONC, LIMITCONC;
    double PSCEDUCTCONC[11][2];
   
    /* Logicals */
    bool SAFEDIV, PSCBOX, STRATBOX, NATSURFACE;
    bool IS_LAND, IS_ICE;

    /* GC_SETHET begins here! */

    /* Zero scalars and arrays */
    ADJUSTEDRATE = 0.0E+00;
    HBr_RTEMP    = 0.0E+00;
    HOBr_RTEMP   = 0.0E+00;
    QICE         = 0.0E+00;
    QLIQ         = 0.0E+00;
    CLDF         = 0.0E+00;
    SPC_BrNO3    = 0.0E+00;
    SPC_ClNO3    = 0.0E+00;
    SPC_H2O      = 0.0E+00;
    SPC_HBr      = 0.0E+00;
    SPC_HCl      = 0.0E+00;
    SPC_HOBr     = 0.0E+00;
    SPC_HOCl     = 0.0E+00;
    SPC_N2O5     = 0.0E+00;
   
    /* Zero variables for UCX */
    EDUCTCONC    = 0.0E+00;
    LIMITCONC    = 0.0E+00;
    for ( PSCIDX = 0; PSCIDX < 11; PSCIDX++ ) {
        PSCEDUCTCONC[PSCIDX][0] = 0.0E+00;
        PSCEDUCTCONC[PSCIDX][1] = 0.0E+00;
    }
    
    /* Initialize logicals */
    SAFEDIV    = 0;
    PSCBOX     = 0;
    STRATBOX   = 0;
    NATSURFACE = 0;

    IS_LAND    = 0;
    IS_ICE     = 0;

    /* Get species concentrations [molec/cm^3] */
#ifndef ind_NIT
    SPC_NIT = 0.0E+00;
#else
    SPC_NIT = SPC[ind_NIT];
#endif

#ifndef ind_SO4
    SPC_SO4 = 0.0E+00;
#else
    SPC_SO4 = SPC[ind_SO4];
#endif

#ifndef ind_HBr
    SPC_HBr = 0.0E+00;
#else
    SPC_HBr = SPC[ind_HBr];
#endif

#ifndef ind_HOBr
    SPC_HOBr = 0.0E+00;
#else
    SPC_HOBr = SPC[ind_HOBr];
#endif

#ifndef ind_N2O5
    SPC_N2O5 = 0.0E+00;
#else
    SPC_N2O5 = SPC[ind_N2O5];
#endif

#ifndef ind_H2O
    SPC_H2O = 0.0E+00;
#else
    SPC_H2O = SPC[ind_H2O];
#endif

#ifndef ind_HCl
    SPC_HCl = 0.0E+00;
#else
    SPC_HCl = SPC[ind_HCl];
#endif

#ifndef ind_ClNO3
    SPC_ClNO3 = 0.0E+00;
#else
    SPC_ClNO3 = SPC[ind_ClNO3];
#endif

#ifndef ind_BrNO3
    SPC_BrNO3 = 0.0E+00;
#else
    SPC_BrNO3 = SPC[ind_BrNO3];
#endif

#ifndef ind_HOCl
    SPC_HOCl = 0.0E+00;
#else
    SPC_HOCl = SPC[ind_HOCl];
#endif

    PSCEDUCTCONC[ 0][0] = SPC_N2O5;
    PSCEDUCTCONC[ 0][1] = SPC_H2O;

    PSCEDUCTCONC[ 1][0] = SPC_N2O5;
    PSCEDUCTCONC[ 1][1] = SPC_HCl;

    PSCEDUCTCONC[ 2][0] = SPC_ClNO3;
    PSCEDUCTCONC[ 2][1] = SPC_H2O;

    PSCEDUCTCONC[ 3][0] = SPC_ClNO3;
    PSCEDUCTCONC[ 3][1] = SPC_HCl;

    PSCEDUCTCONC[ 4][0] = SPC_ClNO3;
    PSCEDUCTCONC[ 4][1] = SPC_HBr;

    PSCEDUCTCONC[ 5][0] = SPC_BrNO3;
    PSCEDUCTCONC[ 5][1] = SPC_H2O;

    PSCEDUCTCONC[ 6][0] = SPC_BrNO3;
    PSCEDUCTCONC[ 6][1] = SPC_HCl;

    PSCEDUCTCONC[ 7][0] = SPC_HOCl;
    PSCEDUCTCONC[ 7][1] = SPC_HCl;

    PSCEDUCTCONC[ 8][0] = SPC_HOCl;
    PSCEDUCTCONC[ 8][1] = SPC_HBr;

    PSCEDUCTCONC[ 9][0] = SPC_HOBr;
    PSCEDUCTCONC[ 9][1] = SPC_HCl;

    /* This is still pseudo-first order - ignore */
    PSCEDUCTCONC[10][0] = SPC_HOBr;
    PSCEDUCTCONC[10][1] = SPC_HBr;


    if ( TEMP <= 248.0 ) 
        QLIQ = 0.0E+00;
    else if ( TEMP >= 268.0 )
        QLIQ = 1.0E-06;
    else 
        QLIQ = 1.0E-06 * ( TEMP - 248.0 ) / double(20.0);

    QICE = 1.0E-06 - QLIQ;


    /* Check surface type of PSCs */
    CHECK_NAT( NATSURFACE, PSCBOX, STRATBOX, STATE_PSC, PATM, tropopausePressure );

    if ( !PSCBOX )
        CLD_BrNO3_RC = CLD1K_BrNO3( AIRDENS, TEMP, QLIQ, CLDF, IS_LAND, IS_ICE );

    if ( !PSCBOX )
        CLDICE_HBrHOBr( AIRDENS, TEMP, QICE, CLDF, SPC_HBr, SPC_HOBr, KI_HBr, KI_HOBr, AREA[0], RADI[0], IWC );
    else {
        KI_HBr  = 0.0E+00;
        KI_HOBr = 0.0E+00;
    }

    for ( PSCIDX = 0; PSCIDX < 11; PSCIDX++ ) {
        HET[PSCIDX][0] = 0.0E+00;
        HET[PSCIDX][1] = 0.0E+00;
        HET[PSCIDX][2] = 0.0E+00;
    }

    /* Calculate and pass het rates to the KPP rate array */
    HET[ind_HO2][0]   = HETHO2(        3.30E+01, 2.00E-01, AREA, RADI, TEMP, AIRDENS);
    HET[ind_NO2][0]   = HETNO2(        4.60E+01, 1.00E-04, AREA, RADI, TEMP, AIRDENS); 
    HET[ind_NO3][0]   = HETNO3(        6.20E+01, 1.00E-01, AREA, RADI, TEMP, AIRDENS);
    HET[ind_N2O5][0]  = HETN2O5(       1.08E+02, 1.00E-01, KHETI_SLA, AREA, RADI, TEMP, AIRDENS, RELHUM, SPC_SO4, SPC_NIT, NATSURFACE);
    HET[ind_BrNO3][0] = HETBrNO3(      1.42E+02, 3.00E-01, KHETI_SLA, AREA, RADI, TEMP, AIRDENS, STRATBOX, PSCBOX, CLD_BrNO3_RC, NATSURFACE); 
    HET[ind_HOBr][0]  = HETHOBr(       0.97E+02, 2.00E-01, KHETI_SLA, AREA, RADI, TEMP, AIRDENS, STRATBOX); 
    HET[ind_HBr][0]   = HETHBr(        0.81E+02, 2.00E-01, KHETI_SLA, AREA, RADI, TEMP, AIRDENS, STRATBOX); 
    HET[ind_HOBr][1]  = HETHOBr_ice( ); 
    HET[ind_HBr][1]   = HETHBr_ice( ); 
    HET[ind_N2O5][1]  = HETN2O5_PSC(   1.08E+02, 0.00E+00, KHETI_SLA, AREA, RADI, TEMP, AIRDENS, STRATBOX, NATSURFACE);
    HET[ind_ClNO3][0] = HETClNO3_PSC1( 0.97E+02, 0.00E+00, KHETI_SLA, AREA, RADI, TEMP, AIRDENS, STRATBOX, NATSURFACE);
    HET[ind_ClNO3][1] = HETClNO3_PSC2( 0.97E+02, 0.00E+00, KHETI_SLA, AREA, RADI, TEMP, AIRDENS, STRATBOX, NATSURFACE); 
    HET[ind_ClNO3][2] = HETClNO3_PSC3( 0.97E+02, 0.00E+00, KHETI_SLA, AREA, RADI, TEMP, AIRDENS, STRATBOX, NATSURFACE); 
    HET[ind_BrNO3][1] = HETBrNO3_PSC(  1.42E+02, 0.00E+00, KHETI_SLA, AREA, RADI, TEMP, AIRDENS, STRATBOX, NATSURFACE); 
    HET[ind_HOCl][0]  = HETHOCl_PSC1(  0.52E+02, 0.00E+00, KHETI_SLA, AREA, RADI, TEMP, AIRDENS, STRATBOX, NATSURFACE); 
    HET[ind_HOCl][1]  = HETHOCl_PSC2(  0.52E+02, 0.00E+00, KHETI_SLA, AREA, RADI, TEMP, AIRDENS, STRATBOX, NATSURFACE); 
    HET[ind_HOBr][2]  = HETHOBr_PSC(   0.97E+02, 0.00E+00, KHETI_SLA, AREA, RADI, TEMP, AIRDENS, STRATBOX, NATSURFACE); 

    /* Kludging the rates to be equal to one another to avoid having
     * to keep setting equality in solver */
    if ( ( HET[ind_HBr][0] > 0 ) && ( HET[ind_HOBr][0] > 0 ) ) {

        /* Select the minimum of the two rates */
        HBr_RTEMP  = HET[ind_HBr][0]  * SPC_HBr;
        HOBr_RTEMP = HET[ind_HOBr][0] * SPC_HOBr;

        /* If HBr rate is the largest */
        if ( HBr_RTEMP > HOBr_RTEMP ) {

            SAFEDIV = SafeDiv( HET[ind_HOBr][0] * SPC_HOBr, SPC_HBr );

            if ( SAFEDIV ) {
                HET[ind_HBr][0] = HET[ind_HOBr][0] * SPC_HOBr / SPC_HBr;
            } else {
                /* If division is not possible, set the rates to really small values */
                HET[ind_HBr][0]  = 1.0E-40;
                HET[ind_HOBr][0] = 1.0E-40;
            }

        } else {

            SAFEDIV = SafeDiv( HET[ind_HBr][0] * SPC_HBr, SPC_HOBr );
            
            if ( SAFEDIV ) {
                HET[ind_HOBr][0] = HET[ind_HBr][0] * SPC_HBr / SPC_HOBr;
            } else {
                /* If division is not possible, set the rates to really small values */
                HET[ind_HBr][0]  = 1.0E-40;
                HET[ind_HOBr][0] = 1.0E-40;
            }

        }

    }

    /* Limit rates to prevent solver failure for PSC heterogeneous chemistry */

    for ( PSCIDX = 0; PSCIDX < 10; PSCIDX++ ) {

        /* Pseudo-first-order reactions - divide by number concentrations
         * of aerosol-phase educt to yield 2nd-order constant */
        EDUCTCONC = PSCEDUCTCONC[PSCIDX][1];
        LIMITCONC = PSCEDUCTCONC[PSCIDX][0];

        /* Initialize adjusted rates */
        if ( PSCIDX == 0 ) {
            /* N2O5 + H2O */
            ADJUSTEDRATE = HET[ind_N2O5][0];
        } else if ( PSCIDX == 1 ) {
            /* N2O5 + HCl */
            ADJUSTEDRATE = HET[ind_N2O5][1];
        } else if ( PSCIDX == 2 ) {
            /* ClNO3 + H2O */
            ADJUSTEDRATE = HET[ind_ClNO3][0];
        } else if ( PSCIDX == 3 ) {
            /* ClNO3 + HCl */
            ADJUSTEDRATE = HET[ind_ClNO3][1];
        } else if ( PSCIDX == 4 ) {
            /* ClNO3 + HBr */
            ADJUSTEDRATE = HET[ind_ClNO3][2];
        } else if ( PSCIDX == 5 ) {
            /* BrNO3 + H2O */
            ADJUSTEDRATE = HET[ind_BrNO3][0];
        } else if ( PSCIDX == 6 ) {
            /* BrNO3 + HCl */
            ADJUSTEDRATE = HET[ind_BrNO3][1];
        } else if ( PSCIDX == 7 ) {
            /* HOCl + HCl */
            ADJUSTEDRATE = HET[ind_HOCl][0];
        } else if ( PSCIDX == 8 ) {
            /* HOCl + HBr */
            ADJUSTEDRATE = HET[ind_HOCl][1];
        } else if ( PSCIDX == 9 ) {
            /* HOBr + HCl */
            ADJUSTEDRATE = HET[ind_HOBr][2];
        }

        /* ---SAFETY-CHECK REACTION---
         * Definition of 2nd order reaction rate:
         * k[A][B] = -d[A]/dt = -d[B]/dt
         *
         * However, here we are using a pseudo-first order
         * reaction rate, ki, and assuming that [B] is
         * abundant. To get k, we will therefore perform:
         * k = ki/[B]
         *
         * This will yield the following when solved:
         * -d[A]/dt = ki[A] = -d[B]/dt
         *
         * This has some problems, especially for small [B]*
         * To get around this, we run the following checks:
         *
         * 1. The lifetime of [A] is 1/ki. If this is below
         *    PSCMINLIFE, limit reaction rate to yield the
         *    specified lifetimedepletion (ki = 1/60)
         * 2. The depletion time of [B] is [B]/(ki[A]). If
         *    this is below PSCMINLIFE, limit reaction rate
         *    (ki = [B]/(T*[A])
         * 3. If [B] is < 100 molec/cm3, or ki/[B] yields
         *    a Nan, force k = 0.
         *
         * If all these checks are passed, we set k = ki/[B].
         * Rxn 11 is first-order - ignore */
        
        if ( PSCIDX == 0 ) {
            /* Convert from 1st-order to 2nd-order */
            SAFEDIV = SafeDiv( EDUCTCONC, LIMITCONC );

            if ( SAFEDIV ) {
                /* Temporarily store [B]/(T*[A]) */
                LIMITCONC = EDUCTCONC/(PSCMINLIFE*LIMITCONC);

                if ( ADJUSTEDRATE > LIMITCONC )
                    ADJUSTEDRATE = LIMITCONC;

            } else {
                ADJUSTEDRATE = 0.0E+00;
            }

            SAFEDIV = SafeDiv( ADJUSTEDRATE, EDUCTCONC );

            if ( ( EDUCTCONC > 1.0E+02 ) && ( SAFEDIV ) ) {
                ADJUSTEDRATE /= EDUCTCONC;
            } else {
                ADJUSTEDRATE = 0.0E+00;
            }

        } else if ( PSCIDX != 10 ) {
            if ( ADJUSTEDRATE > ( 1.0E+00/PSCMINLIFE ) )
                ADJUSTEDRATE = 1.0E+00/PSCMINLIFE;

            /* Convert from 1st-order to 2nd-order */
            SAFEDIV = SafeDiv( EDUCTCONC, LIMITCONC );

            if ( SAFEDIV ) {
                /* Temporarily store [B]/(T*[A]) */
                LIMITCONC = EDUCTCONC/(PSCMINLIFE*LIMITCONC);

                if ( ADJUSTEDRATE > LIMITCONC )
                    ADJUSTEDRATE = LIMITCONC;

            } else {
                ADJUSTEDRATE = 0.0E+00;
            }
            
            SAFEDIV = SafeDiv( ADJUSTEDRATE, EDUCTCONC );
            
            if ( ( EDUCTCONC > 1.0E+02 ) && ( SAFEDIV ) ) {
                ADJUSTEDRATE /= EDUCTCONC;
            } else {
                ADJUSTEDRATE = 0.0E+00;
            }

        }

        /* Copy adjusted rates to HET */
        /* Initialize adjusted rates */
        if ( PSCIDX == 0 ) {
            /* N2O5 + H2O */
            HET[ind_N2O5][0] = ADJUSTEDRATE;
        } else if ( PSCIDX == 1 ) {
            /* N2O5 + HCl */
            HET[ind_N2O5][1] = ADJUSTEDRATE;
        } else if ( PSCIDX == 2 ) {
            /* ClNO3 + H2O */
            HET[ind_ClNO3][0] = ADJUSTEDRATE;
        } else if ( PSCIDX == 3 ) {
            /* ClNO3 + HCl */
            HET[ind_ClNO3][1] = ADJUSTEDRATE;
        } else if ( PSCIDX == 4 ) {
            /* ClNO3 + HBr */
            HET[ind_ClNO3][2] = ADJUSTEDRATE;
        } else if ( PSCIDX == 5 ) {
            /* BrNO3 + H2O */
            HET[ind_BrNO3][0] = ADJUSTEDRATE;
        } else if ( PSCIDX == 6 ) {
            /* BrNO3 + HCl */
            HET[ind_BrNO3][1] = ADJUSTEDRATE;
        } else if ( PSCIDX == 7 ) {
            /* HOCl + HCl */
            HET[ind_HOCl][0] = ADJUSTEDRATE;
        } else if ( PSCIDX == 8 ) {
            /* HOCl + HBr */
            HET[ind_HOCl][1] = ADJUSTEDRATE;
        } else if ( PSCIDX == 9 ) {
            /* HOBr + HCl */
            HET[ind_HOBr][2] = ADJUSTEDRATE;
        }

    }


} /* End of GC_SETHET */

void CHECK_NAT( bool &IS_NAT, bool &IS_PSC, bool &IS_STRAT, \
                const unsigned int STATE_PSC, const double PATM, double tropopausePressure )
{

    /* This function determines whether the solid PSC is composed of ice
     * or NAT (nitric acid trihydrate) (needed for heterogeneous chemistry),
     * or indeed if there is any direct PSC calculation at all */

    bool IS_TROP;

    /* Check if box is in the troposphere */
    IS_TROP  = ( PATM > tropopausePressure );

    /* Check if box is in the stratosphere */
    IS_STRAT = ( !IS_TROP );

    /* Check if there are solid PSCs */
    IS_PSC   = ( ( PSC ) && ( IS_STRAT ) && ( STATE_PSC >= 2.0 ) );

    /* Check if there is surface NAT */
#ifdef SPC_NIT
    IS_NAT   = ( ( IS_PSC ) && ( SPC_NIT > 0 ) );
#else
    IS_NAT   = 0;
#endif

} /* End of CHECK_NAT */


double HETNO3( const double A, const double B, const double AREA[NAERO], const double RADI[NAERO], const double TEMP, const double AIRDENS )
{

    /* DESCRIPTION: Set the heterogeneous chemistry rate for NO3 */

    /* INPUTS:
     * double A           : Molar weight in g/mol
     * double B           : 
     * double AREA[NAERO] : Aerosol surface area in m^2/cm^3 
     * double RADI[NAERO] : Aerosol radius in m 
     * double AIRDENS     : Air density in molec/cm^3 */
    
    bool DO_EDUCT;
    unsigned int N;
    double XSTKCF, ADJUSTEDRATE;
    double HET_NO3;

    /* Initialize */
    HET_NO3 = 0.0E+00;
    ADJUSTEDRATE = 0.0E+00;
    XSTKCF = 0.0E+00;

    /* Don't do PSC rate adjustment */
    DO_EDUCT = 0;

    /* Loop over aerosol types */
    for ( N = 0; N < NAERO; N++ ) {

        XSTKCF = B;

        if ( N == 1 ) {
            /* Calculate for stratospheric liquid aerosol */
            ADJUSTEDRATE = AREA[N] * XSTKCF;
        } else {
            /* Reaction rate for surface of aerosol */
            ADJUSTEDRATE = ARSL1K( AREA[N], RADI[N], AIRDENS, XSTKCF, pow( TEMP, 0.5 ), pow( A, 0.5 ) );
        }

        if ( ( DO_EDUCT ) && ( N < 2 ) ) {
            /* PSC reaction - prevent excessive reaction rate */
            if ( ADJUSTEDRATE > 1.0E+00 / PSCMINLIFE )
                ADJUSTEDRATE = 1.0E+00 / PSCMINLIFE;

        }

        /* Add to overall reaction rate */
        HET_NO3 += ADJUSTEDRATE;

    }

    return HET_NO3;

} /* End of HETNO3 */

double HETNO2( const double A, const double B, const double AREA[NAERO], const double RADI[NAERO], const double TEMP, const double AIRDENS )
{

    /* DESCRIPTION: Set the heterogeneous chemistry rate for NO2 */

    /* INPUTS:
     * double A           : Molar weight in g/mol
     * double B           : 
     * double AREA[NAERO] : Aerosol surface area in m^2/cm^3 
     * double RADI[NAERO] : Aerosol radius in m 
     * double AIRDENS     : Air density in molec/cm^3 */
    
    bool DO_EDUCT;
    unsigned int N;
    double XSTKCF, ADJUSTEDRATE;
    double HET_NO2;

    /* Initialize */
    HET_NO2 = 0.0E+00;
    ADJUSTEDRATE = 0.0E+00;
    XSTKCF = 0.0E+00;

    /* Don't do PSC rate adjustment */
    DO_EDUCT = 0;

    /* Loop over aerosol types */
    for ( N = 0; N < NAERO; N++ ) {

        XSTKCF = B;

        if ( N == 1 ) {
            /* Calculate for stratospheric liquid aerosol */
            ADJUSTEDRATE = AREA[N] * XSTKCF;
        } else {
            /* Reaction rate for surface of aerosol */
            ADJUSTEDRATE = ARSL1K( AREA[N], RADI[N], AIRDENS, XSTKCF, pow( TEMP, 0.5), pow( A, 0.5 ) );
        }

        if ( ( DO_EDUCT ) && ( N < 2 ) ) {
            /* PSC reaction - prevent excessive reaction rate */
            if ( ADJUSTEDRATE > 1.0E+00 / PSCMINLIFE )
                ADJUSTEDRATE = 1.0E+00 / PSCMINLIFE;

        }

        /* Add to overall reaction rate */
        HET_NO2 += ADJUSTEDRATE;

    }

    return HET_NO2;

} /* End of HETNO2 */

double HETHO2( const double A, const double B, const double AREA[NAERO], const double RADI[NAERO], const double TEMP, const double AIRDENS )
{

    /* DESCRIPTION: Set the heterogeneous chemistry rate for HO2 */

    /* INPUTS:
     * double A           : Molar weight in g/mol
     * double B           : 
     * double AREA[NAERO] : Aerosol surface area in m^2/cm^3 
     * double RADI[NAERO] : Aerosol radius in m 
     * double AIRDENS     : Air density in molec/cm^3 */
    
    bool DO_EDUCT;
    unsigned int N;
    double XSTKCF, ADJUSTEDRATE;
    double HET_HO2;

    /* Initialize */
    HET_HO2 = 0.0E+00;
    ADJUSTEDRATE = 0.0E+00;
    XSTKCF = 0.0E+00;

    /* Don't do PSC rate adjustment */
    DO_EDUCT = 0;

    /* Loop over aerosol types */
    for ( N = 0; N < NAERO; N++ ) {

        if ( N < 2 ) 
            XSTKCF = 1.0E-30;
        else
            XSTKCF = GAMMA_HO2;

        if ( N == 1 ) {
            /* Calculate for stratospheric liquid aerosol */
            ADJUSTEDRATE = AREA[N] * XSTKCF;
        } else {
            /* Reaction rate for surface of aerosol */
            ADJUSTEDRATE = ARSL1K( AREA[N], RADI[N], AIRDENS, XSTKCF, pow( TEMP, 0.5 ), pow( A, 0.5 ) );
        }

        if ( ( DO_EDUCT ) && ( N < 2 ) ) {
            /* PSC reaction - prevent excessive reaction rate */
            if ( ADJUSTEDRATE > 1.0E+00 / PSCMINLIFE )
                ADJUSTEDRATE = 1.0E+00 / PSCMINLIFE;

        }

        /* Add to overall reaction rate */
        HET_HO2 += ADJUSTEDRATE;

    }

    return HET_HO2;

} /* End of HETHO2 */

double HETHBr( const double A, const double B, const double KHETI_SLA[11], const double AREA[NAERO], const double RADI[NAERO], const double TEMP, const double AIRDENS, bool IS_STRAT )
{

    /* DESCRIPTION: Set the heterogeneous chemistry rate for HBr */

    /* INPUTS:
     * double A             : Molar weight in g/mol
     * double B             : 
     * double KHETI_SLA[11] : Sticking coefficients 
     * double AREA[NAERO]   : Aerosol surface area in m^2/cm^3 
     * double RADI[NAERO]   : Aerosol radius in m 
     * double AIRDENS       : Air density in molec/cm^3 
     * bool IS_STRAT        : In stratosphere? */
    
    bool DO_EDUCT;
    unsigned int N;
    double XSTKCF, ADJUSTEDRATE;
    double HET_HBr;

    /* Initialize */
    HET_HBr = 0.0E+00;
    ADJUSTEDRATE = 0.0E+00;
    XSTKCF = 0.0E+00;

    /* Only apply PSC rate adjustment if at high altitude */
    DO_EDUCT = IS_STRAT;

    /* Loop over aerosol types */
    for ( N = 0; N < NAERO; N++ ) {

        if ( N == 2 ) 
            XSTKCF = B;
        else if ( N == 1 )
            XSTKCF = KHETI_SLA[10];
        else
            XSTKCF = 1.00E-01;

        if ( N == 1 ) {
            /* Calculate for stratospheric liquid aerosol */
            ADJUSTEDRATE = AREA[N] * XSTKCF;
        } else {
            /* Reaction rate for surface of aerosol */
            ADJUSTEDRATE = ARSL1K( AREA[N], RADI[N], AIRDENS, XSTKCF, pow( TEMP, 0.5 ), pow( A, 0.5 ) );
        }

        if ( ( DO_EDUCT ) && ( N < 2 ) ) {
            /* PSC reaction - prevent excessive reaction rate */
            if ( ADJUSTEDRATE > 1.0E+00 / PSCMINLIFE )
                ADJUSTEDRATE = 1.0E+00 / PSCMINLIFE;

        }

        /* Add to overall reaction rate */
        HET_HBr += ADJUSTEDRATE;

    }

    return HET_HBr;

} /* End of HETHBr */

double HETN2O5( const double A, const double B, const double KHETI_SLA[11], const double AREA[NAERO], const double RADI[NAERO], const double TEMP, const double AIRDENS, \
                const double RH, double SPC_SO4, double SPC_NIT, bool NATSURFACE )
{

    /* DESCRIPTION: Set the heterogeneous chemistry rate for N2O5 */

    /* INPUTS:
     * double A             : Molar weight in g/mol
     * double B             : 
     * double KHETI_SLA[11] : Sticking coefficients 
     * double AREA[NAERO]   : Aerosol surface area in m^2/cm^3 
     * double RADI[NAERO]   : Aerosol radius in m 
     * double AIRDENS       : Air density in molec/cm^3 
     * double SPC_SO4       : SO4 concentration molec/cm^3
     * double SPC_NIT       : NIT concentration molec/cm^3
     * bool NATSURFACE      : Frozen HNO3? */
    
    bool DO_EDUCT;
    unsigned int N;
    double XSTKCF, ADJUSTEDRATE;
    double TMP1, TMP2;
    double HET_N2O5;

    /* Initialize */
    HET_N2O5 = 0.0E+00;
    ADJUSTEDRATE = 0.0E+00;
    XSTKCF = 0.0E+00;
    TMP1 = 0.0E+00;
    TMP2 = 0.0E+00;

    /* Always apply PSC rate adjustment*/
    DO_EDUCT = 1;

    /* Loop over aerosol types */
    for ( N = 0; N < NAERO; N++ ) {

        /* Get gamma for N2O5 hydrolysis, whic is a function of
         * aerosol type, temperature, and relative humidity */
        if ( N == 0 ) {
            if ( NATSURFACE ) 
                XSTKCF = 4.00E-04;
            else
                XSTKCF = 2.00E-02;
        } else if ( N == 1 ) {
            XSTKCF = KHETI_SLA[ 0];
        } else {
            XSTKCF = N2O5( N, TEMP, RH );
        }

        if ( N == 2 ) {
            TMP1 = SPC_NIT + SPC_SO4;
            TMP2 = SPC_NIT;
            if ( TMP1 > 0.0 )
                XSTKCF *= ( 1.0E+00 - 9.0E-01 * TMP2 / TMP1 );
        }
        if ( N == 1 ) {
            ADJUSTEDRATE = AREA[N] * XSTKCF;
        } else {
            ADJUSTEDRATE = ARSL1K( AREA[N], RADI[N], AIRDENS, XSTKCF, pow( TEMP, 0.5 ), pow( A, 0.5 ) );
        }

        if ( ( DO_EDUCT ) && ( N < 2 ) ) {
            /* PSC Reaction - prevent excessive reaction rate */
            if ( ADJUSTEDRATE > 1.0E+00/PSCMINLIFE )
                ADJUSTEDRATE = 1.0E+00/PSCMINLIFE;
        }

        HET_N2O5 += ADJUSTEDRATE;

    }

    return HET_N2O5;

} /* End of HETN2O5 */

double HETBrNO3( const double A, const double B, const double KHETI_SLA[11], const double AREA[NAERO], const double RADI[NAERO], const double TEMP,  const double AIRDENS, bool IS_STRAT, bool IS_PSC, double CLD_BrNO3_RC, bool NATSURFACE )
{

    /* DESCRIPTION: Set the heterogeneous chemistry rate for BrNO3 */

    /* INPUTS:
     * double A             : Molar weight in g/mol
     * double B             : 
     * double KHETI_SLA[11] : Sticking coefficients 
     * double AREA[NAERO]   : Aerosol surface area in m^2/cm^3 
     * double RADI[NAERO]   : Aerosol radius in m 
     * double AIRDENS       : Air density in molec/cm^3 
     * bool IS_STRAT        : In stratosphere? 
     * bool IS_PSC          : Polic stratospheric clouds? 
     * double CLD_BrNO3_RC  : Cloud BrNO3 hydrolysis 
     * bool NATSURFACE      : Frozen HNO3? */
    
    bool DO_EDUCT;
    unsigned int N;
    double XSTKCF, ADJUSTEDRATE;
    double HET_BrNO3;

    /* Initialize */
    HET_BrNO3 = 0.0E+00;
    ADJUSTEDRATE = 0.0E+00;
    XSTKCF = 0.0E+00;

    /* Only apply PSC rate adjustment if at high altitude */
    DO_EDUCT = IS_STRAT;

    /* Loop over aerosol types */
    for ( N = 0; N < NAERO; N++ ) {

        if ( N == 2 ) {
            XSTKCF = 8.0E-01;
        } else if ( N == 1 ) {
            XSTKCF = KHETI_SLA[ 5];
        } else if ( N == 0 ) {
            if ( NATSURFACE )
                XSTKCF = 1.00E-03;
            else
                XSTKCF = 3.00E-01;
        } else {
            XSTKCF = 0.0E+00;
        }

        if ( N == 1 ) {
            /* Calculate for stratospheric liquid aerosol */
            ADJUSTEDRATE = AREA[N] * XSTKCF;
        } else {
            /* Reaction rate for surface of aerosol */
            ADJUSTEDRATE = ARSL1K( AREA[N], RADI[N], AIRDENS, XSTKCF, pow( TEMP, 0.5 ), pow( A, 0.5 ) );
        }

        if ( ( DO_EDUCT ) && ( N < 2 ) ) {
            /* PSC reaction - prevent excessive reaction rate */
            if ( ADJUSTEDRATE > 1.0E+00 / PSCMINLIFE )
                ADJUSTEDRATE = 1.0E+00 / PSCMINLIFE;

        }

        /* Add to overall reaction rate */
        HET_BrNO3 += ADJUSTEDRATE;

    }

    if ( !IS_PSC ) 
        HET_BrNO3 += CLD_BrNO3_RC;

    return HET_BrNO3;

} /* End of HETBrNO3 */

double HETHOBr( const double A, const double B, const double KHETI_SLA[11], const double AREA[NAERO], const double RADI[NAERO], const double TEMP, const double AIRDENS, bool IS_STRAT )
{

    /* DESCRIPTION: Set the heterogeneous chemistry rate for HOBr */

    /* INPUTS:
     * double A             : Molar weight in g/mol
     * double B             : 
     * double KHETI_SLA[11] : Sticking coefficients 
     * double AREA[NAERO]   : Aerosol surface area in m^2/cm^3 
     * double RADI[NAERO]   : Aerosol radius in m 
     * double AIRDENS       : Air density in molec/cm^3 
     * bool IS_STRAT        : In stratosphere? */
    
    bool DO_EDUCT;
    unsigned int N;
    double XSTKCF, ADJUSTEDRATE;
    double HET_HOBr;

    /* Initialize */
    HET_HOBr = 0.0E+00;
    ADJUSTEDRATE = 0.0E+00;
    XSTKCF = 0.0E+00;

    /* Only apply PSC rate adjustment if at high altitude */
    DO_EDUCT = IS_STRAT;

    /* Loop over aerosol types */
    for ( N = 0; N < NAERO; N++ ) {

        if ( N == 2 ) {
            XSTKCF = B;
        } else if ( N == 1 ) {
            XSTKCF = KHETI_SLA[10];
        } else if ( N == 0 ) {
            XSTKCF = 1.00E-01;
        } else {
            XSTKCF = 0.0E+00;
        }

        if ( N == 1 ) {
            /* Calculate for stratospheric liquid aerosol */
            ADJUSTEDRATE = AREA[N] * XSTKCF;
        } else {
            /* Reaction rate for surface of aerosol */
            ADJUSTEDRATE = ARSL1K( AREA[N], RADI[N], AIRDENS, XSTKCF, pow( TEMP, 0.5 ), pow( A, 0.5 ) );
        }

        if ( ( DO_EDUCT ) && ( N < 2 ) ) {
            /* PSC reaction - prevent excessive reaction rate */
            if ( ADJUSTEDRATE > 1.0E+00 / PSCMINLIFE )
                ADJUSTEDRATE = 1.0E+00 / PSCMINLIFE;

        }

        /* Add to overall reaction rate */
        HET_HOBr += ADJUSTEDRATE;

    }

    return HET_HOBr;

} /* End of HETHOBr */

double HETHOBr_ice( )
{

    /* DESCRIPTION: Set the heterogeneous chemistry rate for HOBr_ice */

    /* INPUTS:
     * double A             : Molar weight in g/mol
     * double B             : */

    double HET_HOBr_ice;

    /* Initialize */
    HET_HOBr_ice = 0.0E+00;

    // !@#$
    HET_HOBr_ice = 0.0E+00; //KI_HOBr;

    return HET_HOBr_ice;

} /* End of HETHOBr_ice */

double HETHBr_ice( )
{

    /* DESCRIPTION: Set the heterogeneous chemistry rate for HBr_ice */

    /* INPUTS:
     * double A             : Molar weight in g/mol
     * double B             : */

    double HET_HBr_ice;

    /* Initialize */
    HET_HBr_ice = 0.0E+00;

    // !@#$
    HET_HBr_ice = 0.0E+00; //KI_HBr;

    return HET_HBr_ice;

} /* End of HETHBr_ice */

double ARSL1K( const double AREA, const double RADI, const double AIRDENS, const double STKCF, const double XTEMP, const double SQM )
{

    /* DESCRIPTION: Returns the 1st-order loss rate of species on wet aerosol surface */

    /* INPUTS:
     * double AREA    : Aerosol surface area in m^2/cm^3 
     * double RADI    : Aerosol radius in m 
     * double AIRDENS : Air density in molec/cm^3
     * double STKCF   : Sticking coefficient 
     * double XTEMP   : Square root of the temperature in K 
     * double SQM     : Square root of the molecular weight in g/mole */

    /* The 1st-order loss rate on wet aerosol is computed as:
     *
     *  ARSL1K [1/s] = AREA / [ RADI / DFKG + 4.0 / ( STKCF * XMMS ) ]
     *
     *  where XMMS = Mean molecular speed [m/s] = sqrt(8R*T/PI/M)
     *        DFKG = Gas phase diffusion coeff [m^2/s] (order of 0.1) */

    double ARS_L1K;

    double DFKG;

    /* ARSL1K begins here! */

    if ( ( AREA < 0.0E+00 ) || ( AIRDENS < 1.0E-30 ) || ( RADI  < 1.0E-30 ) || \
         ( SQM  < 1.0E-30 ) || ( XTEMP   < 1.0E-30 ) || ( STKCF < 1.0E-30 ) ) {

        /* Use default value if any of the above values are zero. This will prevent
         * division by 0 in the equation below */

        ARS_L1K = 1.0E-30;

    } else {

        /* DFKG = Gas phase diffusion coeff in m^2/s. ~ 0.1 */
        DFKG = 9.45E+13 / AIRDENS * XTEMP * pow( 3.472E-02 + 1.0E+00 / ( SQM * SQM ), 0.5 );

        /* [m^2/cm^3]*[cm^3/m^3] / ( [m]/[m^2/s] + [s/m]) = [m^2/m^3] / ([s/m]) = [1/s] */
        ARS_L1K = AREA * 1.0E+06 / ( RADI / DFKG + 2.749064E-02 * SQM / ( STKCF * XTEMP ));

    }

    return ARS_L1K;

} /* End of ARSL1K */


double N2O5( unsigned int N, const double TEMP, const double RH )
{

    /* DESCRIPTION: Internal function N2O5 computes the gamma sticking 
     * factor for N2O5 hydrolysis */

    /* INPUTS: 
     * unsigned int N : Aerosol type
     * double TEMP    : Temperature in K
     * double RH      : Relative humidity in [0-1] */

    double RH_P, FACT, TTEMP;

    /* N2O5 begins here! */

    /* Convert RH to % ( max = 100% ) */
    RH_P = MIN( RH * 100.0E+00, 100.0E+00 );

    /* Default value */
    double GAMMA_N2O5 = 1.00E-02;

    switch (N) {
        
        /* Sulfate */
        case 2:

            /*************************************/
            /* RH dependence from Kane et al.,
             * Heterogeneous uptake of gaseous 
             * N2O5 by (NH4)2SO4 and H2SO4 
             * aerosols */
            /*************************************/
            /* No RH dependence above 50.05% */
            RH_P = MIN( RH_P, 5.00E+01 );

            GAMMA_N2O5 = 2.79E-04 + RH_P * (  1.30E-04 + \
                                    RH_P * ( -3.43E-06 + \
                                    RH_P * (  7.52E-08 ) ) );

            /******************************************************/
            /* Temperature dependence factor 
             * (Cox et al) is of the form:
             *          10^( log10(G294) - 0.04 * (TTEMP - 294) )
             * FACT = ---------------------------------------------
             *                  10^( log10(G294) )
             *
             * where G294 = 1.0E-02 and TTEMP is MAX( TEMP, 282 ) */
            /******************************************************/
            TTEMP = MAX( TEMP, 282.00E+00 );
            FACT  = pow( 10.0E+00, \
                    -2.00E+00 - 4.00E-02 * ( TTEMP - 294.0E+00 ) ) / \
                    1.00E-02;

            /* Apply temperature dependence */
            GAMMA_N2O5 *= FACT;

            break;

        /* Black carbon */
        case 3:

            GAMMA_N2O5 = 5.00E-03;

            break;

        default:
            printf("Not a suitable aerosol surface for N2O5 hydrolysis\n");
            printf("AEROSOL TYPE = %d\n", N);

            break;

    } 

    return GAMMA_N2O5;

} /* End of N2O5 */

double CLD1K_BrNO3( const double AIRDENS, const double TEMP, const double QL, const double CLDF, bool IS_LAND, bool IS_ICE )
{

    /* DESCRIPTION: Calculates the rate constant for heterogeneous cycling of BrNO3 off of cloud particles */

    /* INPUTS:
     * double AIRDENS : Air density in molec/cm^3 
     * double TEMP    : Temperature in K
     * double QL      : Cloud water mixing ratio
     * double CLDF    : Cloud fraction in box 
     * double IS_LAND : Is land?
     * double IS_ICE  : Is ice? */

    // !@#$ Contrail ice particles?

    /* Cloud droplet radius in continental warm clouds [cm] */
    const double XCLDR_CONT = 6.0E-04;

    /* Cloud droplet radius in marine warm clouds [cm] */
    const double XCLDR_MARI = 1.0E-03;

    const double MW_BrNO3 = 1.42E-03;               /* [kg/mol]  */
    const double ALPHA    = 3.0E-01;                /* Sticking coefficient */

    double CLD1K;
    double RADIUS, AIRVOL, XAIRM3, Vc, AREA;
    double SQM, STK, DFKG;

    CLD1K  = 0.0E+00;
    RADIUS = 0.0E+00;
    AIRVOL = 0.0E+00;
    XAIRM3 = 0.0E+00;
    Vc     = 0.0E+00;
    AREA   = 0.0E+00;
    SQM    = 0.0E+00;
    STK    = 0.0E+00;
    DFKG   = 0.0E+00;

    if ( ( IS_LAND ) || ( !IS_ICE ) ) {
        // !@#$
        if ( ( CLDF > 0 ) && ( TEMP > 258.0 ) ) {

            /* Calculate the surface area of cloud droplets in the grid box, assuming
             * a. marine warm cloud 
             * or 
             * b. continental warm cloud */

            if ( IS_LAND )
                RADIUS = XCLDR_CONT;
            else
                RADIUS = XCLDR_MARI;

            /* Convert the volume of air to cm^3 */
            XAIRM3 = AIRVOL * 1.0E+06;

            /* Get cloud volume */
            Vc     = CLDF * QL * XAIRM3;

            /* Calculate the cloud droplet surface area */
            AREA   = 3.0E+00 * (Vc/XAIRM3) / (RADIUS);

            /* Calculate the 1st-order rate constant for BrNO3 hydrolysis 
             * a. Calculate the gas phase diffusion coefficient
             * b. Calculate the hydrolysis reaction rate*/

            SQM = sqrt( MW_BrNO3 * 1.0E+03 );
            STK = sqrt( TEMP );
            
            /* Gas diffusion coefficient in [cm^2/s] */
            DFKG = 9.45E+17 / AIRDENS * STK * sqrt( 3.472E-02 + 1.0E+00 / ( SQM * SQM ) );

            /* Compute ARSL1K */
            CLD1K = AREA / ( RADIUS/DFKG + 2.749064E-04 * SQM / ( ALPHA * STK ) );

        } else {
            CLD1K = 0.0E+00;
        }
    } else {
        CLD1K = 0.0E+00;
    }

    return CLD1K;

} /* End of CLD1K_BrNO3 */

double CLDICE_HBrHOBr( const double AIRDENS, const double TEMP, const double QI, const double CLDF, const double HBr, const double HOBr, double &K_HBr, double &K_HOBr, const double AREA_ICE, const double EFFRADI_ICE, const double IWC_ICE )
{

    /* DESCRIPTION: Calculates the rate constants for HBr and HOBr pseudo-reactions with ice */

    /* INPUTS:
     * double AIRDENS  : Air density in molec/cm^3 
     * double TEMP     : Temperature in K
     * double IWC      : Ice water content in kg/cm^3
     * double HBr      : HBr concentration in molec/cm^3
     * double HOBr     : HOBr concentration in molec/cm^3
     * double K_HBr    : 
     * double K_HOBr   : 
     * double AREA_ICE : Ice area in m^2/cm^3
     * double RAD_ICE  : Ice radius in m^2/cm^3 */

    const double DENS_ICE = 0.9167E-03; /* [kg/cm^3] */
    const double MW_HBr   = 81.0E-03;   /* [kg/mol] */
    const double MW_HOBr  = 97.0E-03;   /* [kg/mol] */

    double RADI, AREA, IWC;
    double STK, DFKG, SQM_HBr, SQM_HOBr;
    double B_PARAM;
    double GAMMA_ICE;
    // double GAMMA_HBr, GAMMA_HOBr;
    double HBr_RTEMP, HOBr_RTEMP;
    double CLD1K_HBr, CLD1K_HOBr;

    RADI       = 0.00E+00;
    AREA       = 0.00E+00;
    IWC        = 0.00E+00;
    STK        = 0.00E+00;
    DFKG       = 0.00E+00;
    SQM_HBr    = 0.00E+00;
    SQM_HOBr   = 0.00E+00;
    B_PARAM    = 0.00E+00;
    GAMMA_ICE  = 0.00E+00;
    //GAMMA_HBr  = 0.00E+00;
    //GAMMA_HOBr = 0.00E+00;
    HBr_RTEMP  = 0.00E+00;
    HOBr_RTEMP = 0.00E+00;
    CLD1K_HBr  = 0.00E+00;
    CLD1K_HOBr = 0.00E+00;

    if ( IWC_ICE > 0.0E+00 ) {
        IWC = IWC_ICE;
    } else {
        /* Compute our own ICE */
        IWC = CLDF * QI * DENS_ICE;
    }

    if ( IWC <= 0.0E+00 ) {
        GAMMA_ICE  = 0.00E+00;
        K_HBr      = 0.00E+00;
        K_HOBr     = 0.00E+00;
        return GAMMA_ICE;
    }


    /* Calculate the temperature dependent reactive uptake coefficient 
     * for HBr + HOBr + ice */

    if ( ( TEMP >= 180.0 ) && ( TEMP <= 268.0 ) ) {
        GAMMA_ICE = 1.00E-01;
    } else {
        GAMMA_ICE  = 0.00E+00;
        K_HBr      = 0.00E+00;
        K_HOBr     = 0.00E+00;
        return GAMMA_ICE;
    }

    /* Set the sticking coefficients for HBr and HOBr independently */
    // GAMMA_HOBr = 3.00E-03;
    // GAMMA_HBr  = 3.00E-02;

    /* Calculate the surface area of cloud droplets in the given grid box */
    if ( EFFRADI_ICE > 0.00E+00 ) {
        RADI = EFFRADI_ICE;
    } else {
        /* Calculate our own RADI */
        B_PARAM = -2.0E+00 + 1.0E-03 * pow( 273.0 - TEMP, 1.5 ) * log10( IWC/50.0E+00 );

        /* In mum */
        RADI    =  3.774E+02 + B_PARAM * ( 2.033E+02 + \
                               B_PARAM * ( 3.791E+01 + \
                               B_PARAM * ( 2.3696E+00 ) ) );

        /* Convert from mum to m */
        RADI   /= 1.00E+06;

    }

    if ( ( RADI <= 0.0E+00 ) ) {
        GAMMA_ICE = 0.00E+00;
        K_HBr     = 0.00E+00;
        K_HOBr    = 0.00E+00;
        return GAMMA_ICE;
    }

    if ( AREA_ICE > 0.00E+00 ) {
        AREA = AREA_ICE;
    } else {
        /* Calculate our own AREA */

        /* Convert IWC from kg/cm^3 to g/m^3 */
        AREA  = 1.00E-04 * pow( IWC * 1.00E+03 * 1.00E+06, 0.9 );

        /* Convert to m^2/cm^3 */
        AREA *= 2.00E-04;

    }

    SQM_HBr  = sqrt( MW_HBr  * 1.00E+03 );
    SQM_HOBr = sqrt( MW_HOBr * 1.00E+03 );
    STK      = sqrt( TEMP );

    /* Deal with HBr */
    
    /* DFKG = Gas phase diffusion coeff in m^2/s. ~ 0.1 */
    DFKG = 9.45E+13 / AIRDENS * STK * pow( 3.472E-02 + 1.0E+00 / ( SQM_HBr * SQM_HBr ), 0.5 );

    /* [m^2/cm^3]*[cm^3/m^3] / ( [m]/[m^2/s] + [s/m]) = [m^2/m^3] / ([s/m]) = [1/s] */
    CLD1K_HBr = AREA * 1.0E+06 / ( RADI / DFKG + 2.749064E-02 * SQM_HBr / ( GAMMA_ICE * STK ));

    /* Deal with HOBr */

    /* DFKG = Gas phase diffusion coeff in m^2/s. ~ 0.1 */
    DFKG = 9.45E+13 / AIRDENS * STK * pow( 3.472E-02 + 1.0E+00 / ( SQM_HOBr * SQM_HOBr ), 0.5 );

    /* [m^2/cm^3]*[cm^3/m^3] / ( [m]/[m^2/s] + [s/m]) = [m^2/m^3] / ([s/m]) = [1/s] */
    CLD1K_HOBr = AREA * 1.0E+06 / ( RADI / DFKG + 2.749064E-02 * SQM_HOBr / ( GAMMA_ICE * STK ));

    /* Test which loss rate (HOBr or HBr) is limiting */
    HBr_RTEMP  = CLD1K_HBr  * HBr;
    HOBr_RTEMP = CLD1K_HOBr * HOBr;

    if ( HBr_RTEMP > HOBr_RTEMP ) {

        /* Is it safe to divide */
        if ( SafeDiv( CLD1K_HOBr * HOBr, HBr ) ) {
            CLD1K_HBr = HOBr_RTEMP / HBr;
        } else {
            CLD1K_HOBr = 1.00E-30;
            CLD1K_HBr  = 1.00E-30;
        }

    } else {

        /* Is it safe to divide */
        if ( SafeDiv( CLD1K_HBr * HBr, HOBr ) ) {
            CLD1K_HBr = HBr_RTEMP / HOBr;
        } else {
            CLD1K_HOBr = 1.00E-30;
            CLD1K_HBr  = 1.00E-30;
        }

    }

    /* Store rates */
    K_HBr  = CLD1K_HBr;
    K_HOBr = CLD1K_HOBr;

    return GAMMA_ICE;


} /* End of CLDICE_HBrHOBr */

/********************************************************************/
/***                                                              ***/
/*** The following functions are defined for UCX-based mechanisms ***/
/***                                                              ***/
/********************************************************************/

double HETN2O5_PSC( const double A, const double B, const double KHETI_SLA[11], const double AREA[NAERO], const double RADI[NAERO], const double TEMP, const double AIRDENS, bool IS_STRAT, bool NATSURFACE )
{

    /* DESCRIPTION: Set the heterogeneous chemistry rate for N2O5(g) + HCl(l,s) in PSCs */

    /* INPUTS:
     * double A             : Molar weight in g/mol
     * double B             : 
     * double KHETI_SLA[11] : Sticking coefficients 
     * double AREA[NAERO]   : Aerosol surface area in m^2/cm^3 
     * double RADI[NAERO]   : Aerosol radius in m 
     * double AIRDENS       : Air density in molec/cm^3 
     * bool IS_STRAT        : In the stratosphere?
     * bool NATSURFACE      : Frozen HNO3? */

    bool DO_EDUCT;
    unsigned int N;
    double XSTKCF, ADJUSTEDRATE;
    double HET_N2O5_PSC;

    /* Initialize */
    HET_N2O5_PSC = 0.0E+00;
    ADJUSTEDRATE = 0.0E+00;
    XSTKCF = 0.0E+00;

    /* Always apply PSC rate adjustment*/
    DO_EDUCT = 1;
    
    /* Loop over aerosol types */
    for ( N = 0; N < NAERO; N++ ) {

        /* Only consider PSC reactions in strat */
        if ( IS_STRAT ) {
            if ( N == 2 ) { 
                XSTKCF = 1.0E-05;
            } else if ( N == 1 ) {
                XSTKCF = KHETI_SLA[ 1];
            } else if ( N == 0 ) {
                if ( NATSURFACE )
                    XSTKCF = 3.00E-03;
                else
                    XSTKCF = 3.00E-02;
            } else {
                XSTKCF = 0.00E+00;
            }
        } else {
            XSTKCF = B;
        }

        if ( N == 1 ) {
            /* Calculate for stratospheric liquid aerosol */
            ADJUSTEDRATE = AREA[N] * XSTKCF;
        } else {
            /* Reaction rate for surface of aerosol */
            ADJUSTEDRATE = ARSL1K( AREA[N], RADI[N], AIRDENS, XSTKCF, pow( TEMP, 0.5 ), pow( A, 0.5 ) );
        }

        if ( ( DO_EDUCT ) && ( N < 2 ) ) {
            /* PSC reaction - prevent excessive reaction rate */
            if ( ADJUSTEDRATE > 1.0E+00 / PSCMINLIFE )
                ADJUSTEDRATE = 1.0E+00 / PSCMINLIFE;

        }

        /* Add to overall reaction rate */
        HET_N2O5_PSC += ADJUSTEDRATE;

    }

    return HET_N2O5_PSC;

} /* End of HETN2O5_PSC */

double HETClNO3_PSC1( const double A, const double B, const double KHETI_SLA[11], const double AREA[NAERO], const double RADI[NAERO], const double TEMP, const double AIRDENS, bool IS_STRAT, bool NATSURFACE )
{

    /* DESCRIPTION: Set the heterogeneous chemistry rate for ClNO3(g) + H2O(l,s) in PSCs */

    /* INPUTS:
     * double A             : Molar weight in g/mol
     * double B             : 
     * double KHETI_SLA[11] : Sticking coefficients 
     * double AREA[NAERO]   : Aerosol surface area in m^2/cm^3 
     * double RADI[NAERO]   : Aerosol radius in m 
     * double AIRDENS       : Air density in molec/cm^3 
     * bool IS_STRAT        : In the stratosphere?
     * bool NATSURFACE      : Frozen HNO3? */

    bool DO_EDUCT;
    unsigned int N;
    double XSTKCF, ADJUSTEDRATE;
    double HET_ClNO3_PSC1;

    /* Initialize */
    HET_ClNO3_PSC1 = 0.0E+00;
    ADJUSTEDRATE = 0.0E+00;
    XSTKCF = 0.0E+00;

    /* Always apply PSC rate adjustment*/
    DO_EDUCT = 1;
    
    /* Loop over aerosol types */
    for ( N = 0; N < NAERO; N++ ) {

        /* Only consider PSC reactions in strat */
        if ( IS_STRAT ) {
            if ( N == 2 ) { 
                XSTKCF = 1.0E-04;
            } else if ( N == 1 ) {
                XSTKCF = KHETI_SLA[ 2];
            } else if ( N == 0 ) {
                if ( NATSURFACE )
                    XSTKCF = 4.00E-03;
                else
                    XSTKCF = 3.00E-01;
            } else {
                XSTKCF = 0.00E+00;
            }
        } else {
            XSTKCF = B;
        }

        if ( N == 1 ) {
            /* Calculate for stratospheric liquid aerosol */
            ADJUSTEDRATE = AREA[N] * XSTKCF;
        } else {
            /* Reaction rate for surface of aerosol */
            ADJUSTEDRATE = ARSL1K( AREA[N], RADI[N], AIRDENS, XSTKCF, pow( TEMP, 0.5 ), pow( A, 0.5 ) );
        }

        if ( ( DO_EDUCT ) && ( N < 2 ) ) {
            /* PSC reaction - prevent excessive reaction rate */
            if ( ADJUSTEDRATE > 1.0E+00 / PSCMINLIFE )
                ADJUSTEDRATE = 1.0E+00 / PSCMINLIFE;

        }

        /* Add to overall reaction rate */
        HET_ClNO3_PSC1 += ADJUSTEDRATE;

    }
    
    return HET_ClNO3_PSC1;

} /* End of HETClNO3_PSC1 */

double HETClNO3_PSC2( const double A, const double B, const double KHETI_SLA[11], const double AREA[NAERO], const double RADI[NAERO], const double TEMP, const double AIRDENS, bool IS_STRAT, bool NATSURFACE )
{

    /* DESCRIPTION: Set the heterogeneous chemistry rate for ClNO3(g) + HCl(l,s) in PSCs */

    /* INPUTS:
     * double A             : Molar weight in g/mol
     * double B             : 
     * double KHETI_SLA[11] : Sticking coefficients 
     * double AREA[NAERO]   : Aerosol surface area in m^2/cm^3 
     * double RADI[NAERO]   : Aerosol radius in m 
     * double AIRDENS       : Air density in molec/cm^3 
     * bool IS_STRAT        : In the stratosphere?
     * bool NATSURFACE      : Frozen HNO3? */

    bool DO_EDUCT;
    unsigned int N;
    double XSTKCF, ADJUSTEDRATE;
    double HET_ClNO3_PSC2;

    /* Initialize */
    HET_ClNO3_PSC2 = 0.0E+00;
    ADJUSTEDRATE = 0.0E+00;
    XSTKCF = 0.0E+00;

    /* Always apply PSC rate adjustment*/
    DO_EDUCT = 1;
    
    /* Loop over aerosol types */
    for ( N = 0; N < NAERO; N++ ) {

        /* Only consider PSC reactions in strat */
        if ( IS_STRAT ) {
            if ( N == 2 ) { 
                XSTKCF = 1.0E-05;
            } else if ( N == 1 ) {
                XSTKCF = KHETI_SLA[ 3];
            } else if ( N == 0 ) {
                if ( NATSURFACE )
                    XSTKCF = 2.00E-01;
                else
                    XSTKCF = 3.00E-01;
            } else {
                XSTKCF = 0.00E+00;
            }
        } else {
            XSTKCF = B;
        }

        if ( N == 1 ) {
            /* Calculate for stratospheric liquid aerosol */
            ADJUSTEDRATE = AREA[N] * XSTKCF;
        } else {
            /* Reaction rate for surface of aerosol */
            ADJUSTEDRATE = ARSL1K( AREA[N], RADI[N], AIRDENS, XSTKCF, pow( TEMP, 0.5 ), pow( A, 0.5 ) );
        }

        if ( ( DO_EDUCT ) && ( N < 2 ) ) {
            /* PSC reaction - prevent excessive reaction rate */
            if ( ADJUSTEDRATE > 1.0E+00 / PSCMINLIFE )
                ADJUSTEDRATE = 1.0E+00 / PSCMINLIFE;

        }

        /* Add to overall reaction rate */
        HET_ClNO3_PSC2 += ADJUSTEDRATE;

    }

    return HET_ClNO3_PSC2;
    
} /* End of HETClNO3_PSC2 */

double HETClNO3_PSC3( const double A, const double B, const double KHETI_SLA[11], const double AREA[NAERO], const double RADI[NAERO], const double TEMP, const double AIRDENS, bool IS_STRAT, bool NATSURFACE )
{

    /* DESCRIPTION: Set the heterogeneous chemistry rate for ClNO3(g) + HBr(l,s) in PSCs */

    /* INPUTS:
     * double A             : Molar weight in g/mol
     * double B             : 
     * double KHETI_SLA[11] : Sticking coefficients 
     * double AREA[NAERO]   : Aerosol surface area in m^2/cm^3 
     * double RADI[NAERO]   : Aerosol radius in m 
     * double AIRDENS       : Air density in molec/cm^3 
     * bool IS_STRAT        : In the stratosphere?
     * bool NATSURFACE      : Frozen HNO3? */

    bool DO_EDUCT;
    unsigned int N;
    double XSTKCF, ADJUSTEDRATE;
    double HET_ClNO3_PSC3;

    /* Initialize */
    HET_ClNO3_PSC3 = 0.0E+00;
    ADJUSTEDRATE = 0.0E+00;
    XSTKCF = 0.0E+00;

    /* Always apply PSC rate adjustment*/
    DO_EDUCT = 1;
    
    /* Loop over aerosol types */
    for ( N = 0; N < NAERO; N++ ) {

        /* Only consider PSC reactions in strat */
        if ( IS_STRAT ) {
            if ( N == 2 ) { 
                XSTKCF = 0.00E+00;
            } else if ( N == 1 ) {
                XSTKCF = KHETI_SLA[ 4];
            } else if ( N == 0 ) {
                if ( NATSURFACE )
                    XSTKCF = 3.00E-01;
                else
                    XSTKCF = 3.00E-01;
            } else {
                XSTKCF = 0.00E+00;
            }
        } else {
            XSTKCF = B;
        }

        if ( N == 1 ) {
            /* Calculate for stratospheric liquid aerosol */
            ADJUSTEDRATE = AREA[N] * XSTKCF;
        } else {
            /* Reaction rate for surface of aerosol */
            ADJUSTEDRATE = ARSL1K( AREA[N], RADI[N], AIRDENS, XSTKCF, pow( TEMP, 0.5 ), pow( A, 0.5 ) );
        }

        if ( ( DO_EDUCT ) && ( N < 2 ) ) {
            /* PSC reaction - prevent excessive reaction rate */
            if ( ADJUSTEDRATE > 1.0E+00 / PSCMINLIFE )
                ADJUSTEDRATE = 1.0E+00 / PSCMINLIFE;

        }

        /* Add to overall reaction rate */
        HET_ClNO3_PSC3 += ADJUSTEDRATE;

    }

    return HET_ClNO3_PSC3;

} /* End of HETClNO3_PSC3 */

double HETBrNO3_PSC( const double A, const double B, const double KHETI_SLA[11], const double AREA[NAERO], const double RADI[NAERO], const double TEMP, const double AIRDENS, bool IS_STRAT, bool NATSURFACE )
{

    /* DESCRIPTION: Set the heterogeneous chemistry rate for BrNO3(g) + HCl(l,s) in PSCs */

    /* INPUTS:
     * double A             : Molar weight in g/mol
     * double B             : 
     * double KHETI_SLA[11] : Sticking coefficients 
     * double AREA[NAERO]   : Aerosol surface area in m^2/cm^3 
     * double RADI[NAERO]   : Aerosol radius in m 
     * double AIRDENS       : Air density in molec/cm^3 
     * bool IS_STRAT        : In the stratosphere?
     * bool NATSURFACE      : Frozen HNO3? */

    bool DO_EDUCT;
    unsigned int N;
    double XSTKCF, ADJUSTEDRATE;
    double HET_BrNO3_PSC;

    /* Initialize */
    HET_BrNO3_PSC = 0.0E+00;
    ADJUSTEDRATE = 0.0E+00;
    XSTKCF = 0.0E+00;

    /* Always apply PSC rate adjustment*/
    DO_EDUCT = 1;
    
    /* Loop over aerosol types */
    for ( N = 0; N < NAERO; N++ ) {

        /* Only consider PSC reactions in strat */
        if ( IS_STRAT ) {
            if ( N == 2 ) { 
                XSTKCF = 9.0E-01;
            } else if ( N == 1 ) {
                XSTKCF = KHETI_SLA[ 6];
            } else if ( N == 0 ) {
                if ( NATSURFACE )
                    XSTKCF = 3.00E-01;
                else
                    XSTKCF = 3.00E-01;
            } else {
                XSTKCF = 0.00E+00;
            }
        } else {
            XSTKCF = B;
        }

        if ( N == 1 ) {
            /* Calculate for stratospheric liquid aerosol */
            ADJUSTEDRATE = AREA[N] * XSTKCF;
        } else {
            /* Reaction rate for surface of aerosol */
            ADJUSTEDRATE = ARSL1K( AREA[N], RADI[N], AIRDENS, XSTKCF, pow( TEMP, 0.5 ), pow( A, 0.5 ) );
        }

        if ( ( DO_EDUCT ) && ( N < 2 ) ) {
            /* PSC reaction - prevent excessive reaction rate */
            if ( ADJUSTEDRATE > 1.0E+00 / PSCMINLIFE )
                ADJUSTEDRATE = 1.0E+00 / PSCMINLIFE;

        }

        /* Add to overall reaction rate */
        HET_BrNO3_PSC += ADJUSTEDRATE;

    }
    
    return HET_BrNO3_PSC;

} /* End of HETBrNO3_PSC */

double HETHOCl_PSC1( const double A, const double B, const double KHETI_SLA[11], const double AREA[NAERO], const double RADI[NAERO], const double TEMP, const double AIRDENS, bool IS_STRAT, bool NATSURFACE )
{

    /* DESCRIPTION: Set the heterogeneous chemistry rate for HOCl(g) + HCl(l,s) in PSCs */

    /* INPUTS:
     * double A             : Molar weight in g/mol
     * double B             :
     * double KHETI_SLA[11] : Sticking coefficients
     * double AREA[NAERO]   : Aerosol surface area in m^2/cm^3
     * double RADI[NAERO]   : Aerosol radius in m
     * double AIRDENS       : Air density in molec/cm^3
     * bool IS_STRAT        : In the stratosphere?
     * bool NATSURFACE      : Frozen HNO3? */

    bool DO_EDUCT;
    unsigned int N;
    double XSTKCF, ADJUSTEDRATE;
    double HET_HOCl_PSC1;

    /* Initialize */
    HET_HOCl_PSC1 = 0.0E+00;
    ADJUSTEDRATE = 0.0E+00;
    XSTKCF = 0.0E+00;

    /* Always apply PSC rate adjustment*/
    DO_EDUCT = 1;
    
    /* Loop over aerosol types */
    for ( N = 0; N < NAERO; N++ ) {

        /* Only consider PSC reactions in strat */
        if ( IS_STRAT ) {
            if ( N == 2 ) { 
                XSTKCF = 8.0E-01;
            } else if ( N == 1 ) {
                XSTKCF = KHETI_SLA[ 7];
            } else if ( N == 0 ) {
                if ( NATSURFACE )
                    XSTKCF = 1.00E-01;
                else
                    XSTKCF = 2.00E-01;
            } else {
                XSTKCF = 0.00E+00;
            }
        } else {
            XSTKCF = B;
        }

        if ( N == 1 ) {
            /* Calculate for stratospheric liquid aerosol */
            ADJUSTEDRATE = AREA[N] * XSTKCF;
        } else {
            /* Reaction rate for surface of aerosol */
            ADJUSTEDRATE = ARSL1K( AREA[N], RADI[N], AIRDENS, XSTKCF, pow( TEMP, 0.5 ), pow( A, 0.5 ) );
        }

        if ( ( DO_EDUCT ) && ( N < 2 ) ) {
            /* PSC reaction - prevent excessive reaction rate */
            if ( ADJUSTEDRATE > 1.0E+00 / PSCMINLIFE )
                ADJUSTEDRATE = 1.0E+00 / PSCMINLIFE;

        }

        /* Add to overall reaction rate */
        HET_HOCl_PSC1 += ADJUSTEDRATE;

    }
    
    return HET_HOCl_PSC1;

} /* End of HETHOCl_PSC1 */

double HETHOCl_PSC2( const double A, const double B, const double KHETI_SLA[11], const double AREA[NAERO], const double RADI[NAERO], const double TEMP, const double AIRDENS, bool IS_STRAT, bool NATSURFACE )
{

    /* DESCRIPTION: Set the heterogeneous chemistry rate for HOCl(g) + HBr(l,s) in PSCs */

    /* INPUTS:
     * double A             : Molar weight in g/mol
     * double B             : 
     * double KHETI_SLA[11] : Sticking coefficients 
     * double AREA[NAERO]   : Aerosol surface area in m^2/cm^3 
     * double RADI[NAERO]   : Aerosol radius in m 
     * double AIRDENS       : Air density in molec/cm^3 
     * bool IS_STRAT        : In the stratosphere?
     * bool NATSURFACE      : Frozen HNO3? */

    bool DO_EDUCT;
    unsigned int N;
    double XSTKCF, ADJUSTEDRATE;
    double HET_HOCl_PSC2;

    /* Initialize */
    HET_HOCl_PSC2 = 0.0E+00;
    ADJUSTEDRATE = 0.0E+00;
    XSTKCF = 0.0E+00;

    /* Always apply PSC rate adjustment*/
    DO_EDUCT = 1;
    
    /* Loop over aerosol types */
    for ( N = 0; N < NAERO; N++ ) {

        /* Only consider PSC reactions in strat */
        if ( IS_STRAT ) {
            if ( N == 2 ) { 
                XSTKCF = 8.0E-01;
            } else if ( N == 1 ) {
                XSTKCF = KHETI_SLA[ 8];
            } else if ( N == 0 ) {
                if ( NATSURFACE )
                    XSTKCF = 3.00E-01;
                else
                    XSTKCF = 3.00E-01;
            } else {
                XSTKCF = 0.00E+00;
            }
        } else {
            XSTKCF = B;
        }

        if ( N == 1 ) {
            /* Calculate for stratospheric liquid aerosol */
            ADJUSTEDRATE = AREA[N] * XSTKCF;
        } else {
            /* Reaction rate for surface of aerosol */
            ADJUSTEDRATE = ARSL1K( AREA[N], RADI[N], AIRDENS, XSTKCF, pow( TEMP, 0.5 ), pow( A, 0.5 ) );
        }

        if ( ( DO_EDUCT ) && ( N < 2 ) ) {
            /* PSC reaction - prevent excessive reaction rate */
            if ( ADJUSTEDRATE > 1.0E+00 / PSCMINLIFE )
                ADJUSTEDRATE = 1.0E+00 / PSCMINLIFE;

        }

        /* Add to overall reaction rate */
        HET_HOCl_PSC2 += ADJUSTEDRATE;

    }

    return HET_HOCl_PSC2;

} /* End of HETHOCl_PSC2 */

double HETHOBr_PSC( const double A, const double B, const double KHETI_SLA[11], const double AREA[NAERO], const double RADI[NAERO], const double TEMP, const double AIRDENS, bool IS_STRAT, bool NATSURFACE )
{

    /* DESCRIPTION: Set the heterogeneous chemistry rate for HOBr(g) + HCl(l,s) in PSCs */

    /* INPUTS:
     * double A             : Molar weight in g/mol
     * double B             : 
     * double KHETI_SLA[11] : Sticking coefficients 
     * double AREA[NAERO]   : Aerosol surface area in m^2/cm^3 
     * double RADI[NAERO]   : Aerosol radius in m 
     * double AIRDENS       : Air density in molec/cm^3 
     * bool IS_STRAT        : In the stratosphere?
     * bool NATSURFACE      : Frozen HNO3? */

    bool DO_EDUCT;
    unsigned int N;
    double XSTKCF, ADJUSTEDRATE;
    double HET_HOBr_PSC;

    /* Initialize */
    HET_HOBr_PSC = 0.0E+00;
    ADJUSTEDRATE = 0.0E+00;
    XSTKCF = 0.0E+00;

    /* Always apply PSC rate adjustment*/
    DO_EDUCT = 1;
    
    /* Loop over aerosol types */
    for ( N = 0; N < NAERO; N++ ) {

        /* Only consider PSC reactions in strat */
        if ( IS_STRAT ) {
            if ( N == 2 ) { 
                XSTKCF = 8.0E-01;
            } else if ( N == 1 ) {
                XSTKCF = KHETI_SLA[ 9];
            } else if ( N == 0 ) {
                if ( NATSURFACE )
                    XSTKCF = 1.00E-01;
                else
                    XSTKCF = 3.00E-01;
            } else {
                XSTKCF = 0.00E+00;
            }
        } else {
            XSTKCF = B;
        }

        if ( N == 1 ) {
            /* Calculate for stratospheric liquid aerosol */
            ADJUSTEDRATE = AREA[N] * XSTKCF;
        } else {
            /* Reaction rate for surface of aerosol */
            ADJUSTEDRATE = ARSL1K( AREA[N], RADI[N], AIRDENS, XSTKCF, pow( TEMP, 0.5 ), pow( A, 0.5 ) );
        }

        if ( ( DO_EDUCT ) && ( N < 2 ) ) {
            /* PSC reaction - prevent excessive reaction rate */
            if ( ADJUSTEDRATE > 1.0E+00 / PSCMINLIFE )
                ADJUSTEDRATE = 1.0E+00 / PSCMINLIFE;

        }

        /* Add to overall reaction rate */
        HET_HOBr_PSC += ADJUSTEDRATE;

    }
    
    return HET_HOBr_PSC;

} /* End of HETHOBr_PSC */

