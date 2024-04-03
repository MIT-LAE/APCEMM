/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* LiquidAer Program File                                           */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 10/22/2018                                */
/* File                 : LiquidAer.cpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Core/LiquidAer.hpp"

double H2SO4_GASFRAC( const double temperature_K, const double SO4 )
{

    /* DESCRIPTION: Calculates the fraction of SO4 aerosol which can
     * can be considered to be gaseous H2SO4 */

    /* INPUTS:
     * - temperature_K : temperature in K
     * - SO4           : total sulfate (l + g) concentration in molec/cm^3 */

    if ( SO4 < 0.0E+00 ) {
        std::cout << "In LiquidAer.cpp: SO4 is negative: SO4 = " << SO4 << "\n";
        return 0;
    }

    double frac = 0.0E+00;

    const double nSat = physFunc::pSat_H2SO4( temperature_K ) / ( physConst::kB * temperature_K * 1.0E+06 );
    const bool   isSS = ( SO4 > nSat );

    if ( isSS )
        frac = nSat/SO4;
    else
        frac = 1.0;

    return frac;

} /* End of H2SO4_GASFRAC */

unsigned int STRAT_AER( const double temperature_K     , const double pressure_Pa       , const double airDens ,               \
                        const double latitude_deg      , std::vector<double> &Data      , double Area ,                        \
                        std::vector<double> &KHETI_SLA , std::vector<double> &SOLIDFRAC ,                                      \
                        std::vector<double> &AERFRAC   , std::vector<double> &RAD_AER   ,                                      \
                        std::vector<double> &RHO_AER   , std::vector<double> &KG_AER    ,                                      \
                        std::vector<double> &NDENS_AER , std::vector<double> &SAD_AER   , double tropopausePressure, const bool DBG_IN )
{

    /* DESCRIPTION: Calculates areosol properties in the stratosphere using the thermodynamic
     * parameterization described in Kirner et al., 2011. */

    /* INPUTS:
     * - temperature_K : Temperature in K
     * - pressure_Pa   : Pressure in Pa 
     * - airDens       : Air density in molec/cm^3
     * - latitude_deg  : Latitude in degrees
     * - Data          : Ambient conditions
     * - Area          : Grid cell area of the ambient cell in m^2
     *
     * OUTPUTS:
     * - KHETI_SLA
     * - AERFRAC
     * - RAD_AER
     * - RHO_AER
     * - KG_AER
     * - NDENS_AER
     * - SAD_AER
     * - STATE_PSC */
  
    /* Debug? */
    const bool DBG = ( 0 || DBG_IN ) ;

    /* Species indices */
    const unsigned int id_SO4   = 0;
    const unsigned int id_HNO3  = 1;
    const unsigned int id_HCl   = 2;
    const unsigned int id_HOCl  = 3;
    const unsigned int id_HBr   = 4;
    const unsigned int id_HOBr  = 5;
    const unsigned int id_H2O   = 6;
    const unsigned int id_ClNO3 = 7;
    const unsigned int id_BrNO3 = 8;
    const unsigned int id_NIT   = 9;
    const unsigned int id_NAT   = 10;

    /* Limit on PSC formation */
    const double PSC_MAXLAT =  45.0E+00;
    const double PSC_MINLAT = -55.0E+00;
    const double PSC_PMAX   =  18.0E+03;
    const double PSC_PMIN   =   5.0E+02;

    /* Peak pressure at which NAT can form homogeneously */
    const double P_MAXNAT   =  1.40E+04; /* [Pa] */

    /* Maximum temperature for PSC formation */
    const double T_MAX = 215.0E+00; /* [K] */

    /* Limits on NAT/ice formation */
    const double MIN_RAD = 1.0E-07; /* [m] */
    const double MAX_NDENS = 42.0E+03; /* [#/m^3] */

    unsigned int STATE_LOCAL; 

    /* Phase = 0 -> no formation condition for PSC,
     * Phase = 1 -> formation conditions for STS,
     * Phase = 2 -> formation conditions for STS and NAT,
     * Phase = 3 -> formation conditions for STS, NAT and ice. */

    bool IS_POLAR, IS_STRAT, IS_VALID;

    /* Mass of dry air per grid box [m/mol] */
    const double INVAIR = MW_Air / ( pressure_Pa / ( physConst::R_Air * temperature_K ) * Area ); 
    /* Unit check :
     * kg/mol / ( kg/m^3 * m^2 ) = m / mol*/
   
    /* Grid box mixing ratios and partial pressures */
    double HNO3PP, H2OPP;
    double PSATHNO3, PSATH2O, PSATH2O_SUPERSAT, PSATHNO3_SUPERCOOL;
    const double P_ICE_SUPERSAT = 1.2;

    /* Gridbox aerosol and phase data */
    double HNO3SUM     , H2OSUM       , H2SO4SUM;
    double ClNO3SUM    , HOClSUM      , HClSUM;
    double BrNO3SUM    , HOBrSUM      , HBrSUM;
    double HNO3_BOX_S  , HNO3_BOX_G   , HNO3_BOX_L;
    double H2O_BOX_S   , H2O_BOX_G    , H2O_BOX_L;
    double H2SO4_BOX_G , H2SO4_BOX_L;
    double HCl_BOX_G   , HCl_BOX_L;
    double HOCl_BOX_G  , HOCl_BOX_L;
    double HBr_BOX_G   , HBr_BOX_L;
    double HOBr_BOX_G  , HOBr_BOX_L;
    double HNO3_GASFRAC, HCl_GASFRAC, HOCl_GASFRAC, \
           HBr_GASFRAC , HOBr_GASFRAC;
    double VOL_NAT, VOL_ICE, VOL_SLA, VOL_TOT;
    double RAD_AER_BOX, RHO_AER_BOX, KG_AER_BOX, \
           NDENS_AER_BOX, SAD_AER_BOX;
    double KG_NAT, KG_ICE, KG_NO3;
    double H2SO4_gFRAC;

    /* SLA weight fractions */
    double W_H2SO4, W_H2O,  W_HNO3, \
           W_HCl  , W_HOCl, W_HBr , W_HOBr;

    /* Calculate conversion factors for SLA */
    /* Factor to convert volume (m^3/m^3 air) to
     * surface area density (cm^2/cm^3 air) */
    const double SLA_VA = 8.406E-08 * pow( 10.0, 12.0 * 0.751E+00 );

    /* Factor to convert effective radius to liquid radius (unitless) */
    const double SLA_RR = exp( -0.173E+00 );

    /* Factor to convert volume (m^3/m^3) to effective radius (m) */
    const double SLA_VR = 0.357E-06 * pow( 10.0, 12.0 * 0.249E+00 );

    /* Reaction prefactors */
    double KHET_COMMON;
    double KHET_SPECIFIC;
   
    std::vector<double> GAMMA_BOX( 11, 0.0 );

    /* Initialize */
    STATE_LOCAL = 0;

    IS_POLAR = 0;
    IS_STRAT = 0;
    IS_VALID = 0;

    HNO3PP             = 0.0E+00;
    H2OPP              = 0.0E+00;
    PSATHNO3           = 0.0E+00;
    PSATH2O            = 0.0E+00;
    PSATH2O_SUPERSAT   = 0.0E+00;
    PSATHNO3_SUPERCOOL = 0.0E+00;
    HNO3SUM            = 0.0E+00;
    H2OSUM             = 0.0E+00;
    H2SO4SUM           = 0.0E+00;
    ClNO3SUM           = 0.0E+00;
    HOClSUM            = 0.0E+00;
    HClSUM             = 0.0E+00;
    BrNO3SUM           = 0.0E+00;
    HOBrSUM            = 0.0E+00;
    HBrSUM             = 0.0E+00;
    HNO3_BOX_S         = 0.0E+00;
    HNO3_BOX_G         = 0.0E+00;
    HNO3_BOX_L         = 0.0E+00;
    H2O_BOX_S          = 0.0E+00;
    H2O_BOX_G          = 0.0E+00;
    H2O_BOX_L          = 0.0E+00;
    H2SO4_BOX_G        = 0.0E+00;
    H2SO4_BOX_L        = 0.0E+00;
    HCl_BOX_G          = 0.0E+00;
    HCl_BOX_L          = 0.0E+00;
    HOCl_BOX_G         = 0.0E+00;
    HOCl_BOX_L         = 0.0E+00;
    HBr_BOX_G          = 0.0E+00;
    HBr_BOX_L          = 0.0E+00;
    HOBr_BOX_G         = 0.0E+00;
    HOBr_BOX_L         = 0.0E+00;
    HNO3_GASFRAC       = 0.0E+00;
    HCl_GASFRAC        = 0.0E+00;
    HOCl_GASFRAC       = 0.0E+00;
    HBr_GASFRAC        = 0.0E+00;
    HOBr_GASFRAC       = 0.0E+00;
    VOL_NAT            = 0.0E+00;
    VOL_ICE            = 0.0E+00;
    VOL_SLA            = 0.0E+00;
    VOL_TOT            = 0.0E+00;
    RAD_AER_BOX        = 0.0E+00;
    RHO_AER_BOX        = 0.0E+00;
    KG_AER_BOX         = 0.0E+00;
    NDENS_AER_BOX      = 0.0E+00;
    SAD_AER_BOX        = 0.0E+00;
    KG_NAT             = 0.0E+00;
    KG_ICE             = 0.0E+00;
    KG_NO3             = 0.0E+00;
    H2SO4_gFRAC        = 0.0E+00;
    W_H2SO4            = 0.0E+00;
    W_H2O              = 0.0E+00;
    W_HNO3             = 0.0E+00;
    W_HCl              = 0.0E+00;
    W_HOCl             = 0.0E+00;
    W_HBr              = 0.0E+00;
    W_HOBr             = 0.0E+00;
    KHET_COMMON        = 0.0E+00;
    KHET_SPECIFIC      = 0.0E+00;
   

    /* Code begins here ! */
    IS_POLAR = ( ( latitude_deg <= PSC_MINLAT ) || ( latitude_deg >= PSC_MAXLAT ) );

    IS_STRAT = ( pressure_Pa < 300.0E+02 ) ; //TROPP );

    IS_VALID = ( ( IS_POLAR ) && ( IS_STRAT ) && ( pressure_Pa > PSC_PMIN ) && ( pressure_Pa < PSC_PMAX) );
    IS_VALID = ( IS_VALID || PSC_FULL );

    if ( DBG ) {
        std::cout << "IS_POLAR : " << IS_POLAR << "\n";
        std::cout << "IS_STRAT : " << IS_STRAT << "\n";
        std::cout << "IS_VALID : " << IS_VALID << "\n";
    }

    /* Available NO3 mass */
    KG_NO3 = Data[id_HNO3] * 1.0E+06 / physConst::Na * Area * MW_NIT + Data[id_NIT];
    /* Unit check:
     * molec/cm^3 * cm^3/m^3 * mole/molec * m^2 * kg/mol = kg/m */

    /* HNO3 mixing ratio [mol/mol] */
    HNO3SUM  = Data[id_HNO3] / airDens;
    HNO3SUM += Data[id_NIT] / MW_NIT / (airDens * 1.0E+06 / physConst::Na);
    /* Unit check: 
     * kg/m^3 / (kg/mol) / (molec/cm^3 * cm^3/m^3 * mole/molec) = mol NIT/m^3 / mol air/m^3  */

    /* H2O mixing ratio [mol/mol] */
    H2OSUM = Data[id_H2O] / airDens;

    /* Calculate partial pressures */
    HNO3PP = pressure_Pa * HNO3SUM;
    H2OPP  = pressure_Pa * H2OSUM;

    if ( !IS_VALID ) {
        
        /* No PSCs (SSA only) */
        STATE_LOCAL = 0;

    } else if ( temperature_K > T_MAX ) {

        /* STS/SSA only */
        STATE_LOCAL = 1;

    } else {

        /* Calculate saturation pressures */
        PSATHNO3 = physFunc::pSat_HNO3( temperature_K, H2OPP );
        PSATH2O  = physFunc::pSat_H2Os( temperature_K );

        /* Supersaturation requirement for ice */
        PSATH2O_SUPERSAT = PSATH2O * P_ICE_SUPERSAT;

        /* If homogeneous NAT nucleation allowed, calculate 
         * threshold saturation pressure */
        if ( LHOMNUCNAT ) {
            /* Calculate as if temperature is T + T_NAT_SUPERCOOL */
            PSATHNO3_SUPERCOOL = physFunc::pSat_HNO3( temperature_K + \
                                                      T_NAT_SUPERCOOL, H2OPP );

        } else {
            /* Make homogeneous nucleation impossible */
            PSATHNO3_SUPERCOOL = HNO3PP + 1.0E+00;
        }

        /* Only interested in sign */
        if ( STATE_LOCAL > 1 )
            H2O_BOX_S = H2OPP - PSATH2O;
        else
            H2O_BOX_S = H2OPP - PSATH2O_SUPERSAT;

        /* Ice exists/possible? */
        if ( H2O_BOX_S > 0.0E+00 ) {
            STATE_LOCAL = 3;
            HNO3_BOX_S = HNO3PP - PSATHNO3;
        } else {
            /* If ice not possible could still have NAT */
            H2O_BOX_S = 0.0E+00;
            /* 1. Homogeneous nucleation */
            if ( ( LHOMNUCNAT ) && ( STATE_LOCAL == 1 ) )
                HNO3_BOX_S = HNO3PP - PSATHNO3_SUPERCOOL;

            /* 2. Box formerly contained ice or NAT */
            if ( STATE_LOCAL == 2 )
                HNO3_BOX_S = HNO3PP - PSATHNO3;

            if ( HNO3_BOX_S > 0.0E+00 )
                STATE_LOCAL = 2;
            else
                STATE_LOCAL = 1;
        }

    }

    /* Only continue if we want online solid PSCs */
    if ( LSOLIDPSC ) {
        if ( STATE_LOCAL == 3 ) {
            /* Form ice PSCs */
            H2O_BOX_S = ( H2OPP - PSATH2O ) / pressure_Pa;
            H2O_BOX_S = std::max( H2O_BOX_S, 0.0E+00 );
            KG_ICE    = H2O_BOX_S * MW_H2O / physConst::Na * airDens * ( 1.0E+06 ) * Area;
            /* Unit check:
             *        = mole/mole * kg/mol * mole/molec * molec/cm^3 * cm^3/m^3 * m^2
             *        = kg / m */
            VOL_ICE = H2O_BOX_S * airDens * ( 1.0E+06 ) * MW_H2O / ( physConst::RHO_ICE * physConst::Na );
            /* Unit check:
             *      = mole/mole * molec/cm^3 * cm^3/m^3 * kg/mole/ ( kg/m^3 ice         * molec/mole   )
             *      = m^3 ice / m^3 air*/
        } else {
            VOL_ICE   = 0.0E+00;
            H2O_BOX_S = 0.0E+00;
            KG_ICE    = 0.0E+00;
        }

        /* Calculate NAT if relevant */
        if ( ( HNO3_BOX_S > 0.0E+00 ) && ( STATE_LOCAL >= 2 ) ) {
            HNO3_BOX_S = ( HNO3PP - PSATHNO3 ) / pressure_Pa;
            HNO3_BOX_S = std::max( HNO3_BOX_S, 0.0E+00 );

            /* Calculate m^3 NAT / m^3 air
             * HNO3_BOX_S is the number of moles of HNO3
             * which will be frozen into HNO3.3H2O (NAT) */
            VOL_NAT = HNO3_BOX_S * airDens * ( 1.0E+06 ) * MW_NAT / ( physConst::RHO_NAT * physConst::Na );
            /* Unit check:
             *      = mole       * molec/cm^3 * cm^3/m^3 * kg/mol / ( kg/m^3 ice         * molec/mole   )
             *      = m^3 ice / m^3 air*/
            KG_NAT  = HNO3_BOX_S * MW_NAT / physConst::Na * airDens * ( 1.0E+06 ) * Area;
            /* Unit check:
             *        = mole/mole * kg/mol * mole/molec * molec/cm^3 * cm^3/m^3 * m^2
             *        = kg / m */

        } else {
            HNO3_BOX_S = 0.0E+00;
            VOL_NAT    = 0.0E+00;
            KG_NAT     = 0.0E+00;
        }

        /* Calculate particle properties */
        if ( STATE_LOCAL < 2 ) {
            /* Zero'em all */
            KG_AER_BOX    = 0.0E+00;
            RAD_AER_BOX   = 0.0E+00;
            RHO_AER_BOX   = physConst::RHO_ICE;
            NDENS_AER_BOX = 0.0E+00;
        } else {
            VOL_TOT       = VOL_NAT + VOL_ICE;
            KG_AER_BOX    = KG_NAT + KG_ICE;
            RAD_AER_BOX   = MIN_RAD;
            NDENS_AER_BOX = pow( ( 3.0E+00 * VOL_TOT / ( 4.0E+00 * physConst::PI * MAX_NDENS )), 1.0E+00 / 3.0E+00 );
        }

        /* Prevent div zero */
        if ( SafeDiv( 1.0, VOL_TOT ) ) {
            RHO_AER_BOX = ( VOL_ICE * physConst::RHO_ICE + VOL_NAT * physConst::RHO_NAT ) / VOL_TOT;
        } else {
            RHO_AER_BOX = physConst::RHO_ICE;
        }

        /* Calculate SAD in m^2/m^3 */
        SAD_AER_BOX = 4.0E+00 * physConst::PI * RAD_AER_BOX * RAD_AER_BOX * NDENS_AER_BOX;

    } else {
        /* Solid PSCs not active */
        RAD_AER_BOX   = 0.0E+00;
        RHO_AER_BOX   = 1.0E+03;
        KG_AER_BOX    = 0.0E+00;
        NDENS_AER_BOX = 0.0E+00;
        SAD_AER_BOX   = 0.0E+00;
        HNO3_BOX_S    = 0.0E+00;
        H2O_BOX_S     = 0.0E+00;
    }

    /* Store in outer arrays */
    RAD_AER[0]   = RAD_AER_BOX;   /* [m]       */
    RHO_AER[0]   = RHO_AER_BOX;   /* [kg/m^3]  */
    KG_AER[0]    = KG_AER_BOX;    /* [kg/m]    */
    NDENS_AER[0] = NDENS_AER_BOX; /* [#/m^3]   */
    SAD_AER[0]   = SAD_AER_BOX;   /* [m^2/m^3] */

    /* Repartition NIT and HNO3 to kg NO3 */
    if ( LSOLIDPSC && IS_STRAT ) {

        /* Convert NAT from kg NAT to kg NO3 */
        Data[id_NAT]  = KG_NAT * MW_NIT / MW_NAT; /* [kg/m] */

        /* Remove ( kg NO2 as NAT ) from total kg NO3
         * then convert to kg HNO3 */
        Data[id_HNO3] = ( KG_NO3 - Data[id_NAT] ) / MW_NIT * physConst::Na / Area * 1.00E-06;
        /* Unit check:
         * mole/m * molec/mole / Area * m^3/cm^3 */

    }

    SOLIDFRAC[0] = 0.0E+00;
    SOLIDFRAC[1] = 0.0E+00;
    SOLIDFRAC[2] = 0.0E+00;
    SOLIDFRAC[3] = 0.0E+00;
    SOLIDFRAC[4] = 0.0E+00;
    SOLIDFRAC[5] = 0.0E+00;
    SOLIDFRAC[6] = 0.0E+00;

    if ( ( HNO3SUM > 0.0E+00 ) && ( SafeDiv( HNO3_BOX_S, HNO3SUM ) ) )
        SOLIDFRAC[1] = HNO3_BOX_S / HNO3SUM;
    
    if ( ( H2OSUM > 0.0E+00 ) && ( SafeDiv( H2O_BOX_S, H2OSUM ) ) )
        SOLIDFRAC[6] = H2O_BOX_S  / H2OSUM;

    if ( DBG ) {
        std::cout << "H2SO4 SOLID FRAC : " << SOLIDFRAC[0] << "\n";
        std::cout << "HNO3  SOLID FRAC : " << SOLIDFRAC[1] << "\n";
        std::cout << "HCl   SOLID FRAC : " << SOLIDFRAC[2] << "\n";
        std::cout << "HOCl  SOLID FRAC : " << SOLIDFRAC[3] << "\n";
        std::cout << "HBr   SOLID FRAC : " << SOLIDFRAC[4] << "\n";
        std::cout << "HOBr  SOLID FRAC : " << SOLIDFRAC[5] << "\n";
        std::cout << "H2O   SOLID FRAC : " << SOLIDFRAC[6] << "\n";
    }

    /* Now start liquid aerosol consideration */
    /* Start by assuming all non-solid H2O/HNO3 is gaseous */

    HNO3_BOX_G = HNO3SUM - HNO3_BOX_S;
    HNO3_BOX_L = 0.0E+00;
    H2O_BOX_G  = H2OSUM - H2O_BOX_S;
    H2O_BOX_L  = 0.0E+00;

    /* Calculate mixing ratios of other relevant species */
    H2SO4SUM = Data[id_SO4]   / airDens; //* INVAIR / MW_SO4;
    BrNO3SUM = Data[id_BrNO3] / airDens; //* INVAIR / MW_BrNO3;
    ClNO3SUM = Data[id_ClNO3] / airDens; //* INVAIR / MW_ClNO3;
    HOClSUM  = Data[id_HOCl]  / airDens; //* INVAIR / MW_HOCl;
    HClSUM   = Data[id_HCl]   / airDens; //* INVAIR / MW_HCl;
    HOBrSUM  = Data[id_HOBr]  / airDens; //* INVAIR / MW_HOBr;
    HBrSUM   = Data[id_HBr]   / airDens; //* INVAIR / MW_HBr;

    if ( DBG ) {
        std::cout << "H2SO4SUM : " << H2SO4SUM << "\n";
        std::cout << "H2OSUM   : " << H2OSUM   << "\n";
        std::cout << "HNO3SUM  : " << HNO3SUM  << "\n";
        std::cout << "BrNO3SUM : " << BrNO3SUM << "\n";
        std::cout << "ClNO3SUM : " << ClNO3SUM << "\n";
        std::cout << "HOClSUM  : " << HOClSUM  << "\n";
        std::cout << "HClSUM   : " << HClSUM   << "\n";
        std::cout << "HOBrSUM  : " << HOBrSUM  << "\n";
        std::cout << "HBrSUM   : " << HBrSUM   << "\n";
    }

    /* Consider gaseous H2SO4 to be unavailable for SLA */
    H2SO4_gFRAC = H2SO4_GASFRAC( temperature_K, Data[id_SO4] );
    H2SO4_BOX_L = H2SO4SUM * ( 1.0 - H2SO4_gFRAC );
    H2SO4_BOX_G = H2SO4SUM - H2SO4_BOX_L;
    AERFRAC[0]  = 1.0 - H2SO4_gFRAC;

    if ( DBG ) {
        std::cout << "H2SO4_GASFRAC : " << H2SO4_gFRAC << "\n";
        std::cout << "H2SO4_BOX_G   : " << H2SO4_BOX_G << "\n";
        std::cout << "H2SO4_BOX_L   : " << H2SO4_BOX_L << "\n";
    }

    /* Zero'em all */
    RHO_AER_BOX   = 1.0E+03;
    RAD_AER_BOX   = 0.0E+00;
    KG_AER_BOX    = 0.0E+00;
    NDENS_AER_BOX = 0.0E+00;
    SAD_AER_BOX   = 0.0E+00;
    VOL_SLA       = 0.0E+00;
    W_H2O         = 0.0E+00;
    W_H2SO4       = 0.0E+00;

    if ( !IS_STRAT ) {
        
        if ( DBG )
            std::cout << "Using JPL CTM data for GAMMA_BOX\n";

        /* Use JPL 10-06/Oslo CTM data for conventional sulfates/H2SO4 */
        GAMMA_BOX[ 0] = 1.0E-01;
        GAMMA_BOX[ 1] = 0.0E+00;
        GAMMA_BOX[ 2] = 0.0E+00;
        GAMMA_BOX[ 3] = 0.0E+00;
        GAMMA_BOX[ 4] = 3.0E-01;
        GAMMA_BOX[ 5] = 4.0E-01;
        GAMMA_BOX[ 6] = 9.0E-01;
        GAMMA_BOX[ 7] = 0.0E+00;
        GAMMA_BOX[ 8] = 0.0E+00;
        GAMMA_BOX[ 9] = 2.0E-01;
        GAMMA_BOX[10] = 0.0E+00;

    } else if ( H2SO4_BOX_L < 1.0E-15 ) {
        /* No aerosol */

        if ( DBG )
            std::cout << "No SLA aerosol\n";

        for ( unsigned int k = 0; k < 11; k++ )
            GAMMA_BOX[k] = 0.0E+00;

    } else {
        
        if ( DBG )
            std::cout << "Running TERNARY\n";

        if ( STATE_LOCAL == 0 ) {
            /* Allow binary H2SO4.nH2O only */
            TERNARY( temperature_K , pressure_Pa , H2OSUM         , H2SO4_BOX_L , \
                     0.0E+00       , HClSUM      , HOClSUM        , HBrSUM      , HOBrSUM , \
                     W_H2SO4       , W_H2O       , W_HNO3         , W_HCl       , W_HOCl  , W_HBr , W_HOBr , \
                     HNO3_GASFRAC  , HCl_GASFRAC , HOCl_GASFRAC   , HBr_GASFRAC , \
                     HOBr_GASFRAC  , VOL_SLA     , RHO_AER_BOX );

            /* For safety's safe, zero out HNO3 uptake */
            HNO3_GASFRAC = 1.0E+00;
            W_H2O += W_HNO3;
            W_HNO3 = 0.0E+00;
            HNO3_BOX_G = HNO3SUM - HNO3_BOX_S;
            HNO3_BOX_L = 0.0E+00;

        } else {
            /* Use only non-NAT HNO3 for STS */
            HNO3_BOX_G = HNO3SUM - HNO3_BOX_S;
            TERNARY( temperature_K , pressure_Pa , H2OSUM         , H2SO4_BOX_L , \
                     HNO3_BOX_G    , HClSUM      , HOClSUM        , HBrSUM      , HOBrSUM , \
                     W_H2SO4       , W_H2O       , W_HNO3         , W_HCl       , W_HOCl  , W_HBr , W_HOBr , \
                     HNO3_GASFRAC  , HCl_GASFRAC , HOCl_GASFRAC   , HBr_GASFRAC , \
                     HOBr_GASFRAC  , VOL_SLA     , RHO_AER_BOX );

            /* Partition HNO3 here for safety */
            HNO3_BOX_G = HNO3_BOX_G * HNO3_GASFRAC;
            HNO3_BOX_L = HNO3SUM - ( HNO3_BOX_G + HNO3_BOX_S );
        }
    
        if ( DBG ) {
            std::cout << "W_H2SO4 : " << W_H2SO4 << "\n";
            std::cout << "W_H2O   : " << W_H2O   << "\n";
            std::cout << "W_HNO3  : " << W_HNO3  << "\n";
            std::cout << "W_HCl   : " << W_HCl   << "\n";
            std::cout << "W_HOCl  : " << W_HOCl  << "\n";
            std::cout << "W_HBr   : " << W_HBr   << "\n";
            std::cout << "W_HOBr  : " << W_HOBr  << "\n";
        }

        /* Partition minor species */
        HCl_BOX_G  = HClSUM * HCl_GASFRAC;
        HCl_BOX_L  = HClSUM - HCl_BOX_G;
        HOCl_BOX_G = HOClSUM * HOCl_GASFRAC;
        HOCl_BOX_L = HOClSUM - HOCl_BOX_G;
        HBr_BOX_G  = HBrSUM * HBr_GASFRAC;
        HBr_BOX_L  = HBrSUM - HBr_BOX_G;
        HOBr_BOX_G = HOBrSUM * HOBr_GASFRAC;
        HOBr_BOX_L = HOBrSUM - HBr_BOX_G;

        /* Calculate SLA parameters (Grainger 1995) */
        SAD_AER_BOX = 1.0E+02 * SLA_VA * pow( VOL_SLA, 7.51E-01 ); /* [m^2/m^3] */
        /* Effective radius */
        RAD_AER_BOX = SLA_VR * SLA_RR * pow( VOL_SLA, 2.49E-01 );  /* [m]       */
        KG_AER_BOX  = RHO_AER_BOX * VOL_SLA * Area;                /* [kg/m]    */
        /* Unit check:
         *         = kg/m^3      * m^3/m^3 air * m^2 air */

        if ( DBG ) {
            std::cout << "VOL_SLA : " << VOL_SLA     << " [m^3/m^3]\n";
            std::cout << "RAD_SLA : " << RAD_AER_BOX << " [m]\n";
            std::cout << "SAD_SLA : " << SAD_AER_BOX << " [m^2/m^3]\n";
        }

        if ( VOL_SLA > 1.0E-30 ) {
            /* Approximate particles as spherical for calculation of aerosol 
             * number density */
            NDENS_AER_BOX = VOL_SLA * 3.0E+00 / (4.0E+00 * physConst::PI * RAD_AER_BOX * \
                                                            RAD_AER_BOX * RAD_AER_BOX );

            GAMMA_BOX = SLA_GAMMA( temperature_K         , pressure_Pa , \
                                   W_H2SO4               , \
                                   H2OSUM                , HClSUM      , HBrSUM   , \
                                   HOBrSUM               , ClNO3SUM    , BrNO3SUM ,
                                   RHO_AER_BOX * 1.0E-03 ,
                                   RAD_AER_BOX * 1.0E+02 );
        } else {
            /* Ignore SLA */
            for ( unsigned int k = 0; k < 11; k++ ) {
                GAMMA_BOX[k] = 0;
            }

        }

    }

    /* Store liquid fractions */
    /* Liquid H2O is removed from the sum, then it is assumed
     * that the pre-calculated solid H2O is taken out of this
     * liquid total */
    H2O_BOX_L = (9.809E+01/1.802E+01)*H2SO4_BOX_L * (W_H2O/W_H2SO4);
    H2O_BOX_L = std::max(0.0E+00, std::min(H2O_BOX_L - H2O_BOX_S, H2OSUM));
    H2O_BOX_G = std::max(0.0E+00, H2O_BOX_G - (H2O_BOX_L + H2O_BOX_S));

    AERFRAC[1] = 0.0E+00;
    AERFRAC[2] = 0.0E+00;
    AERFRAC[3] = 0.0E+00;
    AERFRAC[4] = 0.0E+00;
    AERFRAC[5] = 0.0E+00;
    AERFRAC[6] = 0.0E+00;

    if ( ( HNO3SUM > 0.0E+00 ) && ( SafeDiv( HNO3_BOX_L, HNO3SUM ) ) )
        AERFRAC[1] = HNO3_BOX_L / HNO3SUM;

    if ( ( HClSUM > 0.0E+00 ) && ( SafeDiv( HCl_BOX_L, HClSUM ) ) )
        AERFRAC[2] = HCl_BOX_L / HClSUM;

    if ( ( HOClSUM > 0.0E+00 ) && ( SafeDiv( HOCl_BOX_L, HOClSUM ) ) )
        AERFRAC[3] = HOCl_BOX_L / HOClSUM;

    if ( ( HBrSUM > 0.0E+00 ) && ( SafeDiv( HBr_BOX_L, HBrSUM ) ) )
        AERFRAC[4] = HBr_BOX_L / HBrSUM;

    if ( ( HOBrSUM > 0.0E+00 ) && ( SafeDiv( HOBr_BOX_L, HOBrSUM ) ) )
        AERFRAC[5] = HOBr_BOX_L / HOBrSUM;

    if ( ( H2OSUM > 0.0E+00 ) && ( SafeDiv( H2O_BOX_L, H2OSUM ) ) )
        AERFRAC[6] = H2O_BOX_L / H2OSUM;

    if ( DBG ) {
        std::cout << "H2SO4 AER FRAC : " << AERFRAC[0] << "\n";
        std::cout << "HNO3  AER FRAC : " << AERFRAC[1] << "\n";
        std::cout << "HCl   AER FRAC : " << AERFRAC[2] << "\n";
        std::cout << "HOCl  AER FRAC : " << AERFRAC[3] << "\n";
        std::cout << "HBr   AER FRAC : " << AERFRAC[4] << "\n";
        std::cout << "HOBr  AER FRAC : " << AERFRAC[5] << "\n";
        std::cout << "H2O   AER FRAC : " << AERFRAC[6] << "\n";
    }

    /* Convert sticking coefficients into premultiplying factors */
    KHET_COMMON = 0.25E+00 * sqrt( 8 * physConst::R * 1.0E+07 * temperature_K / physConst::PI );

    /* N2O5 + H2O/HCl */
    KHET_SPECIFIC= KHET_COMMON / sqrt(MW_N2O5 * 1.0E+03);
    KHETI_SLA[0] = GAMMA_BOX[0]*KHET_SPECIFIC;
    KHETI_SLA[1] = GAMMA_BOX[1]*KHET_SPECIFIC;

    /* ClNO3 + H2O/HCl/HBr */
    KHET_SPECIFIC= KHET_COMMON / sqrt(MW_ClNO3 * 1.0E+03);
    KHETI_SLA[2] = GAMMA_BOX[2]*KHET_SPECIFIC;
    KHETI_SLA[3] = GAMMA_BOX[3]*KHET_SPECIFIC;
    KHETI_SLA[4] = GAMMA_BOX[4]*KHET_SPECIFIC;
 
    /* BrNO3 + H2O/HCl */
    KHET_SPECIFIC= KHET_COMMON / sqrt(MW_BrNO3 * 1.0E+03);
    KHETI_SLA[5] = GAMMA_BOX[5]*KHET_SPECIFIC;
    KHETI_SLA[6] = GAMMA_BOX[6]*KHET_SPECIFIC;

    /* HOCl + HCl/HBr */
    KHET_SPECIFIC= KHET_COMMON / sqrt(MW_HOCl * 1.0E+03);
    KHETI_SLA[7] = GAMMA_BOX[7]*KHET_SPECIFIC;
    KHETI_SLA[8] = GAMMA_BOX[8]*KHET_SPECIFIC;

    /* HOBr + HBr/HCl */
    KHET_SPECIFIC= KHET_COMMON / sqrt(MW_HOBr * 1.0E+03);
    KHETI_SLA[9] = GAMMA_BOX[9]*KHET_SPECIFIC;
    KHETI_SLA[10]= GAMMA_BOX[10]*KHET_SPECIFIC;

    if ( DBG ) {
        std::cout << "KHETI_SLA[ 0] : " << KHETI_SLA[0] << " (GAMMA_BOX[ 0] : " << GAMMA_BOX[0] << ")\n";
        std::cout << "KHETI_SLA[ 1] : " << KHETI_SLA[1] << " (GAMMA_BOX[ 1] : " << GAMMA_BOX[1] << ")\n";
        std::cout << "KHETI_SLA[ 2] : " << KHETI_SLA[2] << " (GAMMA_BOX[ 2] : " << GAMMA_BOX[2] << ")\n";
        std::cout << "KHETI_SLA[ 3] : " << KHETI_SLA[3] << " (GAMMA_BOX[ 3] : " << GAMMA_BOX[3] << ")\n";
        std::cout << "KHETI_SLA[ 4] : " << KHETI_SLA[4] << " (GAMMA_BOX[ 4] : " << GAMMA_BOX[4] << ")\n";
        std::cout << "KHETI_SLA[ 5] : " << KHETI_SLA[5] << " (GAMMA_BOX[ 5] : " << GAMMA_BOX[5] << ")\n";
        std::cout << "KHETI_SLA[ 6] : " << KHETI_SLA[6] << " (GAMMA_BOX[ 6] : " << GAMMA_BOX[6] << ")\n";
        std::cout << "KHETI_SLA[ 7] : " << KHETI_SLA[7] << " (GAMMA_BOX[ 7] : " << GAMMA_BOX[7] << ")\n";
        std::cout << "KHETI_SLA[ 8] : " << KHETI_SLA[8] << " (GAMMA_BOX[ 8] : " << GAMMA_BOX[8] << ")\n";
        std::cout << "KHETI_SLA[ 9] : " << KHETI_SLA[9] << " (GAMMA_BOX[ 9] : " << GAMMA_BOX[9] << ")\n";
        std::cout << "KHETI_SLA[10] : " << KHETI_SLA[10] << " (GAMMA_BOX[10] : " << GAMMA_BOX[10] << ")\n";
    }

    RAD_AER[1]   = RAD_AER_BOX;   /* [m]       */
    RHO_AER[1]   = RHO_AER_BOX;   /* [kg/m^3]  */
    KG_AER[1]    = KG_AER_BOX;    /* [kg/m]    */
    NDENS_AER[1] = NDENS_AER_BOX; /* [#/m^3]   */
    SAD_AER[1]   = SAD_AER_BOX;   /* [m^2/m^3] */

    if ( DBG ) {
        std::cout << "SPA :\n";
        std::cout << "RAD : " << RAD_AER[0]   << "[m]\n";
        std::cout << "RHO : " << RHO_AER[0]   << "[kg/m^3]\n";
        std::cout << "KG  : " << KG_AER[0]    << "[kg/m]\n";
        std::cout << "DENS: " << NDENS_AER[0] << "[#/m^3]\n";
        std::cout << "SAD : " << SAD_AER[0]   << "[m^2/m^3]\n";
        std::cout << "SLA :\n";
        std::cout << "RAD : " << RAD_AER[1]   << "[m]\n";
        std::cout << "RHO : " << RHO_AER[1]   << "[kg/m^3]\n";
        std::cout << "KG  : " << KG_AER[1]    << "[kg/m]\n";
        std::cout << "DENS: " << NDENS_AER[1] << "[#/m^3]\n";
        std::cout << "SAD : " << SAD_AER[1]   << "[m^2/m^3]\n";
    }

    return STATE_LOCAL;

} /* End of STRAT_AER */

std::vector<double> SLA_GAMMA( const double T_K         , const double P_Pa     , \
                               const double WT_FRC      ,                         \
                               const double H2OSUM      , const double HClSUM   , \
                               const double HBrSUM      , const double HOBrSUM  , \
                               const double ClNO3SUM    , const double BrNO3SUM , \
                               const double RHO         , const double ARAD )
{

    /* DESCRIPTION: Calculates 11 different sticking 
     * coefficients on the surface of local stratospheric
     * liquid aerosols, relevant to each of the 11 reactions
     * listed in Kirner's paper. */

    /* INPUTS:
     * - T_K      : Temperature in K
     * - P_Pa     : Pressure in Pa
     * - WT_FRC   : Weight fraction of H2SO4 (kg/kg)
     * - H2OSUM   : H2O mixing ratio
     * - HClSUM   : HCl mixing ratio
     * - HBrSUM   : HBr mixing ratio
     * - HOBrSUM  : HOBr mixing ratio
     * - ClNO3SUM : ClNO3 mixing ratio
     * - BrNO3SUM : BrNO3 mixing ratio
     * - RHO      : STS density in g/cm^3
     * - ARAD     : SLA radius in cm
     *
     * OUTPUTS:
     * RXNGAMMA   : PRemultyplying factors */

    /* Premultiplying factors */
    std::vector<double> RXNGAMMA( 11, 0.0E+00 );

    /* Weight percentage H2SO4 (100*kg/kg) */
    double WT;

    /* Partial pressure of H2O in Pa */
    double H2OPP;
    
    /* Partial pressures of HCl, HBr, HOBr, ClONO2, BrONO3 in atm */
    double HClPP, HBrPP, HOBrPP, ClNO3PP, BrNO3PP;

    /* Water vapor saturation pressure in Pa */
    double PSATH2O;

    /* Activity of water */
    double ACTH2O;

    /* Molality of H2SO4 (mol H2SO4/kg solvent) */
    double MOLAL;

    /* Mass of H2SO4 */
    double M_H2SO4;

    /* HOBr parameters */
    double c_HOBr, SHOBr, HHOBr, DHOBr, kHOBr_HCl, GHOBrrxn, lHOBr, fHOBr, gHOBr_HCl;

    /* HOCl parameters */
    double c_HOCl, SHOCl, HHOCl, DHOCl, kHOCl_HCl, GHOClrxn, lHOCl, fHOCl, gHOCl_HCl;

    /* ClNO3 parameters */
    double c_ClNO3, SClNO3, HClNO3, DClNO3, GClNO3rxn, lClNO3, fClNO3, gClNO3, gClNO3_HCl, gClNO3_H2O;

    /* N2O5 parameters */
    std::vector<double> AK( 3, 0.0 );

    /* Other parameters */
    double kH2O, kH, khdr, GbH2O, HHCl, MHCl, kHCl, GbHCl, Gs, FHCl, Gsp, GbHClp, Gb, khydr, kII, k_dl;

    /* Interim variables */
    double X, A, H, T_THRESHOLD, aH;
    const double MAX_T_DIFF = 6.0E+00;

    /* Control whether to run calculations */
    bool HClOK, HOBrOK;

    /* Initialize */

    WT          = 0.0E+00;
    H2OPP       = 0.0E+00;
    HClPP       = 0.0E+00;
    HBrPP       = 0.0E+00;
    HOBrPP      = 0.0E+00;
    ClNO3PP     = 0.0E+00;
    BrNO3PP     = 0.0E+00;
    PSATH2O     = 0.0E+00;
    ACTH2O      = 0.0E+00;
    MOLAL       = 0.0E+00;
    M_H2SO4     = 0.0E+00;
    c_HOBr      = 0.0E+00;
    SHOBr       = 0.0E+00;
    HHOBr       = 0.0E+00;
    DHOBr       = 0.0E+00;
    kHOBr_HCl   = 0.0E+00;
    GHOBrrxn    = 0.0E+00;
    lHOBr       = 0.0E+00;
    fHOBr       = 0.0E+00;
    gHOBr_HCl   = 0.0E+00;
    c_HOCl      = 0.0E+00;
    SHOCl       = 0.0E+00;
    HHOCl       = 0.0E+00;
    DHOCl       = 0.0E+00;
    kHOCl_HCl   = 0.0E+00;
    GHOClrxn    = 0.0E+00;
    lHOCl       = 0.0E+00;
    fHOCl       = 0.0E+00;
    gHOCl_HCl   = 0.0E+00;
    c_ClNO3     = 0.0E+00;
    SClNO3      = 0.0E+00;
    HClNO3      = 0.0E+00;
    DClNO3      = 0.0E+00;
    GClNO3rxn   = 0.0E+00;
    lClNO3      = 0.0E+00;
    fClNO3      = 0.0E+00;
    gClNO3      = 0.0E+00;
    gClNO3_HCl  = 0.0E+00;
    gClNO3_H2O  = 0.0E+00;
    kH2O        = 0.0E+00;
    kH          = 0.0E+00;
    khdr        = 0.0E+00;
    GbH2O       = 0.0E+00;
    HHCl        = 0.0E+00;
    MHCl        = 0.0E+00;
    kHCl        = 0.0E+00;
    GbHCl       = 0.0E+00;
    Gs          = 0.0E+00;
    FHCl        = 0.0E+00;
    Gsp         = 0.0E+00;
    GbHClp      = 0.0E+00;
    Gb          = 0.0E+00;
    khydr       = 0.0E+00;
    kII         = 0.0E+00;
    k_dl        = 0.0E+00;
    X           = 0.0E+00;
    A           = 0.0E+00;
    H           = 0.0E+00;
    T_THRESHOLD = 0.0E+00;
    aH          = 0.0E+00;

    HClOK  = 0;
    HOBrOK = 0;


    /* Code begins here ! */

    /* Compute the saturation pressure of H2O, water partial 
     * pressure and water acitivity */
    PSATH2O = physFunc::pSat_H2Os( T_K );
    H2OPP   = H2OSUM * P_Pa;
    ACTH2O  = std::max( (H2OPP / PSATH2O), 1.0E+00 );

    /* Calculate molality of solution */
    WT = 1.00E+02 * WT_FRC; /* Convert from fraction to % */
    MOLAL = 1.0E+03 * ( WT / 9.80E+01 / ( 1.00E+02 - WT ) );

    /* Parameters for H2SO4 solution
     * The solution density is calculated earlier, including
     * contributions from HNO3. This code treats it as a binary
     * solution - so far this is just a kludge. Need to update
     * all this code to acknowledge the presence of at least
     * HNO3 (e.g. X is still calculated based on pure H2O
     * solvent!) */

    /* Molality (mol H2SO4/kg solvent) */
    M_H2SO4     = RHO * WT / 9.8E+01;
    X           = WT / ( WT + ( 1.00E+02 - WT ) * 9.80E+01 / 1.80E+01 );
    A           = 1.6950E+02 + 5.18E+00 * WT - 8.25E-02 * WT * WT + 3.27E-03 * WT * WT * WT;
    T_THRESHOLD = 1.4411E+02 + 1.66E-01 * WT - 1.50E-02 * WT * WT + 2.18E-04 * WT * WT * WT;

    if ( ( T_K - T_THRESHOLD ) > MAX_T_DIFF )
         H = A * pow( T_K, -1.43E+00 ) * exp( 4.48E+02 / ( T_K - T_THRESHOLD ) );
    else
         H = A * pow( T_K, -1.43E+00 ) * exp( 4.48E+02 / MAX_T_DIFF );

    aH = exp( 6.051E+01 - 9.50E-02 * WT + 7.70E-03 * WT * WT - 1.61E-05 * WT * WT * WT \
             - (  1.76E+00   + 2.52E-04 * WT * WT ) * sqrt( T_K )                      \
             + ( -8.0589E+02 + 2.5305E+02 * pow( WT, 0.076 ) ) / sqrt( T_K ) );

    /* All the following partial pressures are computed in atm */
    HClPP   = HClSUM   * P_Pa / physConst::ATM;
    ClNO3PP = ClNO3SUM * P_Pa / physConst::ATM;
    BrNO3PP = BrNO3SUM * P_Pa / physConst::ATM;
    HOBrPP  = HOBrSUM  * P_Pa / physConst::ATM;
 
    /* Should we bother running calculations? */
    HClOK  = ( HClPP  > 1.00E-30 );
    HOBrOK = ( HOBrPP > 1.00E-30 );

    /* Reaction 1. N2O5 + H2O (hydrolysis of N2O5) */
    AK[0] = -2.55265E+01 - 1.33188E-01 * WT + 9.30840E-03 * WT * WT - 9.01940E-05 * WT * WT * WT;
    AK[1] =  9.28376E+03 + 1.15345E+02 * WT - 5.19258E+00 * WT * WT + 4.83464E-02 * WT * WT * WT;
    AK[2] = -8.51801E+05 - 2.21912E+04 * WT + 7.66916E+02 * WT * WT - 6.85427E+00 * WT * WT * WT;
    RXNGAMMA[0] = exp( AK[0] + AK[1]/T_K + AK[2]/( T_K * T_K ) );
      
    /* Reaction 2. N2O5 + HCl */
    /* JPL 10-06 suggests near-zero gamma */
    RXNGAMMA[1] = 1.00E-80;

    
    /* Reactions 3/4. ClNO3 + H2O/HCl */
    /* Now od only if HCl concentrations are large enough
     * (to avoid div-by-zero errors) */
    if ( HClOK ) {
         c_ClNO3     =  1.4740E+03 * sqrt( T_K );
         SClNO3      =  3.06E-01 + 2.40E+01 / T_K;
         HClNO3      =  1.60E-06 * exp( 4.71E+03 / T_K ) * exp( -SClNO3 * M_H2SO4 );
         DClNO3      =  5.00E-08 * T_K / H;
         kH2O        =  1.95E+10 * exp(-2.80E+03 / T_K );
         kH          =  1.22E+12 * exp(-6.20E+03 / T_K );
         khydr       =  kH2O * ACTH2O + kH * aH * ACTH2O;
         GbH2O       =  4.0E+00 * HClNO3 * 8.20E-02 * T_K * sqrt( DClNO3 * khydr ) / c_ClNO3;
         HHCl        =  (9.4E-02 - 6.1E-01 * X + 1.2E+00 * X * X ) \
                        * exp(-8.68E+00 + (8.515E+03 - 1.0718E+04 * pow( X, 0.7 ) ) / T_K );
         MHCl        =  HHCl * HClPP;
         kHCl        =  7.90E+11 * aH * DClNO3 * MHCl;
         lClNO3      =  sqrt( DClNO3 / (khydr + kHCl) );
         if ( lClNO3 > (1.00E+05 * ARAD) ) {
            /* Limiting rate */
            fClNO3   =  ARAD / (3.0E+00 * lClNO3);
         }
         else {
            fClNO3   =  1.00E+00 / tanh(ARAD/lClNO3) - lClNO3/ARAD;
         }
         GClNO3rxn   =  fClNO3 * GbH2O * sqrt( 1.00E+00 + kHCl/khydr );
         GbHCl       =  GClNO3rxn * kHCl / (kHCl+ khydr);
         Gs          =  6.612E+01 * exp( -1.374E+03 / T_K ) * HClNO3 * MHCl;
         FHCl        =  1.00E+00 / ( 1.00E+00 + 6.12E-01 * (Gs+GbHCl) * ClNO3PP/ HClPP );
         Gsp         =  FHCl*Gs;
         GbHClp      =  FHCl*GbHCl;
         Gb          =  GbHClp + GClNO3rxn * khydr / (kHCl+ khydr);
         gClNO3      =  1.00E+00 / (1.00E+00 + 1.00E+00 / (Gsp + Gb) );
         gClNO3_HCl  =  gClNO3 * (Gsp + GbHClp) / (Gsp + Gb);
         gClNO3_H2O  =  gClNO3 - gClNO3_HCl;
   
         RXNGAMMA[2] =  gClNO3_H2O;
         RXNGAMMA[3] =  gClNO3_HCl;
    } else {
         RXNGAMMA[2] = 1.00E-80;
         RXNGAMMA[3] = 1.00E-80;
    }

    /* Reaction 5. ClNO3 + HBr */
    /* Not present in JPL 10-06 for H2SO4 */
    RXNGAMMA[4] = 1.00E-80;

    /* Reaction 6. BrNO3 + H2O */
    RXNGAMMA[5] = 1.00E+00 / ( 1.00E+00 / 8.00E-01 + 1.00E+00 / ( exp( 2.92E+01 - 4.00E-01 * WT ) + 1.1E-01 ) );
      
    /* Reaction 7. BrNO3 + HCl */
    RXNGAMMA[6] = 9.00E-01; /* JPL 10-06 */

    /* Reaction 8. HOCl + HCl */
    if ( HClOK ) {
        c_HOCl    =  physFunc::thermalSpeed( T_K, 5.246E+01 );
        SHOCl     =  7.76E-02 + 5.918E+01 / T_K;
        HHOCl     =  1.91E-06 * exp( 5.8624+03 / T_K ) * exp( -SHOCl * M_H2SO4 );
        DHOCl     =  6.40E-08 * T_K / H;
        kHOCl_HCl =  1.25E+09 * aH * DHOCl * MHCl;
        lHOCl     =  sqrt( DHOCl/kHOCl_HCl );
        if ( lHOCl > (1.00E+05*ARAD) ) {
           /* Limiting rate */
           fHOCl  = ARAD / ( 3.00E+00*lHOCl );
        } else {
           fHOCl  =  1.00E+00 / tanh(ARAD/lHOCl) - lHOCl/ARAD;
        }
        GHOClrxn  =  4.00E+00 * HHOCl * 8.20E-02 * T_K * sqrt( DHOCl * kHOCl_HCl ) / c_HOCl;
        if ( fHOCl == 0.0E+00 )
           gHOCl_HCl =  1.00E-80;
        else
           gHOCl_HCl =  1.00E+00 / ( 1.00E+00 + 1.00E+00 / ( fHOCl * GHOClrxn * FHCl ) );
    } else
        gHOCl_HCl = 1.00E-80;
   
     RXNGAMMA[7] = gHOCl_HCl;

     /* Reaction 9. HOCl + HBr */
     /* Not yet implemented for STS; JPL 10-06 suggests complex
      * relationship, not yet sufficiently well understood or
      * parameterized for the purposes of simulation. Ignore for now */
     RXNGAMMA[8] = 1.00E-80;

     /* Reaction 10. HOBr + HCl */
     if ( ( HClOK ) && ( HOBrOK ) ) {
         c_HOBr    =  physFunc::thermalSpeed( T_K , 9.691E+01 );
         SHOBr     =  7.76E-02 + 5.918E+01 / T_K;
         HHOBr     =  exp( -9.86E+00 + 5.427E+03 / T_K );
         DHOBr     =  1.00E-08;
         kII       =  exp( 1.54E+02 - 1.63E+00 * WT ) * exp( -( 3.85E+04 - 4.78E+02 * WT ) / T_K );
         k_dl      =  7.50E+14 * ( DHOBr * ARAD * 1.00E+07 );
         if ( kII > k_dl )
             kII   = k_dl;

         kHOBr_HCl =  kII * HHOBr * HOBrPP;
         
         GHOBrrxn  =  4.00E+00 * HHCl * 8.20E-02 * T_K * sqrt( DHOBr *kHOBr_HCl ) / c_HOBr;
         lHOBr     =  sqrt( DHOBr/kHOBr_HCl );
         if ( lHOBr > ( 1.00E+03 * ARAD ) ) {
            /* Limiting rate */
            fHOBr  = ARAD / ( 3.0E+00 * lHOBr );
         }
         else
            fHOBr  =  1.00E+00 / tanh(ARAD/lHOBr) - lHOBr/ARAD;

         gHOBr_HCl =  1.00E+00 / (1.0E+00 + 1.0E+00 / (fHOBr*GHOBrrxn) );
         RXNGAMMA[9] = gHOBr_HCl;
     } else
         RXNGAMMA[9] = 1.00E-80;

      /* Reaction 11. HOBr + HBr */
     /* Data from JPL limited; ignore for now */
    RXNGAMMA[10] = 1.00E-80;

    return RXNGAMMA;    


} /* End of SLA_GAMMA */

void TERNARY( const double TIN_K    , const double PIN_Pa  , const double H2OSUM_IN , \
              const double H2SO4SUM , const double HNO3SUM , const double HClSUM    , \
              const double HOClSUM  , const double HBrSUM  , const double HOBrSUM   , \
              double &W_H2SO4       , double &W_H2O        , double &W_HNO3         , \
              double &W_HCl         , double &W_HOCl       , double &W_HBr          , \
              double &W_HOBr        , double &HNO3GASFRAC  , double &HClGASFRAC     , \
              double &HOClGASFRAC   , double &HBrGASFRAC   , double &HOBrGASFRAC    , \
              double &SLA_VOL       , double &SLA_RHO )
{

    /* DESCRIPTION: Subroutine TERNARY calculates the composition of SSA/STS 
     *  aerosols using a paramaterization from Carslaw et al. "A Thermodynamic
     *  Model of the System HCl-HNO3-H2SO4-H2O, Including Solubilities of HBr,
     *  from <200 to 328 K". The bulk of this code was taken directly from the
     *  Global Modeling Initiative implementation by David Considine. */

    /* INPUTS:
     * - TIN_K     : Temperature (K)
     * - PIN_Pa    : Pressure (Pa)
     * - H2OSUM_IN : Total H2O mixing ratio
     * - H2SO4SUM  : Liquid H2SO4 mixing ratio
     * - HNO3SUM   : Total HNO3 mixing ratio
     * - HClSUM    : Total HCl mixing ratio
     * - HOClSUM   : Total HOCl mixing ratio
     * - HBrSUM    : Total HBr mixing ratio
     * - HOBrSUM   : Total HOBr mixing ratio
     *
     * OUTPUTS:
     * - W_H2SO4     : kg H2SO4/kg SLA
     * - W_H2O       : kg H2O  /kg SLA
     * - W_HNO3      : kg HNO3 /kg SLA
     * - W_HCl       : kg HCl  /kg SLA
     * - W_HOCl      : kg HOCl /kg SLA
     * - W_HBr       : kg HBr  /kg SLA
     * - W_HOBr      : kg HOBr /kg SLA
     * - HNO3GASFRAC : Gas fraction HNO3
     * - HClGASFRAC  : Gas fraction HCl
     * - HOClGASFRAC : Gas fraction HOCl
     * - HBrGASFRAC  : Gas fraction HBr
     * - HOBrGASFRAC : Gas fraction HOBr
     * - SLA_VOL     : Aerosol volume (m3/m3)
     * - SLA_RHO     : Aer. mass density (kg/m3) */

    /* Derived inputs */
    double H2OSUM, T_K, P_Pa;

    /* Partial pressures */
    double PATMH2O, PATMHNO3, PATMHCl, PATMHOCl, PATMHBr, PATMHOBr;

    /* Molar densities (mol/m3) */
    double MOLDENS_H2SO4;

    /* Mass totals */
    double M_H2SO4, M_HNO3, M_HCl, M_HOCl, M_HBr, M_HOBr;

    /* Binary solutions denoted with BIN */
    /* Mole fractions */
    double X_H2SO4_BIN, X_HNO3_BIN;

    /* Mass fractions */
    double M_H2SO4_BIN, M_HNO3_BIN;

    /* Effective Henry's Law coefficients */
    double H_H2SO4_BIN, H_HNO3_BIN, H_HCl, H_HOCl, H_HBr, H_HOBr;

    /* Frost point */
    double T_ICE;

    /* Equilibrium vapor pressures */
    double PVAP_HNO3, PVAP_HCl, PVAP_HBr, PVAP_HOBr;

    /* Transitional variables */
    double PR, TR, TT;

    /* Derived parameters */
    double A, B, C, PHI;


    const double R_ATM = physConst::R / physConst::ATM;

    const std::vector<double> QN{  1.45734E+01,  6.15994E-02, -1.14895E+00,  6.916930E-01, -9.88630E-02, \
                                   5.15790E-03,  1.23472E-01, -1.15574E-01,  1.101130E-02,  9.79140E-03 };
    const std::vector<double> QS{  1.44700E+01,  6.38795E-02, -3.29597E+00,  1.778224E+00, -2.23244E-01, \
                                   8.64860E-03,  5.36695E-01, -3.35164E-01,  2.651530E-02,  1.57550E-02 };
    const std::vector<double> KN{ -3.91360E+01,  6.35840E+03,  8.32900E+01, -1.765000E+04,  1.98530E+02, \
                                  -1.19480E+04, -2.84690E+01 };
    const std::vector<double> KS{ -2.16610E+01,  2.72420E+03,  5.18100E+01, -1.573200E+04,  4.70040E+01, \
                                  -6.96900E+03, -4.61830E+00 };

    /* Routine only valid for certain limits */
    H2OSUM = std::max( H2OSUM_IN, 5.0E-07 );
    P_Pa   = std::max( PIN_Pa,5.0E+02 );
    T_K    = TIN_K;

    /* Calculate partial pressure of H2O & HNO3
     * P_Pa is in Pa */
    PATMH2O  = H2OSUM  * P_Pa / physConst::ATM;

    /* Carslaw only valid for 2e-5 < PPH2O < 2e-3 (hPa) */
    PATMH2O = std::max(PATMH2O,1.9738465e-8);
    PATMH2O = std::min(PATMH2O,1.9738465e-6);

    PATMHNO3 = HNO3SUM * P_Pa / physConst::ATM;
    PATMHCl  = HClSUM  * P_Pa / physConst::ATM;
    PATMHOCl = HOClSUM * P_Pa / physConst::ATM;
    PATMHBr  = HBrSUM  * P_Pa / physConst::ATM;
    PATMHOBr = HOBrSUM * P_Pa / physConst::ATM;

    /* Moles of H2SO4 per m3 air */
    MOLDENS_H2SO4 = P_Pa * H2SO4SUM / ( physConst::R * T_K );

    /* Nucleation temperature of ice */
    T_ICE = 2.66870E+03 / (1.04310E+01 - ( log( PATMH2O ) + log( 7.60E+02) ) / log(1.0E+01) );

    /* Pressure relation */
    PR = log(PATMH2O) + 1.84E+01;
      
    /* Therefore if temperature lower, set to T_ICE - 3 */
    if ( T_K < (T_ICE - 3.0E+00) ) 
         T_K = T_ICE-3.0E+00;

    if ( T_K < 185.0E+00) 
         T_K = 185.0E+00;

    TT = T_K * R_ATM * MOLDENS_H2SO4;
      
    /* Temperature relation */
    TR = 1.00E+04 / T_K - 4.34782608E+01;

    /* Determine H2SO4/H2O pure solution concentration
     * Mole fraction of H2SO4 in binary solution */

    X_H2SO4_BIN = 1.00E+00 / ( 2.00 * (KS[2]+KS[3]/T_K) ) * \
                  (-KS[0]-KS[1]/T_K - pow( (KS[0]+KS[1]/T_K) * (KS[0]+KS[1]/T_K) \
                   - 4.00E+00 * (KS[2]+KS[3]/T_K) * (KS[4]+KS[5]/T_K + KS[6]*log(T_K)-log(PATMH2O)), 0.50 ) );

    /* Molality (mol H2SO4/kg H2O) in binary solution */
    M_H2SO4_BIN = 5.551E+01 * X_H2SO4_BIN / ( 1.0E+00 - X_H2SO4_BIN );

    if ( (T_K <= 215.0E+00) && (PATMHNO3 > 0.0E+00) ) {

        /* Determine HNO3/H2SO4/H2O solution composition */
        H_H2SO4_BIN = exp( QS[0] + QS[1]*TR*TR \
                + ( QS[2] + QS[3]*TR + QS[4]*TR*TR + QS[5]*TR*TR*TR )*PR \
                + ( QS[6] + QS[7]*TR + QS[8]*TR*TR )*PR*PR \
                + ( QS[9]*TR )*PR*PR*PR );
        X_HNO3_BIN = 1.00E+00 / ( 2.00 * (KN[2]+KN[3]/T_K) ) * \
                     (-KN[0]-KN[1]/T_K - pow( (KN[0]+KN[1]/T_K) * (KN[0]+KN[1]/T_K) \
                      -4.0E+00 * (KN[2]+KN[3]/T_K) * (KN[4]+KN[5]/T_K + KN[6]*log(T_K) - log(PATMH2O)), 0.50 ) );
        /* Molality */
        M_HNO3_BIN = 5.551E+01 * X_HNO3_BIN / ( 1.0E+00 - X_HNO3_BIN );
        H_HNO3_BIN = exp( QN[0] + QN[1]*TR*TR \
                + ( QN[2] + QN[3]*TR + QN[4]*TR*TR + QN[5]*TR*TR*TR )*PR \
                + ( QN[6] + QN[7]*TR + QN[8]*TR*TR)*PR*PR \
                + ( QN[9]*TR )*PR*PR*PR );

        A = (  TT*H_HNO3_BIN*M_HNO3_BIN*M_HNO3_BIN - TT*H_H2SO4_BIN*M_HNO3_BIN*M_H2SO4_BIN - \
               2.0*M_HNO3_BIN*M_HNO3_BIN*M_H2SO4_BIN + M_HNO3_BIN*M_H2SO4_BIN*M_H2SO4_BIN + \
               H_HNO3_BIN*M_HNO3_BIN*M_H2SO4_BIN*PATMHNO3 - H_H2SO4_BIN*M_H2SO4_BIN*M_H2SO4_BIN*PATMHNO3 )\
            /( M_HNO3_BIN*M_HNO3_BIN - M_HNO3_BIN*M_H2SO4_BIN );
        B = M_H2SO4_BIN * ( -2.0*TT*H_HNO3_BIN*M_HNO3_BIN + TT*H_H2SO4_BIN*M_H2SO4_BIN + \
                            M_HNO3_BIN*M_H2SO4_BIN - H_HNO3_BIN*M_H2SO4_BIN*PATMHNO3)\
            /( M_HNO3_BIN - M_H2SO4_BIN );
        C = ( TT*H_HNO3_BIN*M_HNO3_BIN*M_H2SO4_BIN*M_H2SO4_BIN )\
            /(M_HNO3_BIN - M_H2SO4_BIN);
        
        PHI = atan( sqrt( 4.0 * ( A*A - 3.0*B )*( A*A - 3.0*B )*( A*A - 3.0*B ) - \
                       ( -2.0*A*A*A + 9.0*A*B - 27.0*C )*( -2.0*A*A*A + 9.0*A*B - 27.0*C ) )\
                /( -2.0*A*A*A + 9.0*A*B - 27.0*C));

        if (PHI < 0.0E+00) 
           PHI += physConst::PI;
         
        M_H2SO4 = -1.0/3.0 * ( A + 2.0 * sqrt(A*A - 3.0*B) * cos( (physConst::PI + PHI)/3.0) );
        M_HNO3  = M_HNO3_BIN * (1.0 - M_H2SO4/M_H2SO4_BIN );
        W_H2SO4 = M_H2SO4 * 9.8076E-02 / (1.00E+00 + M_H2SO4 * 9.8076E-02 + M_HNO3 * 6.3012E-02);

        /* Check for low H2SO4 */
        if (M_H2SO4 < 0.0E+00) {
            M_H2SO4 = 0.0e+00;
            M_HNO3 = M_HNO3_BIN;
            W_H2SO4 = 0.0e+00;
        }

        PVAP_HNO3 = M_HNO3 / ( H_HNO3_BIN*M_HNO3/(M_HNO3 + M_H2SO4) + H_H2SO4_BIN*M_H2SO4/(M_HNO3 + M_H2SO4) );
        W_HNO3    = (M_HNO3*6.3012E-02) / (1.00E+00 + M_H2SO4* 9.8076E-02 + M_HNO3*6.3012E-02);

        HNO3GASFRAC = ( 1.00E+00 - (PATMHNO3 - PVAP_HNO3)/PATMHNO3 );
        M_HNO3      = std::max( M_HNO3,0.0E+00 );
        W_HNO3      = std::max( W_HNO3,0.0E+00 ); 

    } else {

        /* Solution is pure H2SO4/H2O */
        M_H2SO4 = M_H2SO4_BIN;
        M_HNO3 = 0.0E+00;
        W_H2SO4 = M_H2SO4_BIN*9.8076E-02 / (1.00E+00 + M_H2SO4_BIN*9.8076E-02);
        W_HNO3 = 0.0E+00;
        PVAP_HNO3 = 0.0E+00;
        HNO3GASFRAC = 1.0E+00;

    }
      
    /* Handle HCl (Luo et al., Vapor pressures of
     * H2SO4/HNO3/HCl/HBr/H2O solutions to low stratospheric
     * temperatures, 1995) 
     * Fit is in mbar = hPa */
    if ( PATMHCl > 0.00E+00 ) {

        H_HCl = exp( - ( 2.10E+01 + 4.6610E+01*W_HNO3 + 4.0690E+00*W_H2SO4 - 4.8370E+00*sqrt(W_HNO3) \
                       + 2.1860E+00*sqrt(W_H2SO4) - 6.30E+01*W_HNO3*W_HNO3 - 4.017E+01*W_HNO3*W_H2SO4 \
                       - 1.571E+00*W_H2SO4*W_H2SO4 ) - \
                    1.0/T_K * ( -7.437E+03 - 8.3278E+03*W_HNO3 + 1.30090E+03*W_H2SO4 \
                                +1.0872E+03*sqrt(W_HNO3) - 2.4271E+02*sqrt(W_H2SO4)  \
                                +1.8749E+04*W_HNO3*W_HNO3 + 1.85E+04*W_HNO3*W_H2SO4  \
                                +5.632E+03*W_H2SO4*W_H2SO4 ) - \
                    log(W_HNO3 + 6.10E-01*W_H2SO4) - \
                    log(3.6461E+01/(1.00E+03 + 9.8076E+01*M_H2SO4 + 6.3012E+01*M_HNO3)) ) * physConst::ATM * 1.00E-02;
        M_HCl = (1.00E+00/R_ATM/T_K * PATMHCl)/(MOLDENS_H2SO4/M_H2SO4 + 1.00E+00/R_ATM/T_K/H_HCl);
        W_HCl = M_HCl*3.6461E+01/(1.00E+03 + 9.8076E+01*M_H2SO4 + 6.3012E+01*M_HNO3);
        PVAP_HCl = M_HCl/H_HCl;
        HClGASFRAC = 1.00E+00 - (PATMHCl - PVAP_HCl)/PATMHCl;

    } else {

        W_HCl = 0.00e+00;
        HClGASFRAC = 1.00E+00;

    }

    /* Now HOCl */
    if ( PATMHOCl > 0.0E+00 ) {

        H_HOCl = exp( 6.49460E+00 - (-4.1070E-02 + 5.456E+01/T_K)*(M_H2SO4+M_HNO3) - \
                      5.862E+03*(1.00E+00/298.15 - 1.00E+00/T_K) );
        M_HOCl = (1.00E+00/R_ATM/T_K*PATMHOCl)/(MOLDENS_H2SO4/M_H2SO4 + 1.00E+00/R_ATM/T_K/H_HOCl);
        W_HOCl = M_HOCl*5.246E+01/(1.0E+03 + 9.8076E+01*M_H2SO4 + 6.3012E+01*M_HNO3);
        /* Realistically expect no gas phase removal */
        HOClGASFRAC = 1.00E+00 ;

    } else {

       W_HOCl = 0.0E+00;
       HOClGASFRAC = 1.00E+00;

    }
 
    /* Now HBr (Luo et al., Vapor pressures of
     * H2SO4/HNO3/HCl/HBr/H2O solutions to low stratospheric
     * temperatures, 1995) */
    if ( PATMHBr > 0.00E+00 ) {
          
        H_HBr = exp( - ( 1.783E+01 + 1.02E+00*W_HNO3 - 1.08E+00*W_H2SO4 + 3.9E+00*sqrt(W_HNO3) \
                       + 4.38E+00*sqrt(W_H2SO4) - 8.87E+00*W_HNO3*W_HNO3 - 1.70E+01*W_HNO3*W_H2SO4 \
                       + 3.73E+00*W_H2SO4*W_H2SO4 ) - \
                    1.0/T_K * ( -8.2205E+03 - 3.6276E+02*W_HNO3 + 6.5893E+02*W_H2SO4 \
                                -9.140E+02*sqrt(W_HNO3) - 9.553E+02*sqrt(W_H2SO4)    \
                                +9.9766E+03*W_HNO3*W_HNO3 + 1.97785E+05*W_HNO3*W_H2SO4\
                                +7.680E+03*W_H2SO4*W_H2SO4 ) - \
                    log(W_HNO3 + 4.10E-01*W_H2SO4) - \
                    log(3.6461E+01/(1.00E+03 + 9.8076E+01*M_H2SO4 + 6.3012E+01*M_HNO3)) ) * physConst::ATM * 1.00E-02;
        M_HBr = (1.00E+00/R_ATM/T_K*PATMHBr)/(MOLDENS_H2SO4/M_H2SO4 + 1.00E+00/R_ATM/T_K/H_HBr);
        W_HBr = M_HBr*8.091E+01/(1.00E+3 + 9.8076E+01*M_H2SO4+ 6.3012E+01*M_HNO3);
        PVAP_HBr = M_HBr/H_HBr;
        HBrGASFRAC = 1.00E+00 - (PATMHBr - PVAP_HBr)/PATMHBr;

    } else {

       W_HBr = 0.00E+00;
       HBrGASFRAC = 1.00E+00;

    }

    /* Finally HOBr (Hanson and Ravishankara, Heterogeneous
     * chemistry of Bromine species in sulfuric acid under
     * stratospheric conditions, 1995) */
    if ( PATMHOBr > 0.0E+00 ) {
          
        /* Hanson and Ravishankara state that the volume-based
         * Henry's Law coefficient for HOBr in H2SO4 is 10^6 M/atm.
         * The molality-based Henry's law constant, H_HOBr, is
         * therefore: */
        H_HOBr = 1.00E+06*MOLDENS_H2SO4/M_H2SO4;
        M_HOBr = (1.00E+00/R_ATM/T_K*PATMHOBr)/(MOLDENS_H2SO4/M_H2SO4 + 1.00E+00/R_ATM/T_K/H_HOBr);
        W_HOBr = M_HOBr*9.6911E+01/(1.00E+03 + 9.8076E+01*M_H2SO4 + 6.3012E+01*M_HNO3);
        PVAP_HOBr = M_HOBr/H_HOBr;
        HOBrGASFRAC = 1.00E+00 - (PATMHOBr - PVAP_HOBr)/PATMHOBr;

    } else {

       W_HOBr = 0.00E+00;
       HOBrGASFRAC = 1.00E+00;

    }

    /* Take W_H2O as remainder */
    W_H2O = 1.0E+00 - (W_H2SO4 + W_HNO3 + W_HCl + W_HOCl + W_HBr + W_HOBr);

    /* Aerosol mass density in kg/m3 aerosol */
    SLA_RHO = CARSLAW_DENSITY(M_H2SO4,M_HNO3,T_K);

    /* Aerosol volume in m3/m3 air */
    SLA_VOL = (MOLDENS_H2SO4*98.076E+00/W_H2SO4/SLA_RHO)*1.00E-03;


} /* End of TERNARY */

double CARSLAW_DENSITY( const double CS, const double CN, const double T_K )
{

    /* DESCRIPTION: Determines the density of a solution through a relationship
     * from Carslaw et al.. Result is in kg/m3. */

    /* INPUTS:
     * - CS  : H2SO4 molality (mol H2SO4/kg solvent)
     * - CN  : HNO3 molality (mol HNO3/kg solvent)
     * - T_K : Temperature in K */

    double DENSS, DENSN;

    /* CARSLAW_DENSITY begins here! */

    DENSS =   1.000E+03 + 1.2364E+02*CS - 5.600E-04*CS*T_K*T_K - 2.954E+01*pow( CS, 1.5 ) \
            + 1.814E-04*pow( CS, 1.5 )*T_K*T_K + 2.343E+00*CS*CS - 1.487E-03*CS*CS*T_K    \
            - 1.324E-05*CS*CS*T_K*T_K;

    DENSN =   1.000E+03 + 8.5107E+01*CN - 5.043E-04*CN*T_K*T_K - 1.896E+01*pow( CN, 1.5 ) \
            + 1.427E-04*pow( CN, 1.5 )*T_K*T_K + 1.458E+00*CN*CN - 1.198E-03*CN*CN*T_K    \
            - 9.703E-06*CN*CN*T_K*T_K;

    return 1.00E+00/( (1.00E+00/DENSS*CS/(CS+CN) + 1.00E+00/DENSN*CN/(CS+CN)) );


} /* End of CARSLAW_DENSITY */

/* End of LiquidAer.cpp */
