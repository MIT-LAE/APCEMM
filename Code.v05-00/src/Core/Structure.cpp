/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Structure Program File                                           */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Structure.cpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Core/Structure.hpp"

Solution::Solution( ) : \
        liquidAerosol( ), 
        solidAerosol( ),
        LA_Kernel( ), 
        PA_Kernel( ),
        nVariables( NSPEC ), 
        nAer( N_AER ), 
        size_x( NX ), 
        size_y( NY ),
        reducedSize( 0 )
{

    /* Constructor */

} /* End of Solution::Solution */

Solution::~Solution()
{

    /* Destructor */

} /* End of Solution::~Solution */

void Solution::Clear( Vector_2D& vector_2D )
{

    for ( UInt i = 0; i < vector_2D.size(); i++ ) {
        vector_2D[i].clear();
    }
    vector_2D.clear();

} /* End of Solution::Clear */

void Solution::SetShape( Vector_2D& vector_2D, \
                         const UInt n_x,       \
                         const UInt n_y,       \
                         const RealDouble value )
{

    Clear( vector_2D );

    /* Dimensions are transposed! */
    for ( UInt i = 0; i < n_y; i++ ) {
        vector_2D.push_back( Vector_1D( n_x, value ) );
    }

} /* End of Solution::SetShape */

void Solution::SetToValue( Vector_2D& vector_2D, \
                           const RealDouble value )
{

    for ( UInt i = 0; i < vector_2D.size(); i++ ) {
        for ( UInt j = 0; j < vector_2D[0].size(); j++ ) {
            vector_2D[i][j] = value;
        }
    }

} /* End of Solution::SetToValue */

void Solution::Print( const Vector_2D& vector_2D, \
                      const UInt i_max,           \
                      const UInt j_max ) const
{

    for ( UInt i = 0; i < i_max; i++ ) {
        for ( UInt j = 0; j < j_max; j++ ) {
            std::cout << vector_2D[i][j] << ", ";
        }
        std::cout << std::endl;
    }

} /* End of Solution::Print */

void Solution::Initialize( char const *fileName,      \
                           const Input &input,        \
                           const RealDouble airDens,  \
                           const Meteorology &met,    \
                           const OptInput &Input_Opt, \
                           const bool DBG )
{

    Vector_1D amb_Value(nVariables, 0.0);
    Vector_2D aer_Value(nAer, Vector_1D(2, 0.0));
    std::ifstream file;


    /* Read input background conditions */
    file.open( fileName );

    if ( file.is_open() ) {
        std::string line;
        UInt i = 0;

        while ( ( std::getline( file, line ) ) && ( i < nVariables + nAer ) ) {
            if ( ( line.length() > 0 ) && ( line != "\r" ) && ( line != "\n" ) && ( line[0] != '#' ) ) {
                std::istringstream iss(line);
                if ( i < nVariables ) {
                    iss >> amb_Value[i];
                }
                else if ( ( i >= nVariables ) && ( i < nVariables + nAer ) ) {
                    iss >> aer_Value[i - nVariables][0];
                    std::getline( file, line );
                    std::istringstream iss(line);
                    iss >> aer_Value[i - nVariables][1];
                }
                i++;
            }
        }
        file.close();
    }
    else {
        std::string const currFunc("Structure::Initialize");
        std::cout << "ERROR: In " << currFunc << ": Can't read (" << fileName << ")" << std::endl;
        exit(-1);
    }

    SpinUp( amb_Value, input, airDens, \
            /* Time for which ambient file is valid in hr */ (const RealDouble) 8.0 );


    /* Enforce pre-defined values? *
     * Read input defined values for background concentrations */

    /* Split NOx as NO and NO2 using the current NO/NO2 ratio */

    /* Inputs are in ppb */

    if ( input.backgNOx() > 0.0E+00 ) {
        const RealDouble NONO2rat = amb_Value[ind_NO]/amb_Value[ind_NO2];
        /* NOx = NO + NO2 = NO2 * ( r + 1 )
         * NO2 = NOx / ( r + 1 );
         * NO  = NOx - NO2; */
        amb_Value[ind_NO2] = input.backgNOx() / ( NONO2rat + 1 );
        amb_Value[ind_NO]  = input.backgNOx() - amb_Value[ind_NO2];

        /* Convert to mixing ratio */
        amb_Value[ind_NO2] /= 1.0E+09;
        amb_Value[ind_NO]  /= 1.0E+09;
    }

    if ( input.backgHNO3() > 0.0E+00 ) {
        /* Convert from ppb to mixing ratio */
        amb_Value[ind_HNO3] = input.backgHNO3() / 1.0E+09;
    }

    if ( input.backgO3() > 0.0E+00 ) {
        /* Convert from ppb to mixing ratio */
        amb_Value[ind_O3] = input.backgO3() / 1.0E+09;
    }

    if ( input.backgCO() > 0.0E+00 ) {
        /* Convert from ppb to mixing ratio */
        amb_Value[ind_CO] = input.backgCO() / 1.0E+09;
    }

    if ( input.backgCH4() > 0.0E+00 ) {
        /* Convert from ppb to mixing ratio */
        amb_Value[ind_CH4] = input.backgCH4() / 1.0E+09;
    }

    if ( input.backgSO2() > 0.0E+00 ) {
        /* Convert from ppb to mixing ratio */
        amb_Value[ind_SO2] = input.backgSO2() / 1.0E+09;
    }


    /* Initialize and allocate space for species */

    UInt actualX = size_x;
    UInt actualY = size_y;

    if ( !Input_Opt.CHEMISTRY_CHEMISTRY ) {
        actualX = 1;
        actualY = 1;

        reducedSize = 1;
    }

    /* Gaseous species */
    SetShape( H2O      , size_x  , size_y  , amb_Value[112] * airDens );
    SetShape( H2O_met  , size_x  , size_y  , amb_Value[112] * airDens );
    SetShape( H2O_plume, size_x  , size_y  , 0.0E+00 );


    SetShape( CO2      , actualX , actualY , amb_Value[  0] * airDens );
    SetShape( PPN      , actualX , actualY , amb_Value[  1] * airDens );
    SetShape( BrNO2    , actualX , actualY , amb_Value[  2] * airDens );
    SetShape( IEPOX    , actualX , actualY , amb_Value[  3] * airDens );
    SetShape( PMNN     , actualX , actualY , amb_Value[  4] * airDens );
    SetShape( N2O      , actualX , actualY , amb_Value[  5] * airDens );
    SetShape( N        , actualX , actualY , amb_Value[  6] * airDens );
    SetShape( PAN      , actualX , actualY , amb_Value[  7] * airDens );
    SetShape( ALK4     , actualX , actualY , amb_Value[  8] * airDens );
    SetShape( MAP      , actualX , actualY , amb_Value[  9] * airDens );
    SetShape( MPN      , actualX , actualY , amb_Value[ 10] * airDens );
    SetShape( Cl2O2    , actualX , actualY , amb_Value[ 11] * airDens );
    SetShape( ETP      , actualX , actualY , amb_Value[ 12] * airDens );
    SetShape( HNO2     , actualX , actualY , amb_Value[ 13] * airDens );
    SetShape( C3H8     , actualX , actualY , amb_Value[ 14] * airDens );
    SetShape( RA3P     , actualX , actualY , amb_Value[ 15] * airDens );
    SetShape( RB3P     , actualX , actualY , amb_Value[ 16] * airDens );
    SetShape( OClO     , actualX , actualY , amb_Value[ 17] * airDens );
    SetShape( ClNO2    , actualX , actualY , amb_Value[ 18] * airDens );
    SetShape( ISOP     , actualX , actualY , amb_Value[ 19] * airDens );
    SetShape( HNO4     , actualX , actualY , amb_Value[ 20] * airDens );
    SetShape( MAOP     , actualX , actualY , amb_Value[ 21] * airDens );
    SetShape( MP       , actualX , actualY , amb_Value[ 22] * airDens );
    SetShape( ClOO     , actualX , actualY , amb_Value[ 23] * airDens );
    SetShape( RP       , actualX , actualY , amb_Value[ 24] * airDens );
    SetShape( BrCl     , actualX , actualY , amb_Value[ 25] * airDens );
    SetShape( PP       , actualX , actualY , amb_Value[ 26] * airDens );
    SetShape( PRPN     , actualX , actualY , amb_Value[ 27] * airDens );
    SetShape( SO4      , actualX , actualY , amb_Value[ 28] * airDens );
    SetShape( Br2      , actualX , actualY , amb_Value[ 29] * airDens );
    SetShape( ETHLN    , actualX , actualY , amb_Value[ 30] * airDens );
    SetShape( MVKN     , actualX , actualY , amb_Value[ 31] * airDens );
    SetShape( R4P      , actualX , actualY , amb_Value[ 32] * airDens );
    SetShape( C2H6     , actualX , actualY , amb_Value[ 33] * airDens );
    SetShape( RIP      , actualX , actualY , amb_Value[ 34] * airDens );
    SetShape( VRP      , actualX , actualY , amb_Value[ 35] * airDens );
    SetShape( ATOOH    , actualX , actualY , amb_Value[ 36] * airDens );
    SetShape( IAP      , actualX , actualY , amb_Value[ 37] * airDens );
    SetShape( DHMOB    , actualX , actualY , amb_Value[ 38] * airDens );
    SetShape( MOBA     , actualX , actualY , amb_Value[ 39] * airDens );
    SetShape( MRP      , actualX , actualY , amb_Value[ 40] * airDens );
    SetShape( N2O5     , actualX , actualY , amb_Value[ 41] * airDens );
    SetShape( ISNOHOO  , actualX , actualY , amb_Value[ 42] * airDens );
    SetShape( ISNP     , actualX , actualY , amb_Value[ 43] * airDens );
    SetShape( ISOPNB   , actualX , actualY , amb_Value[ 44] * airDens );
    SetShape( IEPOXOO  , actualX , actualY , amb_Value[ 45] * airDens );
    SetShape( MACRNO2  , actualX , actualY , amb_Value[ 46] * airDens );
    SetShape( ROH      , actualX , actualY , amb_Value[ 47] * airDens );
    SetShape( MOBAOO   , actualX , actualY , amb_Value[ 48] * airDens );
    SetShape( DIBOO    , actualX , actualY , amb_Value[ 49] * airDens );
    SetShape( PMN      , actualX , actualY , amb_Value[ 50] * airDens );
    SetShape( ISNOOB   , actualX , actualY , amb_Value[ 51] * airDens );
    SetShape( INPN     , actualX , actualY , amb_Value[ 52] * airDens );
    SetShape( H        , actualX , actualY , amb_Value[ 53] * airDens );
    SetShape( BrNO3    , actualX , actualY , amb_Value[ 54] * airDens );
    SetShape( PRPE     , actualX , actualY , amb_Value[ 55] * airDens );
    SetShape( MVKOO    , actualX , actualY , amb_Value[ 56] * airDens );
    SetShape( Cl2      , actualX , actualY , amb_Value[ 57] * airDens );
    SetShape( ISOPND   , actualX , actualY , amb_Value[ 58] * airDens );
    SetShape( HOBr     , actualX , actualY , amb_Value[ 59] * airDens );
    SetShape( A3O2     , actualX , actualY , amb_Value[ 60] * airDens );
    SetShape( PROPNN   , actualX , actualY , amb_Value[ 61] * airDens );
    SetShape( GLYX     , actualX , actualY , amb_Value[ 62] * airDens );
    SetShape( MAOPO2   , actualX , actualY , amb_Value[ 63] * airDens );
    SetShape( CH4      , actualX , actualY , amb_Value[ 64] * airDens );
    SetShape( GAOO     , actualX , actualY , amb_Value[ 65] * airDens );
    SetShape( B3O2     , actualX , actualY , amb_Value[ 66] * airDens );
    SetShape( ACET     , actualX , actualY , amb_Value[ 67] * airDens );
    SetShape( MACRN    , actualX , actualY , amb_Value[ 68] * airDens );
    SetShape( CH2OO    , actualX , actualY , amb_Value[ 69] * airDens );
    SetShape( MGLYOO   , actualX , actualY , amb_Value[ 70] * airDens );
    SetShape( VRO2     , actualX , actualY , amb_Value[ 71] * airDens );
    SetShape( MGLOO    , actualX , actualY , amb_Value[ 72] * airDens );
    SetShape( MACROO   , actualX , actualY , amb_Value[ 73] * airDens );
    SetShape( PO2      , actualX , actualY , amb_Value[ 74] * airDens );
    SetShape( CH3CHOO  , actualX , actualY , amb_Value[ 75] * airDens );
    SetShape( MAN2     , actualX , actualY , amb_Value[ 76] * airDens );
    SetShape( ISNOOA   , actualX , actualY , amb_Value[ 77] * airDens );
    SetShape( H2O2     , actualX , actualY , amb_Value[ 78] * airDens );
    SetShape( PRN1     , actualX , actualY , amb_Value[ 79] * airDens );
    SetShape( ETO2     , actualX , actualY , amb_Value[ 80] * airDens );
    SetShape( KO2      , actualX , actualY , amb_Value[ 81] * airDens );
    SetShape( RCO3     , actualX , actualY , amb_Value[ 82] * airDens );
    SetShape( HC5OO    , actualX , actualY , amb_Value[ 83] * airDens );
    SetShape( GLYC     , actualX , actualY , amb_Value[ 84] * airDens );
    SetShape( ClNO3    , actualX , actualY , amb_Value[ 85] * airDens );
    SetShape( RIO2     , actualX , actualY , amb_Value[ 86] * airDens );
    SetShape( R4N1     , actualX , actualY , amb_Value[ 87] * airDens );
    SetShape( HOCl     , actualX , actualY , amb_Value[ 88] * airDens );
    SetShape( ATO2     , actualX , actualY , amb_Value[ 89] * airDens );
    SetShape( HNO3     , actualX , actualY , amb_Value[ 90] * airDens );
    SetShape( ISN1     , actualX , actualY , amb_Value[ 91] * airDens );
    SetShape( MAO3     , actualX , actualY , amb_Value[ 92] * airDens );
    SetShape( MRO2     , actualX , actualY , amb_Value[ 93] * airDens );
    SetShape( INO2     , actualX , actualY , amb_Value[ 94] * airDens );
    SetShape( HAC      , actualX , actualY , amb_Value[ 95] * airDens );
    SetShape( HC5      , actualX , actualY , amb_Value[ 96] * airDens );
    SetShape( MGLY     , actualX , actualY , amb_Value[ 97] * airDens );
    SetShape( ISOPNBO2 , actualX , actualY , amb_Value[ 98] * airDens );
    SetShape( ISOPNDO2 , actualX , actualY , amb_Value[ 99] * airDens );
    SetShape( R4O2     , actualX , actualY , amb_Value[100] * airDens );
    SetShape( R4N2     , actualX , actualY , amb_Value[101] * airDens );
    SetShape( BrO      , actualX , actualY , amb_Value[102] * airDens );
    SetShape( RCHO     , actualX , actualY , amb_Value[103] * airDens );
    SetShape( MEK      , actualX , actualY , amb_Value[104] * airDens );
    SetShape( ClO      , actualX , actualY , amb_Value[105] * airDens );
    SetShape( MACR     , actualX , actualY , amb_Value[106] * airDens );
    SetShape( SO2      , actualX , actualY , amb_Value[107] * airDens );
    SetShape( MVK      , actualX , actualY , amb_Value[108] * airDens );
    SetShape( ALD2     , actualX , actualY , amb_Value[109] * airDens );
    SetShape( MCO3     , actualX , actualY , amb_Value[110] * airDens );
    SetShape( CH2O     , actualX , actualY , amb_Value[111] * airDens );
    SetShape( Br       , actualX , actualY , amb_Value[113] * airDens );
    SetShape( NO       , actualX , actualY , amb_Value[114] * airDens );
    SetShape( NO3      , actualX , actualY , amb_Value[115] * airDens );
    SetShape( Cl       , actualX , actualY , amb_Value[116] * airDens );
    SetShape( O        , actualX , actualY , amb_Value[117] * airDens );
    SetShape( O1D      , actualX , actualY , amb_Value[118] * airDens );
    SetShape( O3       , actualX , actualY , amb_Value[119] * airDens );
    SetShape( HO2      , actualX , actualY , amb_Value[120] * airDens );
    SetShape( NO2      , actualX , actualY , amb_Value[121] * airDens );
    SetShape( OH       , actualX , actualY , amb_Value[122] * airDens );
    SetShape( HBr      , actualX , actualY , amb_Value[123] * airDens );
    SetShape( HCl      , actualX , actualY , amb_Value[124] * airDens );
    SetShape( CO       , actualX , actualY , amb_Value[125] * airDens );
    SetShape( MO2      , actualX , actualY , amb_Value[126] * airDens );
    SetShape( ACTA     , actualX , actualY , amb_Value[127] * airDens );
    SetShape( EOH      , actualX , actualY , amb_Value[128] * airDens );
    SetShape( H2       , actualX , actualY , amb_Value[129] * airDens );
    SetShape( HCOOH    , actualX , actualY , amb_Value[130] * airDens );
    SetShape( MOH      , actualX , actualY , amb_Value[131] * airDens );
    SetShape( N2       , actualX , actualY , amb_Value[132] * airDens );
    SetShape( O2       , actualX , actualY , amb_Value[133] * airDens );
    SetShape( RCOOH    , actualX , actualY , amb_Value[134] * airDens );

    /* Solid-liquid species */
    SetShape( NIT      , actualX , actualY , (RealDouble) 0.0E+00     );
    SetShape( NAT      , actualX , actualY , (RealDouble) 0.0E+00     );
    SetShape( SO4L     , actualX , actualY , (RealDouble) 0.0E+00     );
    SetShape( H2OL     , size_x  , size_y  , (RealDouble) 0.0E+00     );
    SetShape( H2OS     , size_x  , size_y  , (RealDouble) 0.0E+00     );
    SetShape( HNO3L    , actualX , actualY , (RealDouble) 0.0E+00     );
    SetShape( HNO3S    , actualX , actualY , (RealDouble) 0.0E+00     );
    SetShape( HClL     , actualX , actualY , (RealDouble) 0.0E+00     );
    SetShape( HOClL    , actualX , actualY , (RealDouble) 0.0E+00     );
    SetShape( HBrL     , actualX , actualY , (RealDouble) 0.0E+00     );
    SetShape( HOBrL    , actualX , actualY , (RealDouble) 0.0E+00     );

    /* Tracers */
    SetShape( SO4T     , actualX , actualY , amb_Value[ 28] * airDens );


    if ( Input_Opt.MET_LOADMET ) {
        /* Use meteorological input? */
        //H2O = met.H2O_;
        H2O_met = met.H2O_;
        /* Update H2O */
        for ( UInt i = 0; i < size_x; i++ ) {
            for ( UInt j = 0; j < size_y; j++ ) {
                H2O[j][i] = H2O_met[j][i] + H2O_plume[j][i];
            }
        }
    } else {
        /* Else use user-defined H2O profile */
        RealDouble H2Oval = (input.relHumidity_w()/((RealDouble) 100.0) * \
                          physFunc::pSat_H2Ol( input.temperature_K() ) / ( physConst::kB * input.temperature_K() )) / 1.00E+06;
        for ( UInt i = 0; i < size_x; i++ ) {
            for ( UInt j = 0; j < size_y; j++ ) {
                //H2O[j][i] = H2Oval;
                H2O_met[j][i] = H2Oval;
                /* RH_w = x_H2O * P / Psat_H2Ol(T) = [H2O](#/cm3) * 1E6 * kB * T / Psat_H2Ol(T) */
                H2O[j][i] = H2O_met[j][i] + H2O_plume[j][i];
            }
        }
    }

    Vector_1D stratData{ SO4[0][0], HNO3[0][0], HCl[0][0], HOCl[0][0],  \
                         HBr[0][0], HOBr[0][0], H2O[0][0], ClNO3[0][0], \
                         BrNO3[0][0], NIT[0][0], NAT[0][0] };

    KHETI_SLA.assign( 11, 0.0 );
    AERFRAC.assign( 7, 0.0 );
    SOLIDFRAC.assign( 7, 0.0 );

    Vector_1D RAD, RHO, KG, NDENS, SAD;
    RAD.assign( 2, 0.0 );
    RHO.assign( 2, 0.0 );
    KG.assign( 2, 0.0 );
    NDENS.assign( 2, 0.0 );
    SAD.assign( 2, 0.0 );

    STATE_PSC = STRAT_AER( input.temperature_K(), input.pressure_Pa(), airDens,  \
                           input.latitude_deg(), stratData,                      \
                           (2.0*XLIM)*(YLIM_UP+YLIM_DOWN), KHETI_SLA, SOLIDFRAC, \
                           AERFRAC, RAD, RHO, KG, NDENS, SAD, DBG );

    /* Liquid/solid species */
    SetToValue( SO4L, (RealDouble) AERFRAC[0]                          * stratData[0] );
    SetToValue( SO4 , (RealDouble) ( 1.0 - AERFRAC[0] )                * stratData[0] );

    AERFRAC[6] = 0.0E+00;
    SOLIDFRAC[6] = 0.0E+00;
    SetToValue( H2OL , (RealDouble) AERFRAC[6]                          * stratData[6] );
    SetToValue( H2OS , (RealDouble) SOLIDFRAC[6]                        * stratData[6] );
    /* Do not overwrite H2O!! */
//    SetToValue( H2O  , (RealDouble) ( 1.0 - AERFRAC[6] - SOLIDFRAC[6] ) * stratData[6] );

    SetToValue( HNO3L, (RealDouble) AERFRAC[1]                          * stratData[1] );
    SetToValue( HNO3S, (RealDouble) SOLIDFRAC[1]                        * stratData[1] );
    SetToValue( HNO3 , (RealDouble) ( 1.0 - AERFRAC[1] - SOLIDFRAC[1] ) * stratData[1] );

    SetToValue( HClL , (RealDouble) AERFRAC[2]                          * stratData[2] );
    SetToValue( HCl  , (RealDouble) ( 1.0 - AERFRAC[2] )                * stratData[2] );

    SetToValue( HOClL, (RealDouble) AERFRAC[3]                          * stratData[3] );
    SetToValue( HOCl , (RealDouble) ( 1.0 - AERFRAC[3] )                * stratData[3] );

    SetToValue( HBrL , (RealDouble) AERFRAC[4]                          * stratData[4] );
    SetToValue( HBr  , (RealDouble) ( 1.0 - AERFRAC[4] )                * stratData[4] );

    SetToValue( HOBrL, (RealDouble) AERFRAC[5]                          * stratData[5] );
    SetToValue( HOBr , (RealDouble) ( 1.0 - AERFRAC[5] )                * stratData[5] );

    SetToValue( NIT  , (RealDouble) stratData[ 9] );
    SetToValue( NAT  , (RealDouble) stratData[10] );


    /* Aerosols */
    /* Assume that soot particles are monodisperse */
    SetShape( sootDens , size_x, size_y, (RealDouble) aer_Value[  0][0] );
    SetShape( sootRadi , size_x, size_y, (RealDouble) aer_Value[  0][1] );
    SetShape( sootArea , size_x, size_y, (RealDouble) 4.0 / RealDouble(3.0) * physConst::PI * aer_Value[  0][0] * aer_Value[  0][1] * aer_Value[  0][1] * aer_Value[  0][1] );


    nBin_LA = std::floor( 1 + log( pow( (LA_R_HIG/LA_R_LOW), 3.0 ) ) / log( LA_VRAT ) );

    if ( DBG ) {
        std::cout << "\n DEBUG : LA_R_LOW  = " << LA_R_LOW * 1.00E+09 << " [nm]\n";
        std::cout << " DEBUG : LA_R_HIG  = "   << LA_R_HIG * 1.00E+09 << " [nm]\n";
        std::cout << " DEBUG : LA_VRAT   = "   << LA_VRAT             << " [-]\n";
        std::cout << " DEBUG : nBin_LA   = "   << nBin_LA             << "\n";
        std::cout << " DEBUG : NDENS     = "   << NDENS[1] * 1.00E-06 << " [#/cm^3]\n";
        std::cout << " DEBUG : REFF      = "   << RAD[1] * 1.00E+06   << " [mum]\n";
    }

    Vector_1D LA_rE( nBin_LA + 1, 0.0 ); /* Bin edges in m */
    Vector_1D LA_rJ( nBin_LA    , 0.0 ); /* Bin center radius in m */
    Vector_1D LA_vJ( nBin_LA    , 0.0 ); /* Bin volume centers in m^3 */

    const RealDouble LA_RRAT = pow( LA_VRAT, 1.0 / RealDouble(3.0) );
    LA_rE[0] = LA_R_LOW;
    for ( UInt iBin_LA = 1; iBin_LA < nBin_LA + 1; iBin_LA++ )                              /* [m] */
        LA_rE[iBin_LA] = LA_rE[iBin_LA-1] * LA_RRAT;

    for ( UInt iBin_LA = 0; iBin_LA < nBin_LA; iBin_LA++ ) {
        LA_rJ[iBin_LA] = 0.5 * ( LA_rE[iBin_LA] + LA_rE[iBin_LA+1] );                       /* [m] */
        LA_vJ[iBin_LA] = 4.0 / RealDouble(3.0) * physConst::PI * \
                         ( LA_rE[iBin_LA] * LA_rE[iBin_LA] * LA_rE[iBin_LA] \
                         + LA_rE[iBin_LA+1] * LA_rE[iBin_LA+1] * LA_rE[iBin_LA+1] ) * 0.5;  /* [m^3] */
    }

    LA_nDens = NDENS[1] * 1.00E-06; /* [#/cm^3]      */
    LA_rEff  = RAD[1]   * 1.00E+09; /* [nm]          */
    LA_SAD   = SAD[1]   * 1.00E+06; /* [\mum^2/cm^3] */

    if ( LA_nDens >= 0.0E+00 ) {
        /* For a lognormal distribution:
         * r_eff = r_m * exp( 5/2 * ln(S)^2 )
         * A     = 4\pi N0 r_m^2 * exp ( 2 * ln(S)^2 )
         * A/r_eff^2 = 4\pi N0 * exp( - 3 * ln(S)^2 )
         *
         * ln(S) = sqrt(-1/3*ln(A/(4\pi r_eff^2 * N0)));
         * r_m = r_eff * exp( -5/2 * ln(S)^2 ); */

        const RealDouble sLA = sqrt( - 1.0 / (3.0) * log(SAD[1]/(4.0 * physConst::PI * RAD[1] * RAD[1] * NDENS[1] ) ) );
        const RealDouble rLA = std::max( RAD[1] * exp( - 2.5 * sLA * sLA ), 1.5 * LA_R_LOW );
        const AIM::Grid_Aerosol LAAerosol( size_x, size_y, LA_rJ, LA_rE, LA_nDens, rLA, exp(sLA), "lognormal" );

        liquidAerosol = LAAerosol;
    }

    const AIM::Coagulation kernel1( "liquid", LA_rJ, LA_vJ, physConst::RHO_SULF, \
                                    input.temperature_K(), input.pressure_Pa() );

    LA_Kernel = kernel1;

    if ( DBG ) {

        std::cout << "\n DEBUG : Comparing PDF's number density to exact number density :\n";
        std::cout << "         " << liquidAerosol.Moment( 0, 0, 0 ) << " v " << LA_nDens << " [#/cm^3]\n";
        std::cout << " DEBUG : Comparing PDF's surface area to Grainger's surface area:\n";
        std::cout << "         " << 4.0 * physConst::PI * liquidAerosol.Moment( 2, 0, 0 ) * 1.00E+12 << " v " << LA_SAD << " [mum^2/cm^3]\n";
        std::cout << " DEBUG : Comparing PDF's effective radius to Grainger's effective radius:\n";
        std::cout << "         " << liquidAerosol.EffRadius( 0, 0 ) * 1.00E+09 << " v " << LA_rEff << " [nm]\n";

    }

    nBin_PA = std::floor( 1 + log( pow( (PA_R_HIG/PA_R_LOW), 3.0 ) ) / log( PA_VRAT ) );

    if ( DBG ) {
        std::cout << "\n DEBUG : PA_R_LOW  = " << PA_R_LOW * 1.00E+06 << " [mum]\n";
        std::cout << " DEBUG : PA_R_HIG  = "   << PA_R_HIG * 1.00E+06 << " [mum]\n";
        std::cout << " DEBUG : PA_VRAT   = "   << PA_VRAT             << " [-]\n";
        std::cout << " DEBUG : nBin_PA   = "   << nBin_PA             << "\n";
        std::cout << " DEBUG : NDENS     = "   << NDENS[0] * 1.00E-06 << " [#/cm^3]\n";
        std::cout << " DEBUG : REFF      = "   << RAD[0] * 1.00E+06   << " [mum]\n";
    }

    Vector_1D PA_rE( nBin_PA + 1, 0.0 ); /* Bin edges in m */
    Vector_1D PA_rJ( nBin_PA    , 0.0 ); /* Bin center radius in m */
    Vector_1D PA_vJ( nBin_PA    , 0.0 ); /* Bin volume centers in m^3 */

    const RealDouble PA_RRAT = pow( PA_VRAT, 1.0 / RealDouble(3.0) );
    PA_rE[0] = PA_R_LOW;
    for ( UInt iBin_PA = 1; iBin_PA < nBin_PA + 1; iBin_PA++ )
        PA_rE[iBin_PA] = PA_rE[iBin_PA-1] * PA_RRAT;                                        /* [m]   */

    for ( UInt iBin_PA = 0; iBin_PA < nBin_PA; iBin_PA++ ) {
        PA_rJ[iBin_PA] = 0.5 * ( PA_rE[iBin_PA] + PA_rE[iBin_PA+1] );                       /* [m]   */
        PA_vJ[iBin_PA] = 4.0 / RealDouble(3.0) * physConst::PI * \
                         ( PA_rE[iBin_PA] * PA_rE[iBin_PA] * PA_rE[iBin_PA] \
                         + PA_rE[iBin_PA+1] * PA_rE[iBin_PA+1] * PA_rE[iBin_PA+1] ) * 0.5;  /* [m^3] */
    }

    PA_nDens = NDENS[0] * 1.00E-06; /* [#/cm^3]      */
    PA_rEff  = RAD[0]   * 1.00E+09; /* [nm]          */
    PA_SAD   = SAD[0]   * 1.00E+06; /* [\mum^2/cm^3] */

    if ( PA_nDens >= 0.0E+00 ) {
        const RealDouble expsPA = 1.15;
        const RealDouble rPA = std::max( RAD[0] * exp( - 2.5 * log(expsPA) * log(expsPA) ), 1.5 * PA_R_LOW );
        AIM::Grid_Aerosol PAAerosol( size_x, size_y, PA_rJ, PA_rE, PA_nDens, rPA, expsPA, "lognormal" );

        solidAerosol = PAAerosol;
    }

    const AIM::Coagulation kernel2( "ice", PA_rJ, PA_vJ, physConst::RHO_ICE, \
                                    input.temperature_K(), input.pressure_Pa() );

    PA_Kernel = kernel2;

    if ( DBG ) {
        std::cout << "\n DEBUG : Comparing PDF's number density to exact number density :\n";
        std::cout << "         " << solidAerosol.Moment( 0, 0, 0 ) << " v " << PA_nDens << " [#/cm^3]\n";
        std::cout << " DEBUG : Comparing PDF's effective radius to actual effective radius:\n";
        std::cout << "         " << solidAerosol.EffRadius( 0, 0 ) * 1.00E+09 << " v " << PA_rEff << " [nm]\n";
    }

} /* End of Solution::Initialize */

void Solution::getData( const UInt i, \
                        const UInt j )
{

    VAR[  0] = CO2[j][i];
    VAR[  1] = PPN[j][i];
    VAR[  2] = BrNO2[j][i];
    VAR[  3] = IEPOX[j][i];
    VAR[  4] = PMNN[j][i];
    VAR[  5] = N2O[j][i];
    VAR[  6] = N[j][i];
    VAR[  7] = PAN[j][i];
    VAR[  8] = ALK4[j][i];
    VAR[  9] = MAP[j][i];
    VAR[ 10] = MPN[j][i];
    VAR[ 11] = Cl2O2[j][i];
    VAR[ 12] = ETP[j][i];
    VAR[ 13] = HNO2[j][i];
    VAR[ 14] = C3H8[j][i];
    VAR[ 15] = RA3P[j][i];
    VAR[ 16] = RB3P[j][i];
    VAR[ 17] = OClO[j][i];
    VAR[ 18] = ClNO2[j][i];
    VAR[ 19] = ISOP[j][i];
    VAR[ 20] = HNO4[j][i];
    VAR[ 21] = MAOP[j][i];
    VAR[ 22] = MP[j][i];
    VAR[ 23] = ClOO[j][i];
    VAR[ 24] = RP[j][i];
    VAR[ 25] = BrCl[j][i];
    VAR[ 26] = PP[j][i];
    VAR[ 27] = PRPN[j][i];
    VAR[ 28] = SO4[j][i];
    VAR[ 29] = Br2[j][i];
    VAR[ 30] = ETHLN[j][i];
    VAR[ 31] = MVKN[j][i];
    VAR[ 32] = R4P[j][i];
    VAR[ 33] = C2H6[j][i];
    VAR[ 34] = RIP[j][i];
    VAR[ 35] = VRP[j][i];
    VAR[ 36] = ATOOH[j][i];
    VAR[ 37] = IAP[j][i];
    VAR[ 38] = DHMOB[j][i];
    VAR[ 39] = MOBA[j][i];
    VAR[ 40] = MRP[j][i];
    VAR[ 41] = N2O5[j][i];
    VAR[ 42] = ISNOHOO[j][i];
    VAR[ 43] = ISNP[j][i];
    VAR[ 44] = ISOPNB[j][i];
    VAR[ 45] = IEPOXOO[j][i];
    VAR[ 46] = MACRNO2[j][i];
    VAR[ 47] = ROH[j][i];
    VAR[ 48] = MOBAOO[j][i];
    VAR[ 49] = DIBOO[j][i];
    VAR[ 50] = PMN[j][i];
    VAR[ 51] = ISNOOB[j][i];
    VAR[ 52] = INPN[j][i];
    VAR[ 53] = H[j][i];
    VAR[ 54] = BrNO3[j][i];
    VAR[ 55] = PRPE[j][i];
    VAR[ 56] = MVKOO[j][i];
    VAR[ 57] = Cl2[j][i];
    VAR[ 58] = ISOPND[j][i];
    VAR[ 59] = HOBr[j][i];
    VAR[ 60] = A3O2[j][i];
    VAR[ 61] = PROPNN[j][i];
    VAR[ 62] = GLYX[j][i];
    VAR[ 63] = MAOPO2[j][i];
    VAR[ 64] = CH4[j][i];
    VAR[ 65] = GAOO[j][i];
    VAR[ 66] = B3O2[j][i];
    VAR[ 67] = ACET[j][i];
    VAR[ 68] = MACRN[j][i];
    VAR[ 69] = CH2OO[j][i];
    VAR[ 70] = MGLYOO[j][i];
    VAR[ 71] = VRO2[j][i];
    VAR[ 72] = MGLOO[j][i];
    VAR[ 73] = MACROO[j][i];
    VAR[ 74] = PO2[j][i];
    VAR[ 75] = CH3CHOO[j][i];
    VAR[ 76] = MAN2[j][i];
    VAR[ 77] = ISNOOA[j][i];
    VAR[ 78] = H2O2[j][i];
    VAR[ 79] = PRN1[j][i];
    VAR[ 80] = ETO2[j][i];
    VAR[ 81] = KO2[j][i];
    VAR[ 82] = RCO3[j][i];
    VAR[ 83] = HC5OO[j][i];
    VAR[ 84] = GLYC[j][i];
    VAR[ 85] = ClNO3[j][i];
    VAR[ 86] = RIO2[j][i];
    VAR[ 87] = R4N1[j][i];
    VAR[ 88] = HOCl[j][i];
    VAR[ 89] = ATO2[j][i];
    VAR[ 90] = HNO3[j][i];
    VAR[ 91] = ISN1[j][i];
    VAR[ 92] = MAO3[j][i];
    VAR[ 93] = MRO2[j][i];
    VAR[ 94] = INO2[j][i];
    VAR[ 95] = HAC[j][i];
    VAR[ 96] = HC5[j][i];
    VAR[ 97] = MGLY[j][i];
    VAR[ 98] = ISOPNBO2[j][i];
    VAR[ 99] = ISOPNDO2[j][i];
    VAR[100] = R4O2[j][i];
    VAR[101] = R4N2[j][i];
    VAR[102] = BrO[j][i];
    VAR[103] = RCHO[j][i];
    VAR[104] = MEK[j][i];
    VAR[105] = ClO[j][i];
    VAR[106] = MACR[j][i];
    VAR[107] = SO2[j][i];
    VAR[108] = MVK[j][i];
    VAR[109] = ALD2[j][i];
    VAR[110] = MCO3[j][i];
    VAR[111] = CH2O[j][i];
    VAR[112] = H2O[j][i];
    VAR[113] = Br[j][i];
    VAR[114] = NO[j][i];
    VAR[115] = NO3[j][i];
    VAR[116] = Cl[j][i];
    VAR[117] = O[j][i];
    VAR[118] = O1D[j][i];
    VAR[119] = O3[j][i];
    VAR[120] = HO2[j][i];
    VAR[121] = NO2[j][i];
    VAR[122] = OH[j][i];
    VAR[123] = HBr[j][i];
    VAR[124] = HCl[j][i];
    VAR[125] = CO[j][i];
    VAR[126] = MO2[j][i];
    FIX[  0] = ACTA[j][i];
    FIX[  1] = EOH[j][i];
    FIX[  2] = H2[j][i];
    FIX[  3] = HCOOH[j][i];
    FIX[  4] = MOH[j][i];
    FIX[  5] = N2[j][i];
    FIX[  6] = O2[j][i];
    FIX[  7] = RCOOH[j][i];

} /* End of Solution::getData */

void Solution::applyData( const UInt i, \
                          const UInt j )
{

    CO2[j][i]      = VAR[  0];
    PPN[j][i]      = VAR[  1];
    BrNO2[j][i]    = VAR[  2];
    IEPOX[j][i]    = VAR[  3];
    PMNN[j][i]     = VAR[  4];
    N2O[j][i]      = VAR[  5];
    N[j][i]        = VAR[  6];
    PAN[j][i]      = VAR[  7];
    ALK4[j][i]     = VAR[  8];
    MAP[j][i]      = VAR[  9];
    MPN[j][i]      = VAR[ 10];
    Cl2O2[j][i]    = VAR[ 11];
    ETP[j][i]      = VAR[ 12];
    HNO2[j][i]     = VAR[ 13];
    C3H8[j][i]     = VAR[ 14];
    RA3P[j][i]     = VAR[ 15];
    RB3P[j][i]     = VAR[ 16];
    OClO[j][i]     = VAR[ 17];
    ClNO2[j][i]    = VAR[ 18];
    ISOP[j][i]     = VAR[ 19];
    HNO4[j][i]     = VAR[ 20];
    MAOP[j][i]     = VAR[ 21];
    MP[j][i]       = VAR[ 22];
    ClOO[j][i]     = VAR[ 23];
    RP[j][i]       = VAR[ 24];
    BrCl[j][i]     = VAR[ 25];
    PP[j][i]       = VAR[ 26];
    PRPN[j][i]     = VAR[ 27];
    SO4[j][i]      = VAR[ 28];
    Br2[j][i]      = VAR[ 29];
    ETHLN[j][i]    = VAR[ 30];
    MVKN[j][i]     = VAR[ 31];
    R4P[j][i]      = VAR[ 32];
    C2H6[j][i]     = VAR[ 33];
    RIP[j][i]      = VAR[ 34];
    VRP[j][i]      = VAR[ 35];
    ATOOH[j][i]    = VAR[ 36];
    IAP[j][i]      = VAR[ 37];
    DHMOB[j][i]    = VAR[ 38];
    MOBA[j][i]     = VAR[ 39];
    MRP[j][i]      = VAR[ 40];
    N2O5[j][i]     = VAR[ 41];
    ISNOHOO[j][i]  = VAR[ 42];
    ISNP[j][i]     = VAR[ 43];
    ISOPNB[j][i]   = VAR[ 44];
    IEPOXOO[j][i]  = VAR[ 45];
    MACRNO2[j][i]  = VAR[ 46];
    ROH[j][i]      = VAR[ 47];
    MOBAOO[j][i]   = VAR[ 48];
    DIBOO[j][i]    = VAR[ 49];
    PMN[j][i]      = VAR[ 50];
    ISNOOB[j][i]   = VAR[ 51];
    INPN[j][i]     = VAR[ 52];
    H[j][i]        = VAR[ 53];
    BrNO3[j][i]    = VAR[ 54];
    PRPE[j][i]     = VAR[ 55];
    MVKOO[j][i]    = VAR[ 56];
    Cl2[j][i]      = VAR[ 57];
    ISOPND[j][i]   = VAR[ 58];
    HOBr[j][i]     = VAR[ 59];
    A3O2[j][i]     = VAR[ 60];
    PROPNN[j][i]   = VAR[ 61];
    GLYX[j][i]     = VAR[ 62];
    MAOPO2[j][i]   = VAR[ 63];
    CH4[j][i]      = VAR[ 64];
    GAOO[j][i]     = VAR[ 65];
    B3O2[j][i]     = VAR[ 66];
    ACET[j][i]     = VAR[ 67];
    MACRN[j][i]    = VAR[ 68];
    CH2OO[j][i]    = VAR[ 69];
    MGLYOO[j][i]   = VAR[ 70];
    VRO2[j][i]     = VAR[ 71];
    MGLOO[j][i]    = VAR[ 72];
    MACROO[j][i]   = VAR[ 73];
    PO2[j][i]      = VAR[ 74];
    CH3CHOO[j][i]  = VAR[ 75];
    MAN2[j][i]     = VAR[ 76];
    ISNOOA[j][i]   = VAR[ 77];
    H2O2[j][i]     = VAR[ 78];
    PRN1[j][i]     = VAR[ 79];
    ETO2[j][i]     = VAR[ 80];
    KO2[j][i]      = VAR[ 81];
    RCO3[j][i]     = VAR[ 82];
    HC5OO[j][i]    = VAR[ 83];
    GLYC[j][i]     = VAR[ 84];
    ClNO3[j][i]    = VAR[ 85];
    RIO2[j][i]     = VAR[ 86];
    R4N1[j][i]     = VAR[ 87];
    HOCl[j][i]     = VAR[ 88];
    ATO2[j][i]     = VAR[ 89];
    HNO3[j][i]     = VAR[ 90];
    ISN1[j][i]     = VAR[ 91];
    MAO3[j][i]     = VAR[ 92];
    MRO2[j][i]     = VAR[ 93];
    INO2[j][i]     = VAR[ 94];
    HAC[j][i]      = VAR[ 95];
    HC5[j][i]      = VAR[ 96];
    MGLY[j][i]     = VAR[ 97];
    ISOPNBO2[j][i] = VAR[ 98];
    ISOPNDO2[j][i] = VAR[ 99];
    R4O2[j][i]     = VAR[100];
    R4N2[j][i]     = VAR[101];
    BrO[j][i]      = VAR[102];
    RCHO[j][i]     = VAR[103];
    MEK[j][i]      = VAR[104];
    ClO[j][i]      = VAR[105];
    MACR[j][i]     = VAR[106];
    SO2[j][i]      = VAR[107];
    MVK[j][i]      = VAR[108];
    ALD2[j][i]     = VAR[109];
    MCO3[j][i]     = VAR[110];
    CH2O[j][i]     = VAR[111];
    H2O[j][i]      = VAR[112];
    Br[j][i]       = VAR[113];
    NO[j][i]       = VAR[114];
    NO3[j][i]      = VAR[115];
    Cl[j][i]       = VAR[116];
    O[j][i]        = VAR[117];
    O1D[j][i]      = VAR[118];
    O3[j][i]       = VAR[119];
    HO2[j][i]      = VAR[120];
    NO2[j][i]      = VAR[121];
    OH[j][i]       = VAR[122];
    HBr[j][i]      = VAR[123];
    HCl[j][i]      = VAR[124];
    CO[j][i]       = VAR[125];
    MO2[j][i]      = VAR[126];


} /* End of Solution::applyData */

void Solution::applyRing( RealDouble tempArray[],        \
                          const Vector_2Dui &mapIndices, \
                          const UInt iRing )
{

    UInt iNx = 0;
    UInt jNy = 0;

    for ( jNy = 0; jNy < NY; jNy++ ) {
        for ( iNx = 0; iNx < NX; iNx++ ) {
            if ( mapIndices[jNy][iNx] == iRing ) {

                CO2[jNy][iNx]      *= VAR[  0] / tempArray[  0];
                PPN[jNy][iNx]      *= VAR[  1] / tempArray[  1];
                BrNO2[jNy][iNx]    *= VAR[  2] / tempArray[  2];
                IEPOX[jNy][iNx]    *= VAR[  3] / tempArray[  3];
                PMNN[jNy][iNx]     *= VAR[  4] / tempArray[  4];
                N2O[jNy][iNx]      *= VAR[  5] / tempArray[  5];

                /* Make sure that N does not become NaN */
                if ( isinf( VAR[  6] / tempArray[  6] ) ||
                     VAR[  6] / tempArray[  6] >= 1.0E+20 )
                    N[jNy][iNx]     = VAR[  6];
                else
                    N[jNy][iNx]    *= VAR[  6] / tempArray[  6];

                PAN[jNy][iNx]      *= VAR[  7] / tempArray[  7];
                ALK4[jNy][iNx]     *= VAR[  8] / tempArray[  8];
                MAP[jNy][iNx]      *= VAR[  9] / tempArray[  9];
                MPN[jNy][iNx]      *= VAR[ 10] / tempArray[ 10];
                Cl2O2[jNy][iNx]    *= VAR[ 11] / tempArray[ 11];
                ETP[jNy][iNx]      *= VAR[ 12] / tempArray[ 12];
                HNO2[jNy][iNx]     *= VAR[ 13] / tempArray[ 13];
                C3H8[jNy][iNx]     *= VAR[ 14] / tempArray[ 14];
                RA3P[jNy][iNx]     *= VAR[ 15] / tempArray[ 15];
                RB3P[jNy][iNx]     *= VAR[ 16] / tempArray[ 16];
                OClO[jNy][iNx]     *= VAR[ 17] / tempArray[ 17];
                ClNO2[jNy][iNx]    *= VAR[ 18] / tempArray[ 18];
                ISOP[jNy][iNx]     *= VAR[ 19] / tempArray[ 19];
                HNO4[jNy][iNx]     *= VAR[ 20] / tempArray[ 20];
                MAOP[jNy][iNx]     *= VAR[ 21] / tempArray[ 21];
                MP[jNy][iNx]       *= VAR[ 22] / tempArray[ 22];
                ClOO[jNy][iNx]     *= VAR[ 23] / tempArray[ 23];
                RP[jNy][iNx]       *= VAR[ 24] / tempArray[ 24];
                BrCl[jNy][iNx]     *= VAR[ 25] / tempArray[ 25];
                PP[jNy][iNx]       *= VAR[ 26] / tempArray[ 26];
                PRPN[jNy][iNx]     *= VAR[ 27] / tempArray[ 27];
                SO4[jNy][iNx]      *= VAR[ 28] / tempArray[ 28];
                Br2[jNy][iNx]      *= VAR[ 29] / tempArray[ 29];
                ETHLN[jNy][iNx]    *= VAR[ 30] / tempArray[ 30];
                MVKN[jNy][iNx]     *= VAR[ 31] / tempArray[ 31];
                R4P[jNy][iNx]      *= VAR[ 32] / tempArray[ 32];
                C2H6[jNy][iNx]     *= VAR[ 33] / tempArray[ 33];
                RIP[jNy][iNx]      *= VAR[ 34] / tempArray[ 34];
                VRP[jNy][iNx]      *= VAR[ 35] / tempArray[ 35];
                ATOOH[jNy][iNx]    *= VAR[ 36] / tempArray[ 36];
                IAP[jNy][iNx]      *= VAR[ 37] / tempArray[ 37];
                DHMOB[jNy][iNx]    *= VAR[ 38] / tempArray[ 38];
                MOBA[jNy][iNx]     *= VAR[ 39] / tempArray[ 39];
                MRP[jNy][iNx]      *= VAR[ 40] / tempArray[ 40];
                N2O5[jNy][iNx]     *= VAR[ 41] / tempArray[ 41];
                ISNOHOO[jNy][iNx]  *= VAR[ 42] / tempArray[ 42];
                ISNP[jNy][iNx]     *= VAR[ 43] / tempArray[ 43];
                ISOPNB[jNy][iNx]   *= VAR[ 44] / tempArray[ 44];
                IEPOXOO[jNy][iNx]  *= VAR[ 45] / tempArray[ 45];
                MACRNO2[jNy][iNx]  *= VAR[ 46] / tempArray[ 46];
                ROH[jNy][iNx]      *= VAR[ 47] / tempArray[ 47];
                MOBAOO[jNy][iNx]   *= VAR[ 48] / tempArray[ 48];
                DIBOO[jNy][iNx]    *= VAR[ 49] / tempArray[ 49];
                PMN[jNy][iNx]      *= VAR[ 50] / tempArray[ 50];
                ISNOOB[jNy][iNx]   *= VAR[ 51] / tempArray[ 51];
                INPN[jNy][iNx]     *= VAR[ 52] / tempArray[ 52];
                H[jNy][iNx]        *= VAR[ 53] / tempArray[ 53];
                BrNO3[jNy][iNx]    *= VAR[ 54] / tempArray[ 54];
                PRPE[jNy][iNx]     *= VAR[ 55] / tempArray[ 55];
                MVKOO[jNy][iNx]    *= VAR[ 56] / tempArray[ 56];
                Cl2[jNy][iNx]      *= VAR[ 57] / tempArray[ 57];
                ISOPND[jNy][iNx]   *= VAR[ 58] / tempArray[ 58];
                HOBr[jNy][iNx]     *= VAR[ 59] / tempArray[ 59];
                A3O2[jNy][iNx]     *= VAR[ 60] / tempArray[ 60];
                PROPNN[jNy][iNx]   *= VAR[ 61] / tempArray[ 61];
                GLYX[jNy][iNx]     *= VAR[ 62] / tempArray[ 62];
                MAOPO2[jNy][iNx]   *= VAR[ 63] / tempArray[ 63];
                CH4[jNy][iNx]      *= VAR[ 64] / tempArray[ 64];
                GAOO[jNy][iNx]     *= VAR[ 65] / tempArray[ 65];
                B3O2[jNy][iNx]     *= VAR[ 66] / tempArray[ 66];
                ACET[jNy][iNx]     *= VAR[ 67] / tempArray[ 67];
                MACRN[jNy][iNx]    *= VAR[ 68] / tempArray[ 68];
                CH2OO[jNy][iNx]    *= VAR[ 69] / tempArray[ 69];
                MGLYOO[jNy][iNx]   *= VAR[ 70] / tempArray[ 70];
                VRO2[jNy][iNx]     *= VAR[ 71] / tempArray[ 71];
                MGLOO[jNy][iNx]    *= VAR[ 72] / tempArray[ 72];
                MACROO[jNy][iNx]   *= VAR[ 73] / tempArray[ 73];
                PO2[jNy][iNx]      *= VAR[ 74] / tempArray[ 74];
                CH3CHOO[jNy][iNx]  *= VAR[ 75] / tempArray[ 75];
                MAN2[jNy][iNx]     *= VAR[ 76] / tempArray[ 76];
                ISNOOA[jNy][iNx]   *= VAR[ 77] / tempArray[ 77];
                H2O2[jNy][iNx]     *= VAR[ 78] / tempArray[ 78];
                PRN1[jNy][iNx]     *= VAR[ 79] / tempArray[ 79];
                ETO2[jNy][iNx]     *= VAR[ 80] / tempArray[ 80];
                KO2[jNy][iNx]      *= VAR[ 81] / tempArray[ 81];
                RCO3[jNy][iNx]     *= VAR[ 82] / tempArray[ 82];
                HC5OO[jNy][iNx]    *= VAR[ 83] / tempArray[ 83];
                GLYC[jNy][iNx]     *= VAR[ 84] / tempArray[ 84];
                ClNO3[jNy][iNx]    *= VAR[ 85] / tempArray[ 85];
                RIO2[jNy][iNx]     *= VAR[ 86] / tempArray[ 86];
                R4N1[jNy][iNx]     *= VAR[ 87] / tempArray[ 87];
                HOCl[jNy][iNx]     *= VAR[ 88] / tempArray[ 88];
                ATO2[jNy][iNx]     *= VAR[ 89] / tempArray[ 89];
                HNO3[jNy][iNx]     *= VAR[ 90] / tempArray[ 90];
                ISN1[jNy][iNx]     *= VAR[ 91] / tempArray[ 91];
                MAO3[jNy][iNx]     *= VAR[ 92] / tempArray[ 92];
                MRO2[jNy][iNx]     *= VAR[ 93] / tempArray[ 93];
                INO2[jNy][iNx]     *= VAR[ 94] / tempArray[ 94];
                HAC[jNy][iNx]      *= VAR[ 95] / tempArray[ 95];
                HC5[jNy][iNx]      *= VAR[ 96] / tempArray[ 96];
                MGLY[jNy][iNx]     *= VAR[ 97] / tempArray[ 97];
                ISOPNBO2[jNy][iNx] *= VAR[ 98] / tempArray[ 98];
                ISOPNDO2[jNy][iNx] *= VAR[ 99] / tempArray[ 99];
                R4O2[jNy][iNx]     *= VAR[100] / tempArray[100];
                R4N2[jNy][iNx]     *= VAR[101] / tempArray[101];
                BrO[jNy][iNx]      *= VAR[102] / tempArray[102];
                RCHO[jNy][iNx]     *= VAR[103] / tempArray[103];
                MEK[jNy][iNx]      *= VAR[104] / tempArray[104];
                ClO[jNy][iNx]      *= VAR[105] / tempArray[105];
                MACR[jNy][iNx]     *= VAR[106] / tempArray[106];
                SO2[jNy][iNx]      *= VAR[107] / tempArray[107];
                MVK[jNy][iNx]      *= VAR[108] / tempArray[108];
                ALD2[jNy][iNx]     *= VAR[109] / tempArray[109];
                MCO3[jNy][iNx]     *= VAR[110] / tempArray[110];
                CH2O[jNy][iNx]     *= VAR[111] / tempArray[111];
                H2O[jNy][iNx]      *= VAR[112] / tempArray[112];
                Br[jNy][iNx]       *= VAR[113] / tempArray[113];
                NO[jNy][iNx]       *= VAR[114] / tempArray[114];
                NO3[jNy][iNx]      *= VAR[115] / tempArray[115];
                Cl[jNy][iNx]       *= VAR[116] / tempArray[116];

                /* Make sure that O does not become NaN */
                if ( isinf( VAR[117] / tempArray[117] ) ||
                     VAR[117] / tempArray[117] >= 1.0E+20 )
                    O[jNy][iNx]     = VAR[117];
                else
                    O[jNy][iNx]    *= VAR[117] / tempArray[117];

                /* Make sure that O1D does not become NaN */
                if ( isinf( VAR[118] / tempArray[118] ) ||
                     VAR[118] / tempArray[118] >= 1.0E+20 )
                    O1D[jNy][iNx]   = VAR[118];
                else
                    O1D[jNy][iNx]  *= VAR[118] / tempArray[118];

                O3[jNy][iNx]       *= VAR[119] / tempArray[119];
                HO2[jNy][iNx]      *= VAR[120] / tempArray[120];
                NO2[jNy][iNx]      *= VAR[121] / tempArray[121];
                OH[jNy][iNx]       *= VAR[122] / tempArray[122];
                HBr[jNy][iNx]      *= VAR[123] / tempArray[123];
                HCl[jNy][iNx]      *= VAR[124] / tempArray[124];
                CO[jNy][iNx]       *= VAR[125] / tempArray[125];
                MO2[jNy][iNx]      *= VAR[126] / tempArray[126];

            }
        }
    }

} /* End of Solution::applyRing */

void Solution::applyAmbient( const Vector_2Dui &mapIndices, \
                             const UInt iRing )
{

    UInt iNx = 0;
    UInt jNy = 0;

    for ( jNy = 0; jNy < NY; jNy++ ) {
        for ( iNx = 0; iNx < NX; iNx++ ) {
            if ( mapIndices[jNy][iNx] == iRing ) {

                CO2[jNy][iNx]      = VAR[  0];
                PPN[jNy][iNx]      = VAR[  1];
                BrNO2[jNy][iNx]    = VAR[  2];
                IEPOX[jNy][iNx]    = VAR[  3];
                PMNN[jNy][iNx]     = VAR[  4];
                N2O[jNy][iNx]      = VAR[  5];
                N[jNy][iNx]        = VAR[  6];
                PAN[jNy][iNx]      = VAR[  7];
                ALK4[jNy][iNx]     = VAR[  8];
                MAP[jNy][iNx]      = VAR[  9];
                MPN[jNy][iNx]      = VAR[ 10];
                Cl2O2[jNy][iNx]    = VAR[ 11];
                ETP[jNy][iNx]      = VAR[ 12];
                HNO2[jNy][iNx]     = VAR[ 13];
                C3H8[jNy][iNx]     = VAR[ 14];
                RA3P[jNy][iNx]     = VAR[ 15];
                RB3P[jNy][iNx]     = VAR[ 16];
                OClO[jNy][iNx]     = VAR[ 17];
                ClNO2[jNy][iNx]    = VAR[ 18];
                ISOP[jNy][iNx]     = VAR[ 19];
                HNO4[jNy][iNx]     = VAR[ 20];
                MAOP[jNy][iNx]     = VAR[ 21];
                MP[jNy][iNx]       = VAR[ 22];
                ClOO[jNy][iNx]     = VAR[ 23];
                RP[jNy][iNx]       = VAR[ 24];
                BrCl[jNy][iNx]     = VAR[ 25];
                PP[jNy][iNx]       = VAR[ 26];
                PRPN[jNy][iNx]     = VAR[ 27];
                SO4[jNy][iNx]      = VAR[ 28];
                Br2[jNy][iNx]      = VAR[ 29];
                ETHLN[jNy][iNx]    = VAR[ 30];
                MVKN[jNy][iNx]     = VAR[ 31];
                R4P[jNy][iNx]      = VAR[ 32];
                C2H6[jNy][iNx]     = VAR[ 33];
                RIP[jNy][iNx]      = VAR[ 34];
                VRP[jNy][iNx]      = VAR[ 35];
                ATOOH[jNy][iNx]    = VAR[ 36];
                IAP[jNy][iNx]      = VAR[ 37];
                DHMOB[jNy][iNx]    = VAR[ 38];
                MOBA[jNy][iNx]     = VAR[ 39];
                MRP[jNy][iNx]      = VAR[ 40];
                N2O5[jNy][iNx]     = VAR[ 41];
                ISNOHOO[jNy][iNx]  = VAR[ 42];
                ISNP[jNy][iNx]     = VAR[ 43];
                ISOPNB[jNy][iNx]   = VAR[ 44];
                IEPOXOO[jNy][iNx]  = VAR[ 45];
                MACRNO2[jNy][iNx]  = VAR[ 46];
                ROH[jNy][iNx]      = VAR[ 47];
                MOBAOO[jNy][iNx]   = VAR[ 48];
                DIBOO[jNy][iNx]    = VAR[ 49];
                PMN[jNy][iNx]      = VAR[ 50];
                ISNOOB[jNy][iNx]   = VAR[ 51];
                INPN[jNy][iNx]     = VAR[ 52];
                H[jNy][iNx]        = VAR[ 53];
                BrNO3[jNy][iNx]    = VAR[ 54];
                PRPE[jNy][iNx]     = VAR[ 55];
                MVKOO[jNy][iNx]    = VAR[ 56];
                Cl2[jNy][iNx]      = VAR[ 57];
                ISOPND[jNy][iNx]   = VAR[ 58];
                HOBr[jNy][iNx]     = VAR[ 59];
                A3O2[jNy][iNx]     = VAR[ 60];
                PROPNN[jNy][iNx]   = VAR[ 61];
                GLYX[jNy][iNx]     = VAR[ 62];
                MAOPO2[jNy][iNx]   = VAR[ 63];
                CH4[jNy][iNx]      = VAR[ 64];
                GAOO[jNy][iNx]     = VAR[ 65];
                B3O2[jNy][iNx]     = VAR[ 66];
                ACET[jNy][iNx]     = VAR[ 67];
                MACRN[jNy][iNx]    = VAR[ 68];
                CH2OO[jNy][iNx]    = VAR[ 69];
                MGLYOO[jNy][iNx]   = VAR[ 70];
                VRO2[jNy][iNx]     = VAR[ 71];
                MGLOO[jNy][iNx]    = VAR[ 72];
                MACROO[jNy][iNx]   = VAR[ 73];
                PO2[jNy][iNx]      = VAR[ 74];
                CH3CHOO[jNy][iNx]  = VAR[ 75];
                MAN2[jNy][iNx]     = VAR[ 76];
                ISNOOA[jNy][iNx]   = VAR[ 77];
                H2O2[jNy][iNx]     = VAR[ 78];
                PRN1[jNy][iNx]     = VAR[ 79];
                ETO2[jNy][iNx]     = VAR[ 80];
                KO2[jNy][iNx]      = VAR[ 81];
                RCO3[jNy][iNx]     = VAR[ 82];
                HC5OO[jNy][iNx]    = VAR[ 83];
                GLYC[jNy][iNx]     = VAR[ 84];
                ClNO3[jNy][iNx]    = VAR[ 85];
                RIO2[jNy][iNx]     = VAR[ 86];
                R4N1[jNy][iNx]     = VAR[ 87];
                HOCl[jNy][iNx]     = VAR[ 88];
                ATO2[jNy][iNx]     = VAR[ 89];
                HNO3[jNy][iNx]     = VAR[ 90];
                ISN1[jNy][iNx]     = VAR[ 91];
                MAO3[jNy][iNx]     = VAR[ 92];
                MRO2[jNy][iNx]     = VAR[ 93];
                INO2[jNy][iNx]     = VAR[ 94];
                HAC[jNy][iNx]      = VAR[ 95];
                HC5[jNy][iNx]      = VAR[ 96];
                MGLY[jNy][iNx]     = VAR[ 97];
                ISOPNBO2[jNy][iNx] = VAR[ 98];
                ISOPNDO2[jNy][iNx] = VAR[ 99];
                R4O2[jNy][iNx]     = VAR[100];
                R4N2[jNy][iNx]     = VAR[101];
                BrO[jNy][iNx]      = VAR[102];
                RCHO[jNy][iNx]     = VAR[103];
                MEK[jNy][iNx]      = VAR[104];
                ClO[jNy][iNx]      = VAR[105];
                MACR[jNy][iNx]     = VAR[106];
                SO2[jNy][iNx]      = VAR[107];
                MVK[jNy][iNx]      = VAR[108];
                ALD2[jNy][iNx]     = VAR[109];
                MCO3[jNy][iNx]     = VAR[110];
                CH2O[jNy][iNx]     = VAR[111];
                H2O[jNy][iNx]      = VAR[112];
                Br[jNy][iNx]       = VAR[113];
                NO[jNy][iNx]       = VAR[114];
                NO3[jNy][iNx]      = VAR[115];
                Cl[jNy][iNx]       = VAR[116];
                O[jNy][iNx]        = VAR[117];
                O1D[jNy][iNx]      = VAR[118];
                O3[jNy][iNx]       = VAR[119];
                HO2[jNy][iNx]      = VAR[120];
                NO2[jNy][iNx]      = VAR[121];
                OH[jNy][iNx]       = VAR[122];
                HBr[jNy][iNx]      = VAR[123];
                HCl[jNy][iNx]      = VAR[124];
                CO[jNy][iNx]       = VAR[125];
                MO2[jNy][iNx]      = VAR[126];

            }
        }
    }

} /* End of Solution::applyAmbient */


void Solution::addEmission( const Emission &EI, const Aircraft &AC,        \
                            const Mesh &m,                                 \
                            bool halfRing,                                 \
                            const double temperature, bool set2Saturation, \
                            AIM::Aerosol liqAer, AIM::Aerosol iceAer,      \
                            const double Soot_Den,                         \
                            const Meteorology &met, const RealDouble areaPlume )
{
    /* TODO: Release as Gaussian instead of top-hat? */

    UInt innerRing;
    UInt iNx, jNy;
    RealDouble w;

    RealDouble E_CO2, E_H2O, E_NO, E_NO2, E_HNO2, E_SO2, E_CO, E_CH4, E_C2H6, \
               E_PRPE, E_ALK4, E_CH2O, E_ALD2, E_GLYX, E_MGLY;
    RealDouble E_Soot;
    const RealDouble rad = EI.getSootRad();
    const RealDouble fuelPerDist = AC.FuelFlow() / AC.VFlight();

    /* Unit check:  [kg/m]   =   [kg fuel/s] /     [m/s] */
    E_CO2  = EI.getCO2()  / ( MW_CO2  * 1.0E+03 ) * fuelPerDist * physConst::Na;
    /*     = [g(CO2)/kg f]/ ( [kg/mol]* [g/kg]  ) * [kg fuel/m] * [molec/mol]
     *     = [molec/m]
     */
    E_NO   = EI.getNO()   / ( MW_NO   * 1.0E+03 ) * fuelPerDist * physConst::Na;
    /*     = [g(NO)/kg f] / ( g(NO)/mol         ) * [kg fuel/m] * [molec/mol]
     *     = [molec/m] */
    E_NO2  = EI.getNO2()  / ( MW_NO2  * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_HNO2 = EI.getHNO2() / ( MW_HNO2 * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_CO   = EI.getCO()   / ( MW_CO   * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_CH4  = EI.getCH4()  / ( MW_CH4  * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_C2H6 = EI.getC2H6() / ( MW_C2H6 * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_PRPE = EI.getPRPE() / ( MW_PRPE * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_ALK4 = EI.getALK4() / ( MW_ALK4 * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_CH2O = EI.getCH2O() / ( MW_CH2O * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_ALD2 = EI.getALD2() / ( MW_ALD2 * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_GLYX = EI.getGLYX() / ( MW_GLYX * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_MGLY = EI.getMGLY() / ( MW_MGLY * 1.0E+03 ) * fuelPerDist * physConst::Na;
    if ( !set2Saturation ) {
        E_H2O  = EI.getH2O()  / ( MW_H2O  * 1.0E+03 ) * fuelPerDist * physConst::Na;
    }
    E_SO2  = ( 1.0 - SO2TOSO4 ) * \
             EI.getSO2()  / ( MW_SO2  * 1.0E+03 ) * fuelPerDist * physConst::Na;

    E_Soot = EI.getSoot() / ( 4.0 / 3.0 * physConst::PI * physConst::RHO_SOOT * 1.00E+03 * rad * rad * rad ) * fuelPerDist;
    /*     = [g_soot/kg_fuel]/ (                        * [kg_soot/m^3]       * [g/kg]   * [m^3]           ) * [kg_fuel/m]
     *     = [part/m]
     */

    Vector_2D cellAreas  = m.areas();
    Vector_3D weights    = m.weights;
    Vector_1Dui nCellMap = m.nMap();
    UInt nCell   = 0;

    /* Apply area scaling for aerosols before releasing into the grid
     * *Aer are initially in #/cm^3
     * Scaling them the initial plume area such that they are now in #*m^2/cm^3
     * When applying them to the grid, they are scaled by the total ring area
     * in Grid_Aerosol::addPDF */
    liqAer.scalePdf( areaPlume );
    iceAer.scalePdf( areaPlume );


    if ( !halfRing ) {
        /* Full rings */
        innerRing = 0;
        nCell     = nCellMap[innerRing];
        for ( jNy = 0; jNy < NY; jNy++ ) {
            for ( iNx = 0; iNx < NX; iNx++ ) {

                H2O_met[jNy][iNx] = met.H2O_[jNy][iNx];

                w = weights[innerRing][jNy][iNx];

                /* Initially weights are either 0 or 1 */
                if ( w != 0.0E+00 ) {

                    if ( !reducedSize ) {
                        /* Emissions will only be added if chemistry is turned on! */

                        /* Each of the following arrays are in [# / cm^3], where
                         * # can either be expressed in molecules for gaseous
                         * species or in number of particles for aerosols */
                        CO2[jNy][iNx]  += ( E_CO2  * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        NO[jNy][iNx]   += ( E_NO   * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        NO2[jNy][iNx]  += ( E_NO2  * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        HNO2[jNy][iNx] += ( E_HNO2 * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        CO[jNy][iNx]   += ( E_CO   * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        CH4[jNy][iNx]  += ( E_CH4  * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        C2H6[jNy][iNx] += ( E_C2H6 * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        PRPE[jNy][iNx] += ( E_PRPE * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        ALK4[jNy][iNx] += ( E_ALK4 * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        CH2O[jNy][iNx] += ( E_CH2O * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        ALD2[jNy][iNx] += ( E_ALD2 * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        GLYX[jNy][iNx] += ( E_GLYX * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        MGLY[jNy][iNx] += ( E_MGLY * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        SO2[jNy][iNx]  += ( E_SO2  * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                    }

                    if ( set2Saturation ) {
                        /* If supersaturated, then set water vapor to saturation and no bare soot particles
                         * as they are all covered with ice */
                        H2O_plume[jNy][iNx] += physFunc::pSat_H2Os( met.temp_[jNy][iNx] ) / ( physConst::kB * met.temp_[jNy][iNx] * 1.00E+06 ) - H2O_met[jNy][iNx]; /* [molec / cm^3] */
                    } else {
                        /* If subsaturated, then emit water and soot */
                        H2O_plume[jNy][iNx]      += ( E_H2O * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) ); /* [molec / cm^3] */
                        sootDens[jNy][iNx] += ( Soot_Den * areaPlume / ( nCell * cellAreas[jNy][iNx] ) );
                        sootRadi[jNy][iNx] = rad;
                        sootArea[jNy][iNx] = 4.0 * physConst::PI * rad * rad * sootDens[jNy][iNx];
                    }

                }

            }

        }

        if ( ( std::isfinite(iceAer.Moment()) ) && ( iceAer.Moment() > 0.0E+00 ) )
            solidAerosol.addPDF( iceAer, weights[innerRing], cellAreas, nCell );
        if ( ( std::isfinite(liqAer.Moment()) ) && ( liqAer.Moment() > 0.0E+00 ) )
            liquidAerosol.addPDF( liqAer, weights[innerRing], cellAreas, nCell );


    } else {
        /* Half rings */

        for ( innerRing = 0; innerRing <= 1; innerRing++ ) {
            /* Concentrations should be identical whether it is half-rings or
             * full rings. Therefore, make sure that nCell is doubled when
             * using half-rings! */
            nCell     = 2.0 * nCellMap[innerRing];
            for ( jNy = 0; jNy < NY; jNy++ ) {
                for ( iNx = 0; iNx < NX; iNx++ ) {

                    H2O_met[jNy][iNx] = met.H2O_[jNy][iNx];

                    w = weights[innerRing][jNy][iNx];

                    /* Initially weights are either 0 or 1 */
                    if ( w != 0.0E+00 ) {
                        if ( !reducedSize ) {
                            /* Emissions will only be added if chemistry is turned on! */

                            /* Each of the following arrays are in [# / cm^3], where
                             * # can either be expressed in molecules for gaseous
                             * species or in number of particles for aerosols */
                            CO2[jNy][iNx]  += ( E_CO2  * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            NO[jNy][iNx]   += ( E_NO   * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            NO2[jNy][iNx]  += ( E_NO2  * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            HNO2[jNy][iNx] += ( E_HNO2 * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            CO[jNy][iNx]   += ( E_CO   * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            CH4[jNy][iNx]  += ( E_CH4  * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            C2H6[jNy][iNx] += ( E_C2H6 * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            PRPE[jNy][iNx] += ( E_PRPE * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            ALK4[jNy][iNx] += ( E_ALK4 * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            CH2O[jNy][iNx] += ( E_CH2O * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            ALD2[jNy][iNx] += ( E_ALD2 * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            GLYX[jNy][iNx] += ( E_GLYX * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            MGLY[jNy][iNx] += ( E_MGLY * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            SO2[jNy][iNx]  += ( E_SO2  * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        }

                        if ( set2Saturation ) {
                            /* If supersaturated, then set water vapor to saturation and no
                             * bare soot particles as they are all covered with ice */
                            H2O_plume[jNy][iNx] += physFunc::pSat_H2Os( met.temp_[jNy][iNx] ) / ( physConst::kB * met.temp_[jNy][iNx] * 1.00E+06 ) - H2O_met[jNy][iNx]; /* [molec / cm^3] */
                        } else {
                            /* If subsaturated, then emit water and soot */
                            H2O_plume[jNy][iNx]      += ( E_H2O * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) ); /* [molec / cm^3] */
                            sootDens[jNy][iNx] += ( Soot_Den * areaPlume / ( nCell * cellAreas[jNy][iNx] ) );
                            sootRadi[jNy][iNx] = rad;
                            sootArea[jNy][iNx] = 4.0 * physConst::PI * rad * rad * sootDens[jNy][iNx];
                        }

                    }
                }
            }

            if ( ( std::isfinite(iceAer.Moment()) ) && ( iceAer.Moment() > 0.0E+00 ) )
                solidAerosol.addPDF( iceAer, weights[innerRing], cellAreas, nCell );
            if ( ( std::isfinite(liqAer.Moment()) ) && ( liqAer.Moment() > 0.0E+00 ) )
                liquidAerosol.addPDF( liqAer, weights[innerRing], cellAreas, nCell );

        }

    }

    for ( jNy = 0; jNy < NY; jNy++ ) {
        for ( iNx = 0; iNx < NX; iNx++ ) {
            H2O[jNy][iNx] = H2O_met[jNy][iNx] + H2O_plume[jNy][iNx];
        }
    }

} /* End of Solution::addEmission */

Vector_1D Solution::getAmbient() const
{

    Vector_1D ambVector( NSPEC, 0.0 );

    ambVector[  0] = CO2[0][0];
    ambVector[  1] = PPN[0][0];
    ambVector[  2] = BrNO2[0][0];
    ambVector[  3] = IEPOX[0][0];
    ambVector[  4] = PMNN[0][0];
    ambVector[  5] = N2O[0][0];
    ambVector[  6] = N[0][0];
    ambVector[  7] = PAN[0][0];
    ambVector[  8] = ALK4[0][0];
    ambVector[  9] = MAP[0][0];
    ambVector[ 10] = MPN[0][0];
    ambVector[ 11] = Cl2O2[0][0];
    ambVector[ 12] = ETP[0][0];
    ambVector[ 13] = HNO2[0][0];
    ambVector[ 14] = C3H8[0][0];
    ambVector[ 15] = RA3P[0][0];
    ambVector[ 16] = RB3P[0][0];
    ambVector[ 17] = OClO[0][0];
    ambVector[ 18] = ClNO2[0][0];
    ambVector[ 19] = ISOP[0][0];
    ambVector[ 20] = HNO4[0][0];
    ambVector[ 21] = MAOP[0][0];
    ambVector[ 22] = MP[0][0];
    ambVector[ 23] = ClOO[0][0];
    ambVector[ 24] = RP[0][0];
    ambVector[ 25] = BrCl[0][0];
    ambVector[ 26] = PP[0][0];
    ambVector[ 27] = PRPN[0][0];
    ambVector[ 28] = SO4[0][0];
    ambVector[ 29] = Br2[0][0];
    ambVector[ 30] = ETHLN[0][0];
    ambVector[ 31] = MVKN[0][0];
    ambVector[ 32] = R4P[0][0];
    ambVector[ 33] = C2H6[0][0];
    ambVector[ 34] = RIP[0][0];
    ambVector[ 35] = VRP[0][0];
    ambVector[ 36] = ATOOH[0][0];
    ambVector[ 37] = IAP[0][0];
    ambVector[ 38] = DHMOB[0][0];
    ambVector[ 39] = MOBA[0][0];
    ambVector[ 40] = MRP[0][0];
    ambVector[ 41] = N2O5[0][0];
    ambVector[ 42] = ISNOHOO[0][0];
    ambVector[ 43] = ISNP[0][0];
    ambVector[ 44] = ISOPNB[0][0];
    ambVector[ 45] = IEPOXOO[0][0];
    ambVector[ 46] = MACRNO2[0][0];
    ambVector[ 47] = ROH[0][0];
    ambVector[ 48] = MOBAOO[0][0];
    ambVector[ 49] = DIBOO[0][0];
    ambVector[ 50] = PMN[0][0];
    ambVector[ 51] = ISNOOB[0][0];
    ambVector[ 52] = INPN[0][0];
    ambVector[ 53] = H[0][0];
    ambVector[ 54] = BrNO3[0][0];
    ambVector[ 55] = PRPE[0][0];
    ambVector[ 56] = MVKOO[0][0];
    ambVector[ 57] = Cl2[0][0];
    ambVector[ 58] = ISOPND[0][0];
    ambVector[ 59] = HOBr[0][0];
    ambVector[ 60] = A3O2[0][0];
    ambVector[ 61] = PROPNN[0][0];
    ambVector[ 62] = GLYX[0][0];
    ambVector[ 63] = MAOPO2[0][0];
    ambVector[ 64] = CH4[0][0];
    ambVector[ 65] = GAOO[0][0];
    ambVector[ 66] = B3O2[0][0];
    ambVector[ 67] = ACET[0][0];
    ambVector[ 68] = MACRN[0][0];
    ambVector[ 69] = CH2OO[0][0];
    ambVector[ 70] = MGLYOO[0][0];
    ambVector[ 71] = VRO2[0][0];
    ambVector[ 72] = MGLOO[0][0];
    ambVector[ 73] = MACROO[0][0];
    ambVector[ 74] = PO2[0][0];
    ambVector[ 75] = CH3CHOO[0][0];
    ambVector[ 76] = MAN2[0][0];
    ambVector[ 77] = ISNOOA[0][0];
    ambVector[ 78] = H2O2[0][0];
    ambVector[ 79] = PRN1[0][0];
    ambVector[ 80] = ETO2[0][0];
    ambVector[ 81] = KO2[0][0];
    ambVector[ 82] = RCO3[0][0];
    ambVector[ 83] = HC5OO[0][0];
    ambVector[ 84] = GLYC[0][0];
    ambVector[ 85] = ClNO3[0][0];
    ambVector[ 86] = RIO2[0][0];
    ambVector[ 87] = R4N1[0][0];
    ambVector[ 88] = HOCl[0][0];
    ambVector[ 89] = ATO2[0][0];
    ambVector[ 90] = HNO3[0][0];
    ambVector[ 91] = ISN1[0][0];
    ambVector[ 92] = MAO3[0][0];
    ambVector[ 93] = MRO2[0][0];
    ambVector[ 94] = INO2[0][0];
    ambVector[ 95] = HAC[0][0];
    ambVector[ 96] = HC5[0][0];
    ambVector[ 97] = MGLY[0][0];
    ambVector[ 98] = ISOPNBO2[0][0];
    ambVector[ 99] = ISOPNDO2[0][0];
    ambVector[100] = R4O2[0][0];
    ambVector[101] = R4N2[0][0];
    ambVector[102] = BrO[0][0];
    ambVector[103] = RCHO[0][0];
    ambVector[104] = MEK[0][0];
    ambVector[105] = ClO[0][0];
    ambVector[106] = MACR[0][0];
    ambVector[107] = SO2[0][0];
    ambVector[108] = MVK[0][0];
    ambVector[109] = ALD2[0][0];
    ambVector[110] = MCO3[0][0];
    ambVector[111] = CH2O[0][0];
    ambVector[112] = H2O[0][0];
    ambVector[113] = Br[0][0];
    ambVector[114] = NO[0][0];
    ambVector[115] = NO3[0][0];
    ambVector[116] = Cl[0][0];
    ambVector[117] = O[0][0];
    ambVector[118] = O1D[0][0];
    ambVector[119] = O3[0][0];
    ambVector[120] = HO2[0][0];
    ambVector[121] = NO2[0][0];
    ambVector[122] = OH[0][0];
    ambVector[123] = HBr[0][0];
    ambVector[124] = HCl[0][0];
    ambVector[125] = CO[0][0];
    ambVector[126] = MO2[0][0];
    ambVector[127] = ACTA[0][0];
    ambVector[128] = EOH[0][0];
    ambVector[129] = H2[0][0];
    ambVector[130] = HCOOH[0][0];
    ambVector[131] = MOH[0][0];
    ambVector[132] = N2[0][0];
    ambVector[133] = O2[0][0];
    ambVector[134] = RCOOH[0][0];

    return ambVector;

} /* End of Solution::getAmbient */

Vector_1D Solution::getLiqSpecies( ) const
{

    Vector_1D liqAerVector( 9, 0.0 );

    liqAerVector[0] = SO4L[0][0];
    liqAerVector[1] = H2OL[0][0];
    liqAerVector[2] = HNO3L[0][0];
    liqAerVector[3] = HClL[0][0];
    liqAerVector[4] = HOClL[0][0];
    liqAerVector[5] = HBrL[0][0];
    liqAerVector[6] = HOBrL[0][0];
    liqAerVector[7] = H2OS[0][0];
    liqAerVector[8] = HNO3S[0][0];

    return liqAerVector;

} /* End of Solution::getLiqSpecies */

Vector_2D Solution::getAerosol( ) const
{

    Vector_2D aerVector( nAer, Vector_1D( 3, 0.0 ) );

    aerVector[  0][0] = sootDens[0][0];
    aerVector[  0][1] = sootRadi[0][0];
    aerVector[  0][2] = sootArea[0][0];
    aerVector[  1][0] = PA_nDens;
    aerVector[  1][1] = PA_rEff;
    aerVector[  1][2] = PA_SAD;
    aerVector[  2][0] = LA_nDens;
    aerVector[  2][1] = LA_rEff;
    aerVector[  2][2] = LA_SAD;

    return aerVector;

} /* End of Solution::getAerosol */

Vector_1D Solution::getAerosolDens( ) const
{

    Vector_1D aerVector( nAer, 0.0 );

    aerVector[  0] = sootDens[0][0];
    aerVector[  1] = PA_nDens;
    aerVector[  2] = LA_nDens;

    return aerVector;

} /* End of Solution::getAerosolDens */

Vector_1D Solution::getAerosolRadi( ) const
{

    Vector_1D aerVector( nAer, 0.0 );

    aerVector[  0] = sootRadi[0][0];
    aerVector[  1] = PA_rEff;
    aerVector[  2] = LA_rEff;

    return aerVector;

} /* End of Solution::getAerosolRadi */

Vector_1D Solution::getAerosolArea( ) const
{

    Vector_1D aerVector( nAer, 0.0 );

    aerVector[  0] = sootArea[0][0];
    aerVector[  1] = PA_SAD;
    aerVector[  2] = LA_SAD;

    return aerVector;

} /* End of Solution::getAerosolArea */

void Solution::getAerosolProp( RealDouble ( &radi )[4], \
                               RealDouble ( &area )[4], \
                               RealDouble &IWC,         \
                               const Vector_2D &weights ) const
{

    UInt jNy, iNx;

    Vector_1D aerosolProp( 4, 0.0E+00 );
    RealDouble totalWeight = 0.0E+00;

    for ( jNy = 0; jNy < NY; jNy++ ) {
        for ( iNx = 0; iNx < NX; iNx++ )
            totalWeight += weights[jNy][iNx];
    }

    /* Compute aerosol microphysical properties for ice/NAT */
    aerosolProp = solidAerosol.Average( weights, totalWeight );

    radi[0] = aerosolProp[1];                      /* [m]        */
    area[0] = aerosolProp[2];                      /* [m^2/cm^3] */
    IWC     = physConst::RHO_ICE * aerosolProp[3]; /* [kg/cm^3] */

    /* Compute aerosol microphysical properties for stratospheric liquid aerosols */
    aerosolProp = liquidAerosol.Average( weights, totalWeight );

    radi[1] = aerosolProp[1]; /* [m]        */
    area[1] = aerosolProp[2]; /* [m^2/cm^3] */

    /* Compute aerosol microphysical properties for tropospheric sulfates (near ground) */

    radi[2] = 0.0E+00;
    area[2] = 0.0E+00;

    /* Compute aerosol microphysical properties for BC particles */

    radi[3] = 0.0E+00;
    area[3] = 0.0E+00;

    for ( jNy = 0; jNy < NY; jNy++ ) {
        for ( iNx = 0; iNx < NX; iNx++ ) {
            area[3] += sootDens[jNy][iNx] * 4.0 * physConst::PI * \
                       sootRadi[jNy][iNx] * sootRadi[jNy][iNx] *  \
                       weights[jNy][iNx] / totalWeight;
            radi[3] += sootRadi[jNy][iNx] * \
                       weights[jNy][iNx] / totalWeight;
        }
    }

} /* End of Solution::getAerosolProp */

int Solution::SpinUp( Vector_1D &amb_Value,       \
                      const Input &input,         \
                      const RealDouble airDens,   \
                      const RealDouble startTime, \
                      const bool DBG )
{

    /* Chemistry timestep
     * DT_CHEM               = 10 mins */
    const RealDouble DT_CHEM = 10.0 * 60.0;
    RealDouble curr_Time_s   = startTime * 3600.0;

    /* Integrate chemistry from startTime until endTime
     * Make sure than endTime is greater than startTime,
     * if not integrate until next day at the same time */
    RealDouble RunUntil = input.emissionTime();

//    /* If emission time corresponds to data from file exit here */
//    if ( startTime == RunUntil )
//        return 1;

//    if ( RunUntil < startTime )
//        RunUntil += 24.0;
    RunUntil += 24.0;

    /* Convert to seconds */
    RunUntil *= 3600.0;

    /* Allocate arrays for KPP */
    int IERR = 0;
    RealDouble STEPMIN = (RealDouble)0.0;

    RealDouble RTOL[NVAR];
    RealDouble ATOL[NVAR];

    for( UInt i = 0; i < NVAR; i++ ) {
        RTOL[i] = KPP_RTOLS;
        ATOL[i] = KPP_ATOLS;
    }


    /* Initialize arrays */
    for ( UInt iVar = 0; iVar < NVAR; iVar++ )
        VAR[iVar] = amb_Value[iVar] * airDens;

    for ( UInt iFix = 0; iFix < NFIX; iFix++ )
        FIX[iFix] = amb_Value[NVAR+iFix] * airDens;

    /* Define sun parameters */
    SZA *sun = new SZA( input.latitude_deg(), input.emissionDOY() );

    if ( DBG )
        std::cout << "\n Running spin-up from " << curr_Time_s / 3600.0 << " to " << RunUntil / 3600.0 << " [hr]\n";

    while ( curr_Time_s < RunUntil ) {

        /* Compute the cosize of solar zenith angle midway through the integration step */
        sun->Update( curr_Time_s + DT_CHEM/2 );

        for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
            PHOTOL[iPhotol] = 0.0E+00;

        if ( sun->CSZA > 0.0E+00 )
            Update_JRates( PHOTOL, sun->CSZA );

        if ( DBG ) {
            std::cout << "\n DEBUG : (In SpinUp)\n";
            for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
                std::cout << "         PHOTOL[" << iPhotol << "] = " << PHOTOL[iPhotol] << "\n";
        }

        /* Update reaction rates */
        for ( UInt iReact = 0; iReact < NREACT; iReact++ )
            RCONST[iReact] = 0.0E+00;

        Update_RCONST( input.temperature_K(), input.pressure_Pa(), airDens, VAR[ind_H2O] );

        /* ~~~~~~~~~~~~~~~~~~~~~~~~ */
        /* ~~~~~ Integration ~~~~~~ */
        /* ~~~~~~~~~~~~~~~~~~~~~~~~ */

        IERR = INTEGRATE( VAR, curr_Time_s, curr_Time_s + DT_CHEM, \
                          ATOL, RTOL, STEPMIN );

        if ( IERR < 0 ) {
            /* Integration failed */

            std::cout << " SpinUp Integration failed";
            #ifdef OMP
                std::cout << " on " << omp_get_thread_num();
            #endif /* OMP */
            std::cout << " at time t = " << curr_Time_s/3600.0 << "\n";

            if ( DBG ) {
                std::cout << " ~~~ Printing reaction rates:\n";
                for ( UInt iReact = 0; iReact < NREACT; iReact++ ) {
                    std::cout << "Reaction " << iReact << ": " << RCONST[iReact] << " [molec/cm^3/s]\n";
                }
                std::cout << " ~~~ Printing concentrations:\n";
                for ( UInt iSpec = 0; iSpec < NVAR; iSpec++ ) {
                    std::cout << "Species " << iSpec << ": " << VAR[iSpec]/airDens*1.0E+09 << " [ppb]\n";
                }
            }

            return KPP_FAIL;
        }

        curr_Time_s += DT_CHEM;

    }

    for ( UInt iVar = 0; iVar < NVAR; iVar++ )
        amb_Value[iVar] = VAR[iVar] / airDens;


    /* Clear dynamically allocated variable(s) */
    sun->~SZA();

    return IERR;

} /* End of Solution::SpinUp */

void Solution::Debug( const RealDouble airDens )
{

    UInt iNx, jNy;
    iNx = 0;
    jNy = 0;

    std::cout << std::endl;
    std::cout << "**** Input Debugger ****" << std::endl;
    std::cout << "Background concentrations: " << std::endl;
    std::cout << std::endl;
    std::cout << "   Species " << "    Value" << " Units " << std::endl;

    std::cout << std::setw(9);
    std::cout << "CO2" << ": ";
    std::cout << std::setw(9);
    std::cout << CO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PPN" << ": ";
    std::cout << std::setw(9);
    std::cout << PPN[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "BrNO2" << ": ";
    std::cout << std::setw(9);
    std::cout << BrNO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "IEPOX" << ": ";
    std::cout << std::setw(9);
    std::cout << IEPOX[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PMNN" << ": ";
    std::cout << std::setw(9);
    std::cout << PMNN[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "N2O" << ": ";
    std::cout << std::setw(9);
    std::cout << N2O[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "N" << ": ";
    std::cout << std::setw(9);
    std::cout << N[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PAN" << ": ";
    std::cout << std::setw(9);
    std::cout << PAN[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ALK4" << ": ";
    std::cout << std::setw(9);
    std::cout << ALK4[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MAP" << ": ";
    std::cout << std::setw(9);
    std::cout << MAP[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MPN" << ": ";
    std::cout << std::setw(9);
    std::cout << MPN[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "Cl2O2" << ": ";
    std::cout << std::setw(9);
    std::cout << Cl2O2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ETP" << ": ";
    std::cout << std::setw(9);
    std::cout << ETP[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HNO2" << ": ";
    std::cout << std::setw(9);
    std::cout << HNO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "C3H8" << ": ";
    std::cout << std::setw(9);
    std::cout << C3H8[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RA3P" << ": ";
    std::cout << std::setw(9);
    std::cout << RA3P[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RB3P" << ": ";
    std::cout << std::setw(9);
    std::cout << RB3P[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "OClO" << ": ";
    std::cout << std::setw(9);
    std::cout << OClO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ClNO2" << ": ";
    std::cout << std::setw(9);
    std::cout << ClNO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISOP" << ": ";
    std::cout << std::setw(9);
    std::cout << ISOP[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HNO4" << ": ";
    std::cout << std::setw(9);
    std::cout << HNO4[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MAOP" << ": ";
    std::cout << std::setw(9);
    std::cout << MAOP[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MP" << ": ";
    std::cout << std::setw(9);
    std::cout << MP[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ClOO" << ": ";
    std::cout << std::setw(9);
    std::cout << ClOO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RP" << ": ";
    std::cout << std::setw(9);
    std::cout << RP[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "BrCl" << ": ";
    std::cout << std::setw(9);
    std::cout << BrCl[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PP" << ": ";
    std::cout << std::setw(9);
    std::cout << PP[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PRPN" << ": ";
    std::cout << std::setw(9);
    std::cout << PRPN[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "SO4" << ": ";
    std::cout << std::setw(9);
    std::cout << SO4[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "Br2" << ": ";
    std::cout << std::setw(9);
    std::cout << Br2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ETHLN" << ": ";
    std::cout << std::setw(9);
    std::cout << ETHLN[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MVKN" << ": ";
    std::cout << std::setw(9);
    std::cout << MVKN[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "R4P" << ": ";
    std::cout << std::setw(9);
    std::cout << R4P[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "C2H6" << ": ";
    std::cout << std::setw(9);
    std::cout << C2H6[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RIP" << ": ";
    std::cout << std::setw(9);
    std::cout << RIP[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "VRP" << ": ";
    std::cout << std::setw(9);
    std::cout << VRP[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ATOOH" << ": ";
    std::cout << std::setw(9);
    std::cout << ATOOH[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "IAP" << ": ";
    std::cout << std::setw(9);
    std::cout << IAP[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "DHMOB" << ": ";
    std::cout << std::setw(9);
    std::cout << DHMOB[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MOBA" << ": ";
    std::cout << std::setw(9);
    std::cout << MOBA[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MRP" << ": ";
    std::cout << std::setw(9);
    std::cout << MRP[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "N2O5" << ": ";
    std::cout << std::setw(9);
    std::cout << N2O5[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISNOHOO" << ": ";
    std::cout << std::setw(9);
    std::cout << ISNOHOO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISNP" << ": ";
    std::cout << std::setw(9);
    std::cout << ISNP[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISOPNB" << ": ";
    std::cout << std::setw(9);
    std::cout << ISOPNB[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "IEPOXOO" << ": ";
    std::cout << std::setw(9);
    std::cout << IEPOXOO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MACRNO2" << ": ";
    std::cout << std::setw(9);
    std::cout << MACRNO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ROH" << ": ";
    std::cout << std::setw(9);
    std::cout << ROH[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MOBAOO" << ": ";
    std::cout << std::setw(9);
    std::cout << MOBAOO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "DIBOO" << ": ";
    std::cout << std::setw(9);
    std::cout << DIBOO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PMN" << ": ";
    std::cout << std::setw(9);
    std::cout << PMN[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISNOOB" << ": ";
    std::cout << std::setw(9);
    std::cout << ISNOOB[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "INPN" << ": ";
    std::cout << std::setw(9);
    std::cout << INPN[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "H" << ": ";
    std::cout << std::setw(9);
    std::cout << H[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "BrNO3" << ": ";
    std::cout << std::setw(9);
    std::cout << BrNO3[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PRPE" << ": ";
    std::cout << std::setw(9);
    std::cout << PRPE[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MVKOO" << ": ";
    std::cout << std::setw(9);
    std::cout << MVKOO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "Cl2" << ": ";
    std::cout << std::setw(9);
    std::cout << Cl2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISOPND" << ": ";
    std::cout << std::setw(9);
    std::cout << ISOPND[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HOBr" << ": ";
    std::cout << std::setw(9);
    std::cout << HOBr[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "A3O2" << ": ";
    std::cout << std::setw(9);
    std::cout << A3O2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PROPNN" << ": ";
    std::cout << std::setw(9);
    std::cout << PROPNN[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "GLYX" << ": ";
    std::cout << std::setw(9);
    std::cout << GLYX[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MAOPO2" << ": ";
    std::cout << std::setw(9);
    std::cout << MAOPO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "CH4" << ": ";
    std::cout << std::setw(9);
    std::cout << CH4[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "GAOO" << ": ";
    std::cout << std::setw(9);
    std::cout << GAOO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "B3O2" << ": ";
    std::cout << std::setw(9);
    std::cout << B3O2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ACET" << ": ";
    std::cout << std::setw(9);
    std::cout << ACET[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MACRN" << ": ";
    std::cout << std::setw(9);
    std::cout << MACRN[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "CH2OO" << ": ";
    std::cout << std::setw(9);
    std::cout << CH2OO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MGLYOO" << ": ";
    std::cout << std::setw(9);
    std::cout << MGLYOO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "VRO2" << ": ";
    std::cout << std::setw(9);
    std::cout << VRO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MGLOO" << ": ";
    std::cout << std::setw(9);
    std::cout << MGLOO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MACROO" << ": ";
    std::cout << std::setw(9);
    std::cout << MACROO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PO2" << ": ";
    std::cout << std::setw(9);
    std::cout << PO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "CH3CHOO" << ": ";
    std::cout << std::setw(9);
    std::cout << CH3CHOO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MAN2" << ": ";
    std::cout << std::setw(9);
    std::cout << MAN2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISNOOA" << ": ";
    std::cout << std::setw(9);
    std::cout << ISNOOA[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "H2O2" << ": ";
    std::cout << std::setw(9);
    std::cout << H2O2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PRN1" << ": ";
    std::cout << std::setw(9);
    std::cout << PRN1[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ETO2" << ": ";
    std::cout << std::setw(9);
    std::cout << ETO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "KO2" << ": ";
    std::cout << std::setw(9);
    std::cout << KO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RCO3" << ": ";
    std::cout << std::setw(9);
    std::cout << RCO3[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HC5OO" << ": ";
    std::cout << std::setw(9);
    std::cout << HC5OO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "GLYC" << ": ";
    std::cout << std::setw(9);
    std::cout << GLYC[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ClNO3" << ": ";
    std::cout << std::setw(9);
    std::cout << ClNO3[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RIO2" << ": ";
    std::cout << std::setw(9);
    std::cout << RIO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "R4N1" << ": ";
    std::cout << std::setw(9);
    std::cout << R4N1[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HOCl" << ": ";
    std::cout << std::setw(9);
    std::cout << HOCl[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ATO2" << ": ";
    std::cout << std::setw(9);
    std::cout << ATO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HNO3" << ": ";
    std::cout << std::setw(9);
    std::cout << HNO3[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISN1" << ": ";
    std::cout << std::setw(9);
    std::cout << ISN1[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MAO3" << ": ";
    std::cout << std::setw(9);
    std::cout << MAO3[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MRO2" << ": ";
    std::cout << std::setw(9);
    std::cout << MRO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "INO2" << ": ";
    std::cout << std::setw(9);
    std::cout << INO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HAC" << ": ";
    std::cout << std::setw(9);
    std::cout << HAC[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HC5" << ": ";
    std::cout << std::setw(9);
    std::cout << HC5[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MGLY" << ": ";
    std::cout << std::setw(9);
    std::cout << MGLY[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISOPNBO2" << ": ";
    std::cout << std::setw(9);
    std::cout << ISOPNBO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISOPNDO2" << ": ";
    std::cout << std::setw(9);
    std::cout << ISOPNDO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "R4O2" << ": ";
    std::cout << std::setw(9);
    std::cout << R4O2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "R4N2" << ": ";
    std::cout << std::setw(9);
    std::cout << R4N2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "BrO" << ": ";
    std::cout << std::setw(9);
    std::cout << BrO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RCHO" << ": ";
    std::cout << std::setw(9);
    std::cout << RCHO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MEK" << ": ";
    std::cout << std::setw(9);
    std::cout << MEK[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ClO" << ": ";
    std::cout << std::setw(9);
    std::cout << ClO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MACR" << ": ";
    std::cout << std::setw(9);
    std::cout << MACR[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "SO2" << ": ";
    std::cout << std::setw(9);
    std::cout << SO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MVK" << ": ";
    std::cout << std::setw(9);
    std::cout << MVK[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ALD2" << ": ";
    std::cout << std::setw(9);
    std::cout << ALD2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MCO3" << ": ";
    std::cout << std::setw(9);
    std::cout << MCO3[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "CH2O" << ": ";
    std::cout << std::setw(9);
    std::cout << CH2O[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "H2O" << ": ";
    std::cout << std::setw(9);
    std::cout << H2O[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "Br" << ": ";
    std::cout << std::setw(9);
    std::cout << Br[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "NO" << ": ";
    std::cout << std::setw(9);
    std::cout << NO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "NO3" << ": ";
    std::cout << std::setw(9);
    std::cout << NO3[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "Cl" << ": ";
    std::cout << std::setw(9);
    std::cout << Cl[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "O" << ": ";
    std::cout << std::setw(9);
    std::cout << O[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "O1D" << ": ";
    std::cout << std::setw(9);
    std::cout << O1D[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "O3" << ": ";
    std::cout << std::setw(9);
    std::cout << O3[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HO2" << ": ";
    std::cout << std::setw(9);
    std::cout << HO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "NO2" << ": ";
    std::cout << std::setw(9);
    std::cout << NO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "OH" << ": ";
    std::cout << std::setw(9);
    std::cout << OH[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HBr" << ": ";
    std::cout << std::setw(9);
    std::cout << HBr[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HCl" << ": ";
    std::cout << std::setw(9);
    std::cout << HCl[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "CO" << ": ";
    std::cout << std::setw(9);
    std::cout << CO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MO2" << ": ";
    std::cout << std::setw(9);
    std::cout << MO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ACTA" << ": ";
    std::cout << std::setw(9);
    std::cout << ACTA[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "EOH" << ": ";
    std::cout << std::setw(9);
    std::cout << EOH[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "H2" << ": ";
    std::cout << std::setw(9);
    std::cout << H2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HCOOH" << ": ";
    std::cout << std::setw(9);
    std::cout << HCOOH[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MOH" << ": ";
    std::cout << std::setw(9);
    std::cout << MOH[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "N2" << ": ";
    std::cout << std::setw(9);
    std::cout << N2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "O2" << ": ";
    std::cout << std::setw(9);
    std::cout << O2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RCOOH" << ": ";
    std::cout << std::setw(9);
    std::cout << RCOOH[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

} /* End of Solution::Debug */

/* End of Structure.cpp */
