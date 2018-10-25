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
        liquidAerosol( ), solidAerosol( ),
        nVariables( NSPEC ), nAer( N_AER ), size_x( NX ), size_y( NY )
{
    /* Constructor */

} /* End of Solution::Solution */

Solution::~Solution()
{
    /* Destructor */

} /* End of Solution::~Solution */

void Solution::Clear( std::vector<std::vector<double> >& vector_2D )
{

    for ( unsigned int i = 0; i < vector_2D.size(); i++ ) {
        vector_2D[i].clear();
    }
    vector_2D.clear();

} /* End of Solution::Clear */

void Solution::SetShape( std::vector<std::vector<double> >& vector_2D, unsigned int n_x, unsigned int n_y, double value )
{
    
    Clear( vector_2D );

    /* Dimensions are transposed! */
    for ( unsigned int i = 0; i < n_y; i++ ) {
        vector_2D.push_back( std::vector<double>( n_x, value ) );
    }

} /* End of Solution::SetShape */

void Solution::SetToValue( std::vector<std::vector<double> >& vector_2D, double value )
{
    
    for ( unsigned int i = 0; i < vector_2D.size(); i++ ) {
        for ( unsigned int j = 0; j < vector_2D[0].size(); j++ ) {
            vector_2D[i][j] = (double) value;
        }
    }

} /* End of Solution::SetToValue */

void Solution::Print( std::vector<std::vector<double> >& vector_2D, unsigned int i_max, unsigned int j_max )
{
    
    for ( unsigned int i = 0; i < i_max; i++ ) {
        for ( unsigned int j = 0; j < j_max; j++ ) {
            std::cout << vector_2D[i][j];
        }
        std::cout << "" << std::endl;
    }

} /* End of Solution::Print */

void Solution::Initialize( char const *fileName, const double temperature, const double pressure, const double airDens, const double relHum, const double lat )
{

    const bool DBG = 0;

    std::vector<double> amb_Value(nVariables, 0.0);
    std::vector<std::vector<double> > aer_Value(nAer, std::vector<double>(2, 0.0));
    std::ifstream file;

    file.open( fileName );

    if ( file.is_open() )
    {
        std::string line;
        unsigned int i = 0;

        while ( ( std::getline( file, line ) ) && ( i < nVariables + nAer ) ) {
            if ( ( line != "\r" ) && ( line != "\n" ) && ( line[0] != '#' ) ) {
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
    else
    {
        std::string const currFunc("Structure::Initialize");
        std::cout << "ERROR: In " << currFunc << ": Can't read (" << fileName << ")" << std::endl;
    }

    /* Gaseous species */
    SetShape( CO2      , size_x , size_y, amb_Value[  0] * airDens );
    SetShape( PPN      , size_x , size_y, amb_Value[  1] * airDens );
    SetShape( BrNO2    , size_x , size_y, amb_Value[  2] * airDens );
    SetShape( IEPOX    , size_x , size_y, amb_Value[  3] * airDens );
    SetShape( PMNN     , size_x , size_y, amb_Value[  4] * airDens );
    SetShape( N2O      , size_x , size_y, amb_Value[  5] * airDens );
    SetShape( N        , size_x , size_y, amb_Value[  6] * airDens );
    SetShape( PAN      , size_x , size_y, amb_Value[  7] * airDens );
    SetShape( ALK4     , size_x , size_y, amb_Value[  8] * airDens );
    SetShape( MAP      , size_x , size_y, amb_Value[  9] * airDens );
    SetShape( MPN      , size_x , size_y, amb_Value[ 10] * airDens );
    SetShape( Cl2O2    , size_x , size_y, amb_Value[ 11] * airDens );
    SetShape( ETP      , size_x , size_y, amb_Value[ 12] * airDens );
    SetShape( HNO2     , size_x , size_y, amb_Value[ 13] * airDens );
    SetShape( C3H8     , size_x , size_y, amb_Value[ 14] * airDens );
    SetShape( RA3P     , size_x , size_y, amb_Value[ 15] * airDens );
    SetShape( RB3P     , size_x , size_y, amb_Value[ 16] * airDens );
    SetShape( OClO     , size_x , size_y, amb_Value[ 17] * airDens );
    SetShape( ClNO2    , size_x , size_y, amb_Value[ 18] * airDens );
    SetShape( ISOP     , size_x , size_y, amb_Value[ 19] * airDens );
    SetShape( HNO4     , size_x , size_y, amb_Value[ 20] * airDens );
    SetShape( MAOP     , size_x , size_y, amb_Value[ 21] * airDens );
    SetShape( MP       , size_x , size_y, amb_Value[ 22] * airDens );
    SetShape( ClOO     , size_x , size_y, amb_Value[ 23] * airDens );
    SetShape( RP       , size_x , size_y, amb_Value[ 24] * airDens );
    SetShape( BrCl     , size_x , size_y, amb_Value[ 25] * airDens );
    SetShape( PP       , size_x , size_y, amb_Value[ 26] * airDens );
    SetShape( PRPN     , size_x , size_y, amb_Value[ 27] * airDens );
    SetShape( SO4      , size_x , size_y, amb_Value[ 28] * airDens );
    SetShape( Br2      , size_x , size_y, amb_Value[ 29] * airDens );
    SetShape( ETHLN    , size_x , size_y, amb_Value[ 30] * airDens );
    SetShape( MVKN     , size_x , size_y, amb_Value[ 31] * airDens );
    SetShape( R4P      , size_x , size_y, amb_Value[ 32] * airDens );
    SetShape( C2H6     , size_x , size_y, amb_Value[ 33] * airDens );
    SetShape( RIP      , size_x , size_y, amb_Value[ 34] * airDens );
    SetShape( VRP      , size_x , size_y, amb_Value[ 35] * airDens );
    SetShape( ATOOH    , size_x , size_y, amb_Value[ 36] * airDens );
    SetShape( IAP      , size_x , size_y, amb_Value[ 37] * airDens );
    SetShape( DHMOB    , size_x , size_y, amb_Value[ 38] * airDens );
    SetShape( MOBA     , size_x , size_y, amb_Value[ 39] * airDens );
    SetShape( MRP      , size_x , size_y, amb_Value[ 40] * airDens );
    SetShape( N2O5     , size_x , size_y, amb_Value[ 41] * airDens );
    SetShape( ISNOHOO  , size_x , size_y, amb_Value[ 42] * airDens );
    SetShape( ISNP     , size_x , size_y, amb_Value[ 43] * airDens );
    SetShape( ISOPNB   , size_x , size_y, amb_Value[ 44] * airDens );
    SetShape( IEPOXOO  , size_x , size_y, amb_Value[ 45] * airDens );
    SetShape( MACRNO2  , size_x , size_y, amb_Value[ 46] * airDens );
    SetShape( ROH      , size_x , size_y, amb_Value[ 47] * airDens );
    SetShape( MOBAOO   , size_x , size_y, amb_Value[ 48] * airDens );
    SetShape( DIBOO    , size_x , size_y, amb_Value[ 49] * airDens );
    SetShape( PMN      , size_x , size_y, amb_Value[ 50] * airDens );
    SetShape( ISNOOB   , size_x , size_y, amb_Value[ 51] * airDens );
    SetShape( INPN     , size_x , size_y, amb_Value[ 52] * airDens );
    SetShape( H        , size_x , size_y, amb_Value[ 53] * airDens );
    SetShape( BrNO3    , size_x , size_y, amb_Value[ 54] * airDens );
    SetShape( PRPE     , size_x , size_y, amb_Value[ 55] * airDens );
    SetShape( MVKOO    , size_x , size_y, amb_Value[ 56] * airDens );
    SetShape( Cl2      , size_x , size_y, amb_Value[ 57] * airDens );
    SetShape( ISOPND   , size_x , size_y, amb_Value[ 58] * airDens );
    SetShape( HOBr     , size_x , size_y, amb_Value[ 59] * airDens );
    SetShape( A3O2     , size_x , size_y, amb_Value[ 60] * airDens );
    SetShape( PROPNN   , size_x , size_y, amb_Value[ 61] * airDens );
    SetShape( GLYX     , size_x , size_y, amb_Value[ 62] * airDens );
    SetShape( MAOPO2   , size_x , size_y, amb_Value[ 63] * airDens );
    SetShape( CH4      , size_x , size_y, amb_Value[ 64] * airDens );
    SetShape( GAOO     , size_x , size_y, amb_Value[ 65] * airDens );
    SetShape( B3O2     , size_x , size_y, amb_Value[ 66] * airDens );
    SetShape( ACET     , size_x , size_y, amb_Value[ 67] * airDens );
    SetShape( MACRN    , size_x , size_y, amb_Value[ 68] * airDens );
    SetShape( CH2OO    , size_x , size_y, amb_Value[ 69] * airDens );
    SetShape( MGLYOO   , size_x , size_y, amb_Value[ 70] * airDens );
    SetShape( VRO2     , size_x , size_y, amb_Value[ 71] * airDens );
    SetShape( MGLOO    , size_x , size_y, amb_Value[ 72] * airDens );
    SetShape( MACROO   , size_x , size_y, amb_Value[ 73] * airDens );
    SetShape( PO2      , size_x , size_y, amb_Value[ 74] * airDens );
    SetShape( CH3CHOO  , size_x , size_y, amb_Value[ 75] * airDens );
    SetShape( MAN2     , size_x , size_y, amb_Value[ 76] * airDens );
    SetShape( ISNOOA   , size_x , size_y, amb_Value[ 77] * airDens );
    SetShape( H2O2     , size_x , size_y, amb_Value[ 78] * airDens );
    SetShape( PRN1     , size_x , size_y, amb_Value[ 79] * airDens );
    SetShape( ETO2     , size_x , size_y, amb_Value[ 80] * airDens );
    SetShape( KO2      , size_x , size_y, amb_Value[ 81] * airDens );
    SetShape( RCO3     , size_x , size_y, amb_Value[ 82] * airDens );
    SetShape( HC5OO    , size_x , size_y, amb_Value[ 83] * airDens );
    SetShape( GLYC     , size_x , size_y, amb_Value[ 84] * airDens );
    SetShape( ClNO3    , size_x , size_y, amb_Value[ 85] * airDens );
    SetShape( RIO2     , size_x , size_y, amb_Value[ 86] * airDens );
    SetShape( R4N1     , size_x , size_y, amb_Value[ 87] * airDens );
    SetShape( HOCl     , size_x , size_y, amb_Value[ 88] * airDens );
    SetShape( ATO2     , size_x , size_y, amb_Value[ 89] * airDens );
    SetShape( HNO3     , size_x , size_y, amb_Value[ 90] * airDens );
    SetShape( ISN1     , size_x , size_y, amb_Value[ 91] * airDens );
    SetShape( MAO3     , size_x , size_y, amb_Value[ 92] * airDens );
    SetShape( MRO2     , size_x , size_y, amb_Value[ 93] * airDens );
    SetShape( INO2     , size_x , size_y, amb_Value[ 94] * airDens );
    SetShape( HAC      , size_x , size_y, amb_Value[ 95] * airDens );
    SetShape( HC5      , size_x , size_y, amb_Value[ 96] * airDens );
    SetShape( MGLY     , size_x , size_y, amb_Value[ 97] * airDens );
    SetShape( ISOPNBO2 , size_x , size_y, amb_Value[ 98] * airDens );
    SetShape( ISOPNDO2 , size_x , size_y, amb_Value[ 99] * airDens );
    SetShape( R4O2     , size_x , size_y, amb_Value[100] * airDens );
    SetShape( R4N2     , size_x , size_y, amb_Value[101] * airDens );
    SetShape( BrO      , size_x , size_y, amb_Value[102] * airDens );
    SetShape( RCHO     , size_x , size_y, amb_Value[103] * airDens );
    SetShape( MEK      , size_x , size_y, amb_Value[104] * airDens );
    SetShape( ClO      , size_x , size_y, amb_Value[105] * airDens );
    SetShape( MACR     , size_x , size_y, amb_Value[106] * airDens );
    SetShape( SO2      , size_x , size_y, amb_Value[107] * airDens );
    SetShape( MVK      , size_x , size_y, amb_Value[108] * airDens );
    SetShape( ALD2     , size_x , size_y, amb_Value[109] * airDens );
    SetShape( MCO3     , size_x , size_y, amb_Value[110] * airDens );
    SetShape( CH2O     , size_x , size_y, amb_Value[111] * airDens );
    SetShape( H2O      , size_x , size_y, amb_Value[112] * airDens );
    SetShape( Br       , size_x , size_y, amb_Value[113] * airDens );
    SetShape( NO       , size_x , size_y, amb_Value[114] * airDens );
    SetShape( NO3      , size_x , size_y, amb_Value[115] * airDens );
    SetShape( Cl       , size_x , size_y, amb_Value[116] * airDens );
    SetShape( O        , size_x , size_y, amb_Value[117] * airDens );
    SetShape( O1D      , size_x , size_y, amb_Value[118] * airDens );
    SetShape( O3       , size_x , size_y, amb_Value[119] * airDens );
    SetShape( HO2      , size_x , size_y, amb_Value[120] * airDens );
    SetShape( NO2      , size_x , size_y, amb_Value[121] * airDens );
    SetShape( OH       , size_x , size_y, amb_Value[122] * airDens );
    SetShape( HBr      , size_x , size_y, amb_Value[123] * airDens );
    SetShape( HCl      , size_x , size_y, amb_Value[124] * airDens );
    SetShape( CO       , size_x , size_y, amb_Value[125] * airDens );
    SetShape( MO2      , size_x , size_y, amb_Value[126] * airDens );
    SetShape( ACTA     , size_x , size_y, amb_Value[127] * airDens );
    SetShape( EOH      , size_x , size_y, amb_Value[128] * airDens );
    SetShape( H2       , size_x , size_y, amb_Value[129] * airDens );
    SetShape( HCOOH    , size_x , size_y, amb_Value[130] * airDens );
    SetShape( MOH      , size_x , size_y, amb_Value[131] * airDens );
    SetShape( N2       , size_x , size_y, amb_Value[132] * airDens );
    SetShape( O2       , size_x , size_y, amb_Value[133] * airDens );
    SetShape( RCOOH    , size_x , size_y, amb_Value[134] * airDens );

    /* Setting up ambient water vapor */
    for ( unsigned int i = 0; i < size_x; i++ ) {
        for ( unsigned int j = 0; j < size_y; j++ ) {
            H2O[j][i] = (relHum/((double) 100.0) * \
                             physFunc::pSat_H2Ol( temperature ) / ( physConst::kB * temperature )) / 1.00E+06;
            /* RH_w = x_H2O * P / Psat_H2Ol(T) = [H2O](#/cm3) * 1E6 * kB * T / Psat_H2Ol(T) */

        }
    }

    SetShape( NIT  , size_x, size_y, (double) 0.0 );
    SetShape( NAT  , size_x, size_y, (double) 0.0 );
    
    std::vector<double> stratData{ SO4[0][0], HNO3[0][0], HCl[0][0], HOCl[0][0], HBr[0][0], HOBr[0][0],
                                   H2O[0][0], ClNO3[0][0], BrNO3[0][0], NIT[0][0], NAT[0][0] };
   
    KHETI_SLA.assign( 11, 0.0 );
    AERFRAC.assign( 7, 0.0 );
    SOLIDFRAC.assign( 7, 0.0 );
    
    std::vector<double> RAD, RHO, KG, NDENS, SAD;
    RAD.assign( 2, 0.0 );
    RHO.assign( 2, 0.0 );
    KG.assign( 2, 0.0 );
    NDENS.assign( 2, 0.0 );
    SAD.assign( 2, 0.0 );

    STATE_PSC = STRAT_AER( temperature, pressure, airDens, lat, stratData, (2.0*XLIM)*(2.0*YLIM), KHETI_SLA, SOLIDFRAC, AERFRAC, RAD, RHO, KG, NDENS, SAD );
  
    /* Liquid/solid species */
    SetShape( SO4L , size_x, size_y, (double) AERFRAC[0]                          * stratData[0] );
    SetShape( SO4  , size_x, size_y, (double) ( 1.0 - AERFRAC[0] )                * stratData[0] );
    
    SetShape( H2OL , size_x, size_y, (double) AERFRAC[6]                          * stratData[6] );
    SetShape( H2OS , size_x, size_y, (double) SOLIDFRAC[6]                        * stratData[6] );
    SetShape( H2O  , size_x, size_y, (double) ( 1.0 - AERFRAC[6] - SOLIDFRAC[6] ) * stratData[6] );
    
    SetShape( HNO3L, size_x, size_y, (double) AERFRAC[1]                          * stratData[1] );
    SetShape( HNO3S, size_x, size_y, (double) SOLIDFRAC[1]                        * stratData[1] );
    SetShape( HNO3 , size_x, size_y, (double) ( 1.0 - AERFRAC[1] - SOLIDFRAC[1] ) * stratData[1] );
    
    SetShape( HClL , size_x, size_y, (double) AERFRAC[2]                          * stratData[2] );
    SetShape( HCl  , size_x, size_y, (double) ( 1.0 - AERFRAC[2] )                * stratData[2] );
    
    SetShape( HOClL, size_x, size_y, (double) AERFRAC[3]                          * stratData[3] );
    SetShape( HOCl , size_x, size_y, (double) ( 1.0 - AERFRAC[3] )                * stratData[3] );
    
    SetShape( HBrL , size_x, size_y, (double) AERFRAC[4]                          * stratData[4] );
    SetShape( HBr  , size_x, size_y, (double) ( 1.0 - AERFRAC[4] )                * stratData[4] );
    
    SetShape( HOBrL, size_x, size_y, (double) AERFRAC[5]                          * stratData[5] );
    SetShape( HOBr , size_x, size_y, (double) ( 1.0 - AERFRAC[5] )                * stratData[5] );
    
    SetShape( NIT  , size_x, size_y, (double) stratData[ 9] );
    SetShape( NAT  , size_x, size_y, (double) stratData[10] );

    
    /* Aerosols */
    /* Assume that soot particles are monodisperse */
    SetShape( sootDens , size_x, size_y, (double) aer_Value[  0][0] );
    SetShape( sootRadi , size_x, size_y, (double) aer_Value[  0][1] );
    SetShape( sootArea , size_x, size_y, (double) 4.0 / double(3.0) * physConst::PI * aer_Value[  0][0] * aer_Value[  0][1] * aer_Value[  0][1] * aer_Value[  0][1] );


    nBin_LA = std::floor( 1 + log( pow( (LA_R_HIG/LA_R_LOW), 3.0 ) ) / log( LA_VRAT ) );

    if ( DBG ) {
        std::cout << "\n DEBUG : LA_R_LOW  = " << LA_R_LOW * 1.00E+09 << " [nm]\n";
        std::cout << " DEBUG : LA_R_HIG  = " << LA_R_HIG * 1.00E+09 << " [nm]\n";
        std::cout << " DEBUG : LA_VRAT   = " << LA_VRAT             << " [-]\n";
        std::cout << " DEBUG : nBin_LA   = " << nBin_LA             << "\n";
        std::cout << " DEBUG : NDENS     = " << NDENS[1] * 1.00E-06 << " [#/cm^3]\n";
        std::cout << " DEBUG : REFF      = " << RAD[1] * 1.00E+06   << " [mum]\n";
    }
    
    std::vector<double> LA_rE( nBin_LA + 1, 0.0 ); /* Bin edges in m */
    std::vector<double> LA_rJ( nBin_LA    , 0.0 ); /* Bin center radius in m */
    
    for ( unsigned int iBin_LA = 0; iBin_LA < nBin_LA + 1; iBin_LA++ ) {
        LA_rE[iBin_LA] = LA_R_LOW * pow( LA_VRAT, iBin_LA / RealDouble(3.0) );              /* [m] */
    }
    for ( unsigned int iBin_LA = 0; iBin_LA < nBin_LA; iBin_LA++ ) {
        LA_rJ[iBin_LA] = 0.5 * ( LA_rE[iBin_LA] + LA_rE[iBin_LA+1] );                       /* [m] */
    }

    LA_nDens = NDENS[1] * 1.00E-06; /* [#/cm^3]      */
    LA_rEff  = RAD[1]   * 1.00E+09; /* [nm]          */
    LA_SAD   = SAD[1]   * 1.00E+06; /* [\mum^2/cm^3] */

    /* For a lognormal distribution:
     * r_eff = r_m * exp( 5/2 * ln(S)^2 ) 
     * A     = 4\pi N0 r_m^2 * exp ( 2 * ln(S)^2 ) 
     * A/r_eff^2 = 4\pi N0 * exp( - 1/2 * ln(S)^2 )
     *
     * ln(S) = sqrt(-2*ln(A/(4\pi r_eff^2 * N0)));
     * r_m = r_eff * exp( -5/2 * ln(S)^2 ); */
    
    const double sLA = sqrt( - 1.0 / (3.0) * log(SAD[1]/(4.0 * physConst::PI * RAD[1] * RAD[1] * NDENS[1] ) ) );
    const double rLA = std::max( RAD[1] * exp( - 2.5 * sLA * sLA ), 1.5 * LA_R_LOW ); 
    AIM::Grid_Aerosol LAAerosol( size_x, size_y, LA_rJ, LA_rE, LA_nDens, rLA, exp(sLA), "lognormal" );

    liquidAerosol = LAAerosol;

    if ( DBG ) {

        std::cout << "\n DEBUG : Comparing PDF's number density to exact number density :\n";
        std::cout << "         " << liquidAerosol.Moment( 0, 0, 0 ) << " v " << LA_nDens << " [#/cm^3]\n";
        std::cout << " DEBUG : Comparing PDF's surface area to Grainger's surface area:\n";
        std::cout << "         " << 4.0 * physConst::PI * liquidAerosol.Moment( 2, 0, 0 ) * 1.00E+12 << " v " << LA_SAD << " [mum^2/cm^3]\n";
        std::cout << " DEBUG : Comparing PDF's effective radius to Grainger's effective radius:\n";
        std::cout << "         " << liquidAerosol.getEffRadius( 0, 0 ) * 1.00E+09 << " v " << LA_rEff << " [nm]\n";

    }

    nBin_PA = std::floor( 1 + log( pow( (PA_R_HIG/PA_R_LOW), 3.0 ) ) / log( PA_VRAT ) );
    
    if ( DBG ) {
        std::cout << "\n DEBUG : PA_R_LOW  = " << PA_R_LOW * 1.00E+06 << " [mum]\n";
        std::cout << " DEBUG : PA_R_HIG  = " << PA_R_HIG * 1.00E+06 << " [mum]\n";
        std::cout << " DEBUG : PA_VRAT   = " << PA_VRAT             << " [-]\n";
        std::cout << " DEBUG : nBin_PA   = " << nBin_PA             << "\n";
        std::cout << " DEBUG : NDENS     = " << NDENS[0] * 1.00E-06 << " [#/cm^3]\n";
        std::cout << " DEBUG : REFF      = " << RAD[0] * 1.00E+06   << " [mum]\n";
    }

    std::vector<double> PA_rE( nBin_PA + 1, 0.0 ); /* Bin edges in m */
    std::vector<double> PA_rJ( nBin_PA    , 0.0 ); /* Bin center radius in m */
    
    for ( unsigned int iBin_PA = 0; iBin_PA < nBin_PA + 1; iBin_PA++ ) {
        PA_rE[iBin_PA] = PA_R_LOW * pow( PA_VRAT, iBin_PA / RealDouble(3.0) );              /* [m] */
    }
    for ( unsigned int iBin_PA = 0; iBin_PA < nBin_PA; iBin_PA++ ) {
        PA_rJ[iBin_PA] = 0.5 * ( PA_rE[iBin_PA] + PA_rE[iBin_PA+1] );                       /* [m] */
    }

    PA_nDens = NDENS[0] * 1.00E-06; /* [#/cm^3]      */
    PA_rEff  = RAD[0]   * 1.00E+09; /* [nm]          */
    PA_SAD   = SAD[0]   * 1.00E+06; /* [\mum^2/cm^3] */

    const double expsPA = 1.05;
    const double rPA = std::max( RAD[0] * exp( - 2.5 * log(expsPA) * log(expsPA) ), 1.5 * PA_R_LOW ); 
    AIM::Grid_Aerosol PAAerosol( size_x, size_y, PA_rJ, PA_rE, PA_nDens, rPA, expsPA, "lognormal" );

    solidAerosol = PAAerosol;

    if ( DBG ) {
        std::cout << "\n DEBUG : Comparing PDF's number density to exact number density :\n";
        std::cout << "         " << solidAerosol.Moment( 0, 0, 0 ) << " v " << PA_nDens << " [#/cm^3]\n";
        std::cout << " DEBUG : Comparing PDF's effective radius to actual effective radius:\n";
        std::cout << "         " << solidAerosol.getEffRadius( 0, 0 ) * 1.00E+09 << " v " << PA_rEff << " [nm]\n";
    }

    /* Tracers */
    SetShape( SO4T , size_x, size_y, (double) amb_Value[ 28] * airDens );

} /* End of Solution::Initialize */

void Solution::getData( double varArray[], double fixArray[], unsigned int i, unsigned int j )
{

    varArray[  0] = CO2[j][i];
    varArray[  1] = PPN[j][i];
    varArray[  2] = BrNO2[j][i];
    varArray[  3] = IEPOX[j][i];
    varArray[  4] = PMNN[j][i];
    varArray[  5] = N2O[j][i];
    varArray[  6] = N[j][i];
    varArray[  7] = PAN[j][i];
    varArray[  8] = ALK4[j][i];
    varArray[  9] = MAP[j][i];
    varArray[ 10] = MPN[j][i];
    varArray[ 11] = Cl2O2[j][i];
    varArray[ 12] = ETP[j][i];
    varArray[ 13] = HNO2[j][i];
    varArray[ 14] = C3H8[j][i];
    varArray[ 15] = RA3P[j][i];
    varArray[ 16] = RB3P[j][i];
    varArray[ 17] = OClO[j][i];
    varArray[ 18] = ClNO2[j][i];
    varArray[ 19] = ISOP[j][i];
    varArray[ 20] = HNO4[j][i];
    varArray[ 21] = MAOP[j][i];
    varArray[ 22] = MP[j][i];
    varArray[ 23] = ClOO[j][i];
    varArray[ 24] = RP[j][i];
    varArray[ 25] = BrCl[j][i];
    varArray[ 26] = PP[j][i];
    varArray[ 27] = PRPN[j][i];
    varArray[ 28] = SO4[j][i];
    varArray[ 29] = Br2[j][i];
    varArray[ 30] = ETHLN[j][i];
    varArray[ 31] = MVKN[j][i];
    varArray[ 32] = R4P[j][i];
    varArray[ 33] = C2H6[j][i];
    varArray[ 34] = RIP[j][i];
    varArray[ 35] = VRP[j][i];
    varArray[ 36] = ATOOH[j][i];
    varArray[ 37] = IAP[j][i];
    varArray[ 38] = DHMOB[j][i];
    varArray[ 39] = MOBA[j][i];
    varArray[ 40] = MRP[j][i];
    varArray[ 41] = N2O5[j][i];
    varArray[ 42] = ISNOHOO[j][i];
    varArray[ 43] = ISNP[j][i];
    varArray[ 44] = ISOPNB[j][i];
    varArray[ 45] = IEPOXOO[j][i];
    varArray[ 46] = MACRNO2[j][i];
    varArray[ 47] = ROH[j][i];
    varArray[ 48] = MOBAOO[j][i];
    varArray[ 49] = DIBOO[j][i];
    varArray[ 50] = PMN[j][i];
    varArray[ 51] = ISNOOB[j][i];
    varArray[ 52] = INPN[j][i];
    varArray[ 53] = H[j][i];
    varArray[ 54] = BrNO3[j][i];
    varArray[ 55] = PRPE[j][i];
    varArray[ 56] = MVKOO[j][i];
    varArray[ 57] = Cl2[j][i];
    varArray[ 58] = ISOPND[j][i];
    varArray[ 59] = HOBr[j][i];
    varArray[ 60] = A3O2[j][i];
    varArray[ 61] = PROPNN[j][i];
    varArray[ 62] = GLYX[j][i];
    varArray[ 63] = MAOPO2[j][i];
    varArray[ 64] = CH4[j][i];
    varArray[ 65] = GAOO[j][i];
    varArray[ 66] = B3O2[j][i];
    varArray[ 67] = ACET[j][i];
    varArray[ 68] = MACRN[j][i];
    varArray[ 69] = CH2OO[j][i];
    varArray[ 70] = MGLYOO[j][i];
    varArray[ 71] = VRO2[j][i];
    varArray[ 72] = MGLOO[j][i];
    varArray[ 73] = MACROO[j][i];
    varArray[ 74] = PO2[j][i];
    varArray[ 75] = CH3CHOO[j][i];
    varArray[ 76] = MAN2[j][i];
    varArray[ 77] = ISNOOA[j][i];
    varArray[ 78] = H2O2[j][i];
    varArray[ 79] = PRN1[j][i];
    varArray[ 80] = ETO2[j][i];
    varArray[ 81] = KO2[j][i];
    varArray[ 82] = RCO3[j][i];
    varArray[ 83] = HC5OO[j][i];
    varArray[ 84] = GLYC[j][i];
    varArray[ 85] = ClNO3[j][i];
    varArray[ 86] = RIO2[j][i];
    varArray[ 87] = R4N1[j][i];
    varArray[ 88] = HOCl[j][i];
    varArray[ 89] = ATO2[j][i];
    varArray[ 90] = HNO3[j][i];
    varArray[ 91] = ISN1[j][i];
    varArray[ 92] = MAO3[j][i];
    varArray[ 93] = MRO2[j][i];
    varArray[ 94] = INO2[j][i];
    varArray[ 95] = HAC[j][i];
    varArray[ 96] = HC5[j][i];
    varArray[ 97] = MGLY[j][i];
    varArray[ 98] = ISOPNBO2[j][i];
    varArray[ 99] = ISOPNDO2[j][i];
    varArray[100] = R4O2[j][i];
    varArray[101] = R4N2[j][i];
    varArray[102] = BrO[j][i];
    varArray[103] = RCHO[j][i];
    varArray[104] = MEK[j][i];
    varArray[105] = ClO[j][i];
    varArray[106] = MACR[j][i];
    varArray[107] = SO2[j][i];
    varArray[108] = MVK[j][i];
    varArray[109] = ALD2[j][i];
    varArray[110] = MCO3[j][i];
    varArray[111] = CH2O[j][i];
    varArray[112] = H2O[j][i];
    varArray[113] = Br[j][i];
    varArray[114] = NO[j][i];
    varArray[115] = NO3[j][i];
    varArray[116] = Cl[j][i];
    varArray[117] = O[j][i];
    varArray[118] = O1D[j][i];
    varArray[119] = O3[j][i];
    varArray[120] = HO2[j][i];
    varArray[121] = NO2[j][i];
    varArray[122] = OH[j][i];
    varArray[123] = HBr[j][i];
    varArray[124] = HCl[j][i];
    varArray[125] = CO[j][i];
    varArray[126] = MO2[j][i];
    fixArray[  0] = ACTA[j][i];
    fixArray[  1] = EOH[j][i];
    fixArray[  2] = H2[j][i];
    fixArray[  3] = HCOOH[j][i];
    fixArray[  4] = MOH[j][i];
    fixArray[  5] = N2[j][i];
    fixArray[  6] = O2[j][i];
    fixArray[  7] = RCOOH[j][i];

} /* End of Solution::getData */

void Solution::applyData( double varArray[], unsigned int i, unsigned int j )
{

    CO2[j][i]      = varArray[  0];
    PPN[j][i]      = varArray[  1];
    BrNO2[j][i]    = varArray[  2];
    IEPOX[j][i]    = varArray[  3];
    PMNN[j][i]     = varArray[  4];
    N2O[j][i]      = varArray[  5];
    N[j][i]        = varArray[  6];
    PAN[j][i]      = varArray[  7];
    ALK4[j][i]     = varArray[  8];
    MAP[j][i]      = varArray[  9];
    MPN[j][i]      = varArray[ 10];
    Cl2O2[j][i]    = varArray[ 11];
    ETP[j][i]      = varArray[ 12];
    HNO2[j][i]     = varArray[ 13];
    C3H8[j][i]     = varArray[ 14];
    RA3P[j][i]     = varArray[ 15];
    RB3P[j][i]     = varArray[ 16];
    OClO[j][i]     = varArray[ 17];
    ClNO2[j][i]    = varArray[ 18];
    ISOP[j][i]     = varArray[ 19];
    HNO4[j][i]     = varArray[ 20];
    MAOP[j][i]     = varArray[ 21];
    MP[j][i]       = varArray[ 22];
    ClOO[j][i]     = varArray[ 23];
    RP[j][i]       = varArray[ 24];
    BrCl[j][i]     = varArray[ 25];
    PP[j][i]       = varArray[ 26];
    PRPN[j][i]     = varArray[ 27];
    SO4[j][i]      = varArray[ 28];
    Br2[j][i]      = varArray[ 29];
    ETHLN[j][i]    = varArray[ 30];
    MVKN[j][i]     = varArray[ 31];
    R4P[j][i]      = varArray[ 32];
    C2H6[j][i]     = varArray[ 33];
    RIP[j][i]      = varArray[ 34];
    VRP[j][i]      = varArray[ 35];
    ATOOH[j][i]    = varArray[ 36];
    IAP[j][i]      = varArray[ 37];
    DHMOB[j][i]    = varArray[ 38];
    MOBA[j][i]     = varArray[ 39];
    MRP[j][i]      = varArray[ 40];
    N2O5[j][i]     = varArray[ 41];
    ISNOHOO[j][i]  = varArray[ 42];
    ISNP[j][i]     = varArray[ 43];
    ISOPNB[j][i]   = varArray[ 44];
    IEPOXOO[j][i]  = varArray[ 45];
    MACRNO2[j][i]  = varArray[ 46];
    ROH[j][i]      = varArray[ 47];
    MOBAOO[j][i]   = varArray[ 48];
    DIBOO[j][i]    = varArray[ 49];
    PMN[j][i]      = varArray[ 50];
    ISNOOB[j][i]   = varArray[ 51];
    INPN[j][i]     = varArray[ 52];
    H[j][i]        = varArray[ 53];
    BrNO3[j][i]    = varArray[ 54];
    PRPE[j][i]     = varArray[ 55];
    MVKOO[j][i]    = varArray[ 56];
    Cl2[j][i]      = varArray[ 57];
    ISOPND[j][i]   = varArray[ 58];
    HOBr[j][i]     = varArray[ 59];
    A3O2[j][i]     = varArray[ 60];
    PROPNN[j][i]   = varArray[ 61];
    GLYX[j][i]     = varArray[ 62];
    MAOPO2[j][i]   = varArray[ 63];
    CH4[j][i]      = varArray[ 64];
    GAOO[j][i]     = varArray[ 65];
    B3O2[j][i]     = varArray[ 66];
    ACET[j][i]     = varArray[ 67];
    MACRN[j][i]    = varArray[ 68];
    CH2OO[j][i]    = varArray[ 69];
    MGLYOO[j][i]   = varArray[ 70];
    VRO2[j][i]     = varArray[ 71];
    MGLOO[j][i]    = varArray[ 72];
    MACROO[j][i]   = varArray[ 73];
    PO2[j][i]      = varArray[ 74];
    CH3CHOO[j][i]  = varArray[ 75];
    MAN2[j][i]     = varArray[ 76];
    ISNOOA[j][i]   = varArray[ 77];
    H2O2[j][i]     = varArray[ 78];
    PRN1[j][i]     = varArray[ 79];
    ETO2[j][i]     = varArray[ 80];
    KO2[j][i]      = varArray[ 81];
    RCO3[j][i]     = varArray[ 82];
    HC5OO[j][i]    = varArray[ 83];
    GLYC[j][i]     = varArray[ 84];
    ClNO3[j][i]    = varArray[ 85];
    RIO2[j][i]     = varArray[ 86];
    R4N1[j][i]     = varArray[ 87];
    HOCl[j][i]     = varArray[ 88];
    ATO2[j][i]     = varArray[ 89];
    HNO3[j][i]     = varArray[ 90];
    ISN1[j][i]     = varArray[ 91];
    MAO3[j][i]     = varArray[ 92];
    MRO2[j][i]     = varArray[ 93];
    INO2[j][i]     = varArray[ 94];
    HAC[j][i]      = varArray[ 95];
    HC5[j][i]      = varArray[ 96];
    MGLY[j][i]     = varArray[ 97];
    ISOPNBO2[j][i] = varArray[ 98];
    ISOPNDO2[j][i] = varArray[ 99];
    R4O2[j][i]     = varArray[100];
    R4N2[j][i]     = varArray[101];
    BrO[j][i]      = varArray[102];
    RCHO[j][i]     = varArray[103];
    MEK[j][i]      = varArray[104];
    ClO[j][i]      = varArray[105];
    MACR[j][i]     = varArray[106];
    SO2[j][i]      = varArray[107];
    MVK[j][i]      = varArray[108];
    ALD2[j][i]     = varArray[109];
    MCO3[j][i]     = varArray[110];
    CH2O[j][i]     = varArray[111];
    H2O[j][i]      = varArray[112];
    Br[j][i]       = varArray[113];
    NO[j][i]       = varArray[114];
    NO3[j][i]      = varArray[115];
    Cl[j][i]       = varArray[116];
    O[j][i]        = varArray[117];
    O1D[j][i]      = varArray[118];
    O3[j][i]       = varArray[119];
    HO2[j][i]      = varArray[120];
    NO2[j][i]      = varArray[121];
    OH[j][i]       = varArray[122];
    HBr[j][i]      = varArray[123];
    HCl[j][i]      = varArray[124];
    CO[j][i]       = varArray[125];
    MO2[j][i]      = varArray[126];


} /* End of Solution::applyData */

void Solution::applyRing( double varArray[], double tempArray[], std::vector<std::vector<std::pair<unsigned int, unsigned int>>> mapRing2Mesh, unsigned int iRing )
{

    unsigned int i, j;
    unsigned int size = mapRing2Mesh[iRing].size();
    for ( unsigned int iList = 0; iList < size; iList++ ) {
        i = mapRing2Mesh[iRing][iList].first;
        j = mapRing2Mesh[iRing][iList].second;

        CO2[j][i]      *= varArray[  0] / tempArray[  0];
        PPN[j][i]      *= varArray[  1] / tempArray[  1];
        BrNO2[j][i]    *= varArray[  2] / tempArray[  2];
        IEPOX[j][i]    *= varArray[  3] / tempArray[  3];
        PMNN[j][i]     *= varArray[  4] / tempArray[  4];
        N2O[j][i]      *= varArray[  5] / tempArray[  5];
        N[j][i]        *= varArray[  6] / tempArray[  6];
        PAN[j][i]      *= varArray[  7] / tempArray[  7];
        ALK4[j][i]     *= varArray[  8] / tempArray[  8];
        MAP[j][i]      *= varArray[  9] / tempArray[  9];
        MPN[j][i]      *= varArray[ 10] / tempArray[ 10];
        Cl2O2[j][i]    *= varArray[ 11] / tempArray[ 11];
        ETP[j][i]      *= varArray[ 12] / tempArray[ 12];
        HNO2[j][i]     *= varArray[ 13] / tempArray[ 13];
        C3H8[j][i]     *= varArray[ 14] / tempArray[ 14];
        RA3P[j][i]     *= varArray[ 15] / tempArray[ 15];
        RB3P[j][i]     *= varArray[ 16] / tempArray[ 16];
        OClO[j][i]     *= varArray[ 17] / tempArray[ 17];
        ClNO2[j][i]    *= varArray[ 18] / tempArray[ 18];
        ISOP[j][i]     *= varArray[ 19] / tempArray[ 19];
        HNO4[j][i]     *= varArray[ 20] / tempArray[ 20];
        MAOP[j][i]     *= varArray[ 21] / tempArray[ 21];
        MP[j][i]       *= varArray[ 22] / tempArray[ 22];
        ClOO[j][i]     *= varArray[ 23] / tempArray[ 23];
        RP[j][i]       *= varArray[ 24] / tempArray[ 24];
        BrCl[j][i]     *= varArray[ 25] / tempArray[ 25];
        PP[j][i]       *= varArray[ 26] / tempArray[ 26];
        PRPN[j][i]     *= varArray[ 27] / tempArray[ 27];
        SO4[j][i]      *= varArray[ 28] / tempArray[ 28];
        Br2[j][i]      *= varArray[ 29] / tempArray[ 29];
        ETHLN[j][i]    *= varArray[ 30] / tempArray[ 30];
        MVKN[j][i]     *= varArray[ 31] / tempArray[ 31];
        R4P[j][i]      *= varArray[ 32] / tempArray[ 32];
        C2H6[j][i]     *= varArray[ 33] / tempArray[ 33];
        RIP[j][i]      *= varArray[ 34] / tempArray[ 34];
        VRP[j][i]      *= varArray[ 35] / tempArray[ 35];
        ATOOH[j][i]    *= varArray[ 36] / tempArray[ 36];
        IAP[j][i]      *= varArray[ 37] / tempArray[ 37];
        DHMOB[j][i]    *= varArray[ 38] / tempArray[ 38];
        MOBA[j][i]     *= varArray[ 39] / tempArray[ 39];
        MRP[j][i]      *= varArray[ 40] / tempArray[ 40];
        N2O5[j][i]     *= varArray[ 41] / tempArray[ 41];
        ISNOHOO[j][i]  *= varArray[ 42] / tempArray[ 42];
        ISNP[j][i]     *= varArray[ 43] / tempArray[ 43];
        ISOPNB[j][i]   *= varArray[ 44] / tempArray[ 44];
        IEPOXOO[j][i]  *= varArray[ 45] / tempArray[ 45];
        MACRNO2[j][i]  *= varArray[ 46] / tempArray[ 46];
        ROH[j][i]      *= varArray[ 47] / tempArray[ 47];
        MOBAOO[j][i]   *= varArray[ 48] / tempArray[ 48];
        DIBOO[j][i]    *= varArray[ 49] / tempArray[ 49];
        PMN[j][i]      *= varArray[ 50] / tempArray[ 50];
        ISNOOB[j][i]   *= varArray[ 51] / tempArray[ 51];
        INPN[j][i]     *= varArray[ 52] / tempArray[ 52];
        H[j][i]        *= varArray[ 53] / tempArray[ 53];
        BrNO3[j][i]    *= varArray[ 54] / tempArray[ 54];
        PRPE[j][i]     *= varArray[ 55] / tempArray[ 55];
        MVKOO[j][i]    *= varArray[ 56] / tempArray[ 56];
        Cl2[j][i]      *= varArray[ 57] / tempArray[ 57];
        ISOPND[j][i]   *= varArray[ 58] / tempArray[ 58];
        HOBr[j][i]     *= varArray[ 59] / tempArray[ 59];
        A3O2[j][i]     *= varArray[ 60] / tempArray[ 60];
        PROPNN[j][i]   *= varArray[ 61] / tempArray[ 61];
        GLYX[j][i]     *= varArray[ 62] / tempArray[ 62];
        MAOPO2[j][i]   *= varArray[ 63] / tempArray[ 63];
        CH4[j][i]      *= varArray[ 64] / tempArray[ 64];
        GAOO[j][i]     *= varArray[ 65] / tempArray[ 65];
        B3O2[j][i]     *= varArray[ 66] / tempArray[ 66];
        ACET[j][i]     *= varArray[ 67] / tempArray[ 67];
        MACRN[j][i]    *= varArray[ 68] / tempArray[ 68];
        CH2OO[j][i]    *= varArray[ 69] / tempArray[ 69];
        MGLYOO[j][i]   *= varArray[ 70] / tempArray[ 70];
        VRO2[j][i]     *= varArray[ 71] / tempArray[ 71];
        MGLOO[j][i]    *= varArray[ 72] / tempArray[ 72];
        MACROO[j][i]   *= varArray[ 73] / tempArray[ 73];
        PO2[j][i]      *= varArray[ 74] / tempArray[ 74];
        CH3CHOO[j][i]  *= varArray[ 75] / tempArray[ 75];
        MAN2[j][i]     *= varArray[ 76] / tempArray[ 76];
        ISNOOA[j][i]   *= varArray[ 77] / tempArray[ 77];
        H2O2[j][i]     *= varArray[ 78] / tempArray[ 78];
        PRN1[j][i]     *= varArray[ 79] / tempArray[ 79];
        ETO2[j][i]     *= varArray[ 80] / tempArray[ 80];
        KO2[j][i]      *= varArray[ 81] / tempArray[ 81];
        RCO3[j][i]     *= varArray[ 82] / tempArray[ 82];
        HC5OO[j][i]    *= varArray[ 83] / tempArray[ 83];
        GLYC[j][i]     *= varArray[ 84] / tempArray[ 84];
        ClNO3[j][i]    *= varArray[ 85] / tempArray[ 85];
        RIO2[j][i]     *= varArray[ 86] / tempArray[ 86];
        R4N1[j][i]     *= varArray[ 87] / tempArray[ 87];
        HOCl[j][i]     *= varArray[ 88] / tempArray[ 88];
        ATO2[j][i]     *= varArray[ 89] / tempArray[ 89];
        HNO3[j][i]     *= varArray[ 90] / tempArray[ 90];
        ISN1[j][i]     *= varArray[ 91] / tempArray[ 91];
        MAO3[j][i]     *= varArray[ 92] / tempArray[ 92];
        MRO2[j][i]     *= varArray[ 93] / tempArray[ 93];
        INO2[j][i]     *= varArray[ 94] / tempArray[ 94];
        HAC[j][i]      *= varArray[ 95] / tempArray[ 95];
        HC5[j][i]      *= varArray[ 96] / tempArray[ 96];
        MGLY[j][i]     *= varArray[ 97] / tempArray[ 97];
        ISOPNBO2[j][i] *= varArray[ 98] / tempArray[ 98];
        ISOPNDO2[j][i] *= varArray[ 99] / tempArray[ 99];
        R4O2[j][i]     *= varArray[100] / tempArray[100];
        R4N2[j][i]     *= varArray[101] / tempArray[101];
        BrO[j][i]      *= varArray[102] / tempArray[102];
        RCHO[j][i]     *= varArray[103] / tempArray[103];
        MEK[j][i]      *= varArray[104] / tempArray[104];
        ClO[j][i]      *= varArray[105] / tempArray[105];
        MACR[j][i]     *= varArray[106] / tempArray[106];
        SO2[j][i]      *= varArray[107] / tempArray[107];
        MVK[j][i]      *= varArray[108] / tempArray[108];
        ALD2[j][i]     *= varArray[109] / tempArray[109];
        MCO3[j][i]     *= varArray[110] / tempArray[110];
        CH2O[j][i]     *= varArray[111] / tempArray[111];
        H2O[j][i]      *= varArray[112] / tempArray[112];
        Br[j][i]       *= varArray[113] / tempArray[113];
        NO[j][i]       *= varArray[114] / tempArray[114];
        NO3[j][i]      *= varArray[115] / tempArray[115];
        Cl[j][i]       *= varArray[116] / tempArray[116];
        O[j][i]        *= varArray[117] / tempArray[117];
        O1D[j][i]      *= varArray[118] / tempArray[118];
        O3[j][i]       *= varArray[119] / tempArray[119];
        HO2[j][i]      *= varArray[120] / tempArray[120];
        NO2[j][i]      *= varArray[121] / tempArray[121];
        OH[j][i]       *= varArray[122] / tempArray[122];
        HBr[j][i]      *= varArray[123] / tempArray[123];
        HCl[j][i]      *= varArray[124] / tempArray[124];
        CO[j][i]       *= varArray[125] / tempArray[125];
        MO2[j][i]      *= varArray[126] / tempArray[126];
    
    }

} /* End of Solution::applyRing */

void Solution::applyAmbient( double varArray[], std::vector<std::vector<std::pair<unsigned int, unsigned int>>> mapRing2Mesh, unsigned int ambIndex )
{

    unsigned int i, j;
    unsigned int size = mapRing2Mesh[ambIndex].size();
    for ( unsigned int iList = 0; iList < size; iList++ ) {
        i = mapRing2Mesh[ambIndex][iList].first;
        j = mapRing2Mesh[ambIndex][iList].second;

        CO2[j][i]      = varArray[  0];
        PPN[j][i]      = varArray[  1];
        BrNO2[j][i]    = varArray[  2];
        IEPOX[j][i]    = varArray[  3];
        PMNN[j][i]     = varArray[  4];
        N2O[j][i]      = varArray[  5];
        N[j][i]        = varArray[  6];
        PAN[j][i]      = varArray[  7];
        ALK4[j][i]     = varArray[  8];
        MAP[j][i]      = varArray[  9];
        MPN[j][i]      = varArray[ 10];
        Cl2O2[j][i]    = varArray[ 11];
        ETP[j][i]      = varArray[ 12];
        HNO2[j][i]     = varArray[ 13];
        C3H8[j][i]     = varArray[ 14];
        RA3P[j][i]     = varArray[ 15];
        RB3P[j][i]     = varArray[ 16];
        OClO[j][i]     = varArray[ 17];
        ClNO2[j][i]    = varArray[ 18];
        ISOP[j][i]     = varArray[ 19];
        HNO4[j][i]     = varArray[ 20];
        MAOP[j][i]     = varArray[ 21];
        MP[j][i]       = varArray[ 22];
        ClOO[j][i]     = varArray[ 23];
        RP[j][i]       = varArray[ 24];
        BrCl[j][i]     = varArray[ 25];
        PP[j][i]       = varArray[ 26];
        PRPN[j][i]     = varArray[ 27];
        SO4[j][i]      = varArray[ 28];
        Br2[j][i]      = varArray[ 29];
        ETHLN[j][i]    = varArray[ 30];
        MVKN[j][i]     = varArray[ 31];
        R4P[j][i]      = varArray[ 32];
        C2H6[j][i]     = varArray[ 33];
        RIP[j][i]      = varArray[ 34];
        VRP[j][i]      = varArray[ 35];
        ATOOH[j][i]    = varArray[ 36];
        IAP[j][i]      = varArray[ 37];
        DHMOB[j][i]    = varArray[ 38];
        MOBA[j][i]     = varArray[ 39];
        MRP[j][i]      = varArray[ 40];
        N2O5[j][i]     = varArray[ 41];
        ISNOHOO[j][i]  = varArray[ 42];
        ISNP[j][i]     = varArray[ 43];
        ISOPNB[j][i]   = varArray[ 44];
        IEPOXOO[j][i]  = varArray[ 45];
        MACRNO2[j][i]  = varArray[ 46];
        ROH[j][i]      = varArray[ 47];
        MOBAOO[j][i]   = varArray[ 48];
        DIBOO[j][i]    = varArray[ 49];
        PMN[j][i]      = varArray[ 50];
        ISNOOB[j][i]   = varArray[ 51];
        INPN[j][i]     = varArray[ 52];
        H[j][i]        = varArray[ 53];
        BrNO3[j][i]    = varArray[ 54];
        PRPE[j][i]     = varArray[ 55];
        MVKOO[j][i]    = varArray[ 56];
        Cl2[j][i]      = varArray[ 57];
        ISOPND[j][i]   = varArray[ 58];
        HOBr[j][i]     = varArray[ 59];
        A3O2[j][i]     = varArray[ 60];
        PROPNN[j][i]   = varArray[ 61];
        GLYX[j][i]     = varArray[ 62];
        MAOPO2[j][i]   = varArray[ 63];
        CH4[j][i]      = varArray[ 64];
        GAOO[j][i]     = varArray[ 65];
        B3O2[j][i]     = varArray[ 66];
        ACET[j][i]     = varArray[ 67];
        MACRN[j][i]    = varArray[ 68];
        CH2OO[j][i]    = varArray[ 69];
        MGLYOO[j][i]   = varArray[ 70];
        VRO2[j][i]     = varArray[ 71];
        MGLOO[j][i]    = varArray[ 72];
        MACROO[j][i]   = varArray[ 73];
        PO2[j][i]      = varArray[ 74];
        CH3CHOO[j][i]  = varArray[ 75];
        MAN2[j][i]     = varArray[ 76];
        ISNOOA[j][i]   = varArray[ 77];
        H2O2[j][i]     = varArray[ 78];
        PRN1[j][i]     = varArray[ 79];
        ETO2[j][i]     = varArray[ 80];
        KO2[j][i]      = varArray[ 81];
        RCO3[j][i]     = varArray[ 82];
        HC5OO[j][i]    = varArray[ 83];
        GLYC[j][i]     = varArray[ 84];
        ClNO3[j][i]    = varArray[ 85];
        RIO2[j][i]     = varArray[ 86];
        R4N1[j][i]     = varArray[ 87];
        HOCl[j][i]     = varArray[ 88];
        ATO2[j][i]     = varArray[ 89];
        HNO3[j][i]     = varArray[ 90];
        ISN1[j][i]     = varArray[ 91];
        MAO3[j][i]     = varArray[ 92];
        MRO2[j][i]     = varArray[ 93];
        INO2[j][i]     = varArray[ 94];
        HAC[j][i]      = varArray[ 95];
        HC5[j][i]      = varArray[ 96];
        MGLY[j][i]     = varArray[ 97];
        ISOPNBO2[j][i] = varArray[ 98];
        ISOPNDO2[j][i] = varArray[ 99];
        R4O2[j][i]     = varArray[100];
        R4N2[j][i]     = varArray[101];
        BrO[j][i]      = varArray[102];
        RCHO[j][i]     = varArray[103];
        MEK[j][i]      = varArray[104];
        ClO[j][i]      = varArray[105];
        MACR[j][i]     = varArray[106];
        SO2[j][i]      = varArray[107];
        MVK[j][i]      = varArray[108];
        ALD2[j][i]     = varArray[109];
        MCO3[j][i]     = varArray[110];
        CH2O[j][i]     = varArray[111];
        H2O[j][i]      = varArray[112];
        Br[j][i]       = varArray[113];
        NO[j][i]       = varArray[114];
        NO3[j][i]      = varArray[115];
        Cl[j][i]       = varArray[116];
        O[j][i]        = varArray[117];
        O1D[j][i]      = varArray[118];
        O3[j][i]       = varArray[119];
        HO2[j][i]      = varArray[120];
        NO2[j][i]      = varArray[121];
        OH[j][i]       = varArray[122];
        HBr[j][i]      = varArray[123];
        HCl[j][i]      = varArray[124];
        CO[j][i]       = varArray[125];
        MO2[j][i]      = varArray[126];

    }

} /* End of Solution::applyAmbient */


void Solution::addEmission( const Emission &EI, const Aircraft &AC, std::vector<std::vector<std::pair<unsigned int, unsigned int>>> &map, std::vector<std::vector<double>> cellAreas, bool halfRing, const double temperature, bool set2Saturation, AIM::Aerosol &liqAer )
{

    unsigned int innerRing, nCell;
    unsigned int i, j;

    double E_CO2, E_H2O, E_NO, E_NO2, E_HNO2, E_SO2, E_CO, E_CH4, E_C2H6, E_PRPE, E_ALK4, E_CH2O, E_ALD2, E_GLYX, E_MGLY;
    E_CO2  = EI.getCO2()  / ( MW_CO2  * 1.0E+03 ) * AC.getFuelFlow()  / AC.getVFlight() * physConst::Na;
    /*     = [g/kg fuel]  / ( [kg/mol]* [g/kg]  ) * [kg fuel/s]       / [m/s]           * [molec/mol]
     *     = [molec/m] 
     */
    E_NO   = EI.getNO()   / ( MW_NO   * 1.0E+03 ) * AC.getFuelFlow() / AC.getVFlight() * physConst::Na;
    E_NO2  = EI.getNO2()  / ( MW_NO2  * 1.0E+03 ) * AC.getFuelFlow() / AC.getVFlight() * physConst::Na;
    E_HNO2 = EI.getHNO2() / ( MW_HNO2 * 1.0E+03 ) * AC.getFuelFlow() / AC.getVFlight() * physConst::Na;
    E_CO   = EI.getCO()   / ( MW_CO   * 1.0E+03 ) * AC.getFuelFlow() / AC.getVFlight() * physConst::Na;
    E_CH4  = EI.getCH4()  / ( MW_CH4  * 1.0E+03 ) * AC.getFuelFlow() / AC.getVFlight() * physConst::Na;
    E_C2H6 = EI.getC2H6() / ( MW_C2H6 * 1.0E+03 ) * AC.getFuelFlow() / AC.getVFlight() * physConst::Na;
    E_PRPE = EI.getPRPE() / ( MW_PRPE * 1.0E+03 ) * AC.getFuelFlow() / AC.getVFlight() * physConst::Na;
    E_ALK4 = EI.getALK4() / ( MW_ALK4 * 1.0E+03 ) * AC.getFuelFlow() / AC.getVFlight() * physConst::Na;
    E_CH2O = EI.getCH2O() / ( MW_CH2O * 1.0E+03 ) * AC.getFuelFlow() / AC.getVFlight() * physConst::Na;
    E_ALD2 = EI.getALD2() / ( MW_ALD2 * 1.0E+03 ) * AC.getFuelFlow() / AC.getVFlight() * physConst::Na;
    E_GLYX = EI.getGLYX() / ( MW_GLYX * 1.0E+03 ) * AC.getFuelFlow() / AC.getVFlight() * physConst::Na;
    E_MGLY = EI.getMGLY() / ( MW_MGLY * 1.0E+03 ) * AC.getFuelFlow() / AC.getVFlight() * physConst::Na;
    if ( !set2Saturation ) {
        E_H2O  = EI.getH2O()  / ( MW_H2O  * 1.0E+03 ) * AC.getFuelFlow() / AC.getVFlight() * physConst::Na;
    }
    E_SO2  = ( 1.0 - SO2TOSO4 ) * \
             EI.getSO2()  / ( MW_SO2  * 1.0E+03 ) * AC.getFuelFlow() / AC.getVFlight() * physConst::Na;

    if ( !halfRing ) {
        /* Full rings */
        innerRing = 0;
        nCell = map[innerRing].size();
        for ( unsigned int iList = 0; iList < map[innerRing].size(); iList++ ) {
            i = map[innerRing][iList].first;
            j = map[innerRing][iList].second;

            CO2[j][i]  += ( E_CO2  / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
            NO[j][i]   += ( E_NO   / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
            NO2[j][i]  += ( E_NO2  / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
            HNO2[j][i] += ( E_HNO2 / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
            CO[j][i]   += ( E_CO   / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
            CH4[j][i]  += ( E_CH4  / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
            C2H6[j][i] += ( E_C2H6 / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
            PRPE[j][i] += ( E_PRPE / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
            ALK4[j][i] += ( E_ALK4 / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
            CH2O[j][i] += ( E_CH2O / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
            ALD2[j][i] += ( E_ALD2 / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
            GLYX[j][i] += ( E_GLYX / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
            MGLY[j][i] += ( E_MGLY / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
            if ( set2Saturation ) {
                H2O[j][i] = physFunc::pSat_H2Os( temperature ) / ( physConst::kB * temperature * 1.00E+06 ); /* [molec / cm^3] */
            } else {
                H2O[j][i]+= ( E_H2O  / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
            }
            SO2[j][i]  += ( E_SO2  / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */

        }

        liquidAerosol.addPDF( liqAer, map[innerRing] );

    }
    else {
        /* Half rings */
        nCell = map[0].size() + map[1].size();
        for ( innerRing = 0; innerRing <= 1; innerRing++ ) {
            for ( unsigned int iList = 0; iList < map[innerRing].size(); iList++ ) {
                i = map[innerRing][iList].first;
                j = map[innerRing][iList].second;

                CO2[j][i]  += ( E_CO2  / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
                NO[j][i]   += ( E_NO   / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
                NO2[j][i]  += ( E_NO2  / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
                HNO2[j][i] += ( E_HNO2 / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
                CO[j][i]   += ( E_CO   / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
                CH4[j][i]  += ( E_CH4  / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
                C2H6[j][i] += ( E_C2H6 / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
                PRPE[j][i] += ( E_PRPE / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
                ALK4[j][i] += ( E_ALK4 / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
                CH2O[j][i] += ( E_CH2O / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
                ALD2[j][i] += ( E_ALD2 / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
                GLYX[j][i] += ( E_GLYX / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
                MGLY[j][i] += ( E_MGLY / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
                if ( set2Saturation ) {
                    H2O[j][i] = physFunc::pSat_H2Os( temperature ) / ( physConst::kB * temperature * 1.00E+06 ); /* [molec / cm^3] */
                } else {
                    H2O[j][i]+= ( E_H2O  / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */
                }
                SO2[j][i]  += ( E_SO2  / cellAreas[j][i] * 1.0E-06 / nCell ); /* [molec / cm^3] */

            }
        
            liquidAerosol.addPDF( liqAer, map[innerRing] );

        }

    }

} /* End of Solution::addEmission */

std::vector<double> Solution::getAmbient() const
{

    std::vector<double> ambVector( NSPEC, 0.0 );
    
    ambVector[  0] = CO2[0][0]      ;
    ambVector[  1] = PPN[0][0]      ;
    ambVector[  2] = BrNO2[0][0]    ;
    ambVector[  3] = IEPOX[0][0]    ;
    ambVector[  4] = PMNN[0][0]     ;
    ambVector[  5] = N2O[0][0]      ;
    ambVector[  6] = N[0][0]        ;
    ambVector[  7] = PAN[0][0]      ;
    ambVector[  8] = ALK4[0][0]     ;
    ambVector[  9] = MAP[0][0]      ;
    ambVector[ 10] = MPN[0][0]      ;
    ambVector[ 11] = Cl2O2[0][0]    ;
    ambVector[ 12] = ETP[0][0]      ;
    ambVector[ 13] = HNO2[0][0]     ;
    ambVector[ 14] = C3H8[0][0]     ;
    ambVector[ 15] = RA3P[0][0]     ;
    ambVector[ 16] = RB3P[0][0]     ;
    ambVector[ 17] = OClO[0][0]     ;
    ambVector[ 18] = ClNO2[0][0]    ;
    ambVector[ 19] = ISOP[0][0]     ;
    ambVector[ 20] = HNO4[0][0]     ;
    ambVector[ 21] = MAOP[0][0]     ;
    ambVector[ 22] = MP[0][0]       ;
    ambVector[ 23] = ClOO[0][0]     ;
    ambVector[ 24] = RP[0][0]       ;
    ambVector[ 25] = BrCl[0][0]     ;
    ambVector[ 26] = PP[0][0]       ;
    ambVector[ 27] = PRPN[0][0]     ;
    ambVector[ 28] = SO4[0][0]      ;
    ambVector[ 29] = Br2[0][0]      ;
    ambVector[ 30] = ETHLN[0][0]    ;
    ambVector[ 31] = MVKN[0][0]     ;
    ambVector[ 32] = R4P[0][0]      ;
    ambVector[ 33] = C2H6[0][0]     ;
    ambVector[ 34] = RIP[0][0]      ;
    ambVector[ 35] = VRP[0][0]      ;
    ambVector[ 36] = ATOOH[0][0]    ;
    ambVector[ 37] = IAP[0][0]      ;
    ambVector[ 38] = DHMOB[0][0]    ;
    ambVector[ 39] = MOBA[0][0]     ;
    ambVector[ 40] = MRP[0][0]      ;
    ambVector[ 41] = N2O5[0][0]     ;
    ambVector[ 42] = ISNOHOO[0][0]  ;
    ambVector[ 43] = ISNP[0][0]     ;
    ambVector[ 44] = ISOPNB[0][0]   ;
    ambVector[ 45] = IEPOXOO[0][0]  ;
    ambVector[ 46] = MACRNO2[0][0]  ;
    ambVector[ 47] = ROH[0][0]      ;
    ambVector[ 48] = MOBAOO[0][0]   ;
    ambVector[ 49] = DIBOO[0][0]    ;
    ambVector[ 50] = PMN[0][0]      ;
    ambVector[ 51] = ISNOOB[0][0]   ;
    ambVector[ 52] = INPN[0][0]     ;
    ambVector[ 53] = H[0][0]        ;
    ambVector[ 54] = BrNO3[0][0]    ;
    ambVector[ 55] = PRPE[0][0]     ;
    ambVector[ 56] = MVKOO[0][0]    ;
    ambVector[ 57] = Cl2[0][0]      ;
    ambVector[ 58] = ISOPND[0][0]   ;
    ambVector[ 59] = HOBr[0][0]     ;
    ambVector[ 60] = A3O2[0][0]     ;
    ambVector[ 61] = PROPNN[0][0]   ;
    ambVector[ 62] = GLYX[0][0]     ;
    ambVector[ 63] = MAOPO2[0][0]   ;
    ambVector[ 64] = CH4[0][0]      ;
    ambVector[ 65] = GAOO[0][0]     ;
    ambVector[ 66] = B3O2[0][0]     ;
    ambVector[ 67] = ACET[0][0]     ;
    ambVector[ 68] = MACRN[0][0]    ;
    ambVector[ 69] = CH2OO[0][0]    ;
    ambVector[ 70] = MGLYOO[0][0]   ;
    ambVector[ 71] = VRO2[0][0]     ;
    ambVector[ 72] = MGLOO[0][0]    ;
    ambVector[ 73] = MACROO[0][0]   ;
    ambVector[ 74] = PO2[0][0]      ;
    ambVector[ 75] = CH3CHOO[0][0]  ;
    ambVector[ 76] = MAN2[0][0]     ;
    ambVector[ 77] = ISNOOA[0][0]   ;
    ambVector[ 78] = H2O2[0][0]     ;
    ambVector[ 79] = PRN1[0][0]     ;
    ambVector[ 80] = ETO2[0][0]     ;
    ambVector[ 81] = KO2[0][0]      ;
    ambVector[ 82] = RCO3[0][0]     ;
    ambVector[ 83] = HC5OO[0][0]    ;
    ambVector[ 84] = GLYC[0][0]     ;
    ambVector[ 85] = ClNO3[0][0]    ;
    ambVector[ 86] = RIO2[0][0]     ;
    ambVector[ 87] = R4N1[0][0]     ;
    ambVector[ 88] = HOCl[0][0]     ;
    ambVector[ 89] = ATO2[0][0]     ;
    ambVector[ 90] = HNO3[0][0]     ;
    ambVector[ 91] = ISN1[0][0]     ;
    ambVector[ 92] = MAO3[0][0]     ;
    ambVector[ 93] = MRO2[0][0]     ;
    ambVector[ 94] = INO2[0][0]     ;
    ambVector[ 95] = HAC[0][0]      ;
    ambVector[ 96] = HC5[0][0]      ;
    ambVector[ 97] = MGLY[0][0]     ;
    ambVector[ 98] = ISOPNBO2[0][0] ;
    ambVector[ 99] = ISOPNDO2[0][0] ;
    ambVector[100] = R4O2[0][0]     ;
    ambVector[101] = R4N2[0][0]     ;
    ambVector[102] = BrO[0][0]      ;
    ambVector[103] = RCHO[0][0]     ;
    ambVector[104] = MEK[0][0]      ;
    ambVector[105] = ClO[0][0]      ;
    ambVector[106] = MACR[0][0]     ;
    ambVector[107] = SO2[0][0]      ;
    ambVector[108] = MVK[0][0]      ;
    ambVector[109] = ALD2[0][0]     ;
    ambVector[110] = MCO3[0][0]     ;
    ambVector[111] = CH2O[0][0]     ;
    ambVector[112] = H2O[0][0]      ;
    ambVector[113] = Br[0][0]       ;
    ambVector[114] = NO[0][0]       ;
    ambVector[115] = NO3[0][0]      ;
    ambVector[116] = Cl[0][0]       ;
    ambVector[117] = O[0][0]        ;
    ambVector[118] = O1D[0][0]      ;
    ambVector[119] = O3[0][0]       ;
    ambVector[120] = HO2[0][0]      ;
    ambVector[121] = NO2[0][0]      ;
    ambVector[122] = OH[0][0]       ;
    ambVector[123] = HBr[0][0]      ;
    ambVector[124] = HCl[0][0]      ;
    ambVector[125] = CO[0][0]       ;
    ambVector[126] = MO2[0][0]      ;
    ambVector[127] = ACTA[0][0]     ;
    ambVector[128] = EOH[0][0]      ;
    ambVector[129] = H2[0][0]       ;
    ambVector[130] = HCOOH[0][0]    ;
    ambVector[131] = MOH[0][0]      ;
    ambVector[132] = N2[0][0]       ;
    ambVector[133] = O2[0][0]       ;
    ambVector[134] = RCOOH[0][0]    ;

    return ambVector;

} /* End of Solution::getAmbient */

std::vector<double> Solution::getLiqSpecies( ) const
{

    std::vector<double> liqAerVector( 9, 0.0 );
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

std::vector<std::vector<double> > Solution::getAerosol( ) const
{

    std::vector<std::vector<double> > aerVector( nAer, std::vector<double>( 3, 0.0 ) );
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

std::vector<double> Solution::getAerosolDens( ) const
{

    std::vector<double> aerVector( nAer, 0.0 );
    aerVector[  0] = sootDens[0][0];
    aerVector[  1] = PA_nDens;
    aerVector[  2] = LA_nDens;

    return aerVector;

} /* End of Solution::getAerosolDens */

std::vector<double> Solution::getAerosolRadi( ) const
{

    std::vector<double> aerVector( nAer, 0.0 );
    aerVector[  0] = sootRadi[0][0];
    aerVector[  1] = PA_rEff;
    aerVector[  2] = LA_rEff;

    return aerVector;

} /* End of Solution::getAerosolRadi */

std::vector<double> Solution::getAerosolArea( ) const
{

    std::vector<double> aerVector( nAer, 0.0 );
    aerVector[  0] = sootArea[0][0];
    aerVector[  1] = PA_SAD;
    aerVector[  2] = LA_SAD;

    return aerVector;

} /* End of Solution::getAerosolArea */

unsigned int Solution::getNx() const
{

    return size_x;

} /* End of Solution::getNx */

unsigned int Solution::getNy() const
{

    return size_y;

} /* End of Solution::getNy */

void Solution::Debug( const double airDens )
{
    unsigned int iNx, jNy;
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
    std::cout << CO2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PPN" << ": ";
    std::cout << std::setw(9);
    std::cout << PPN[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "BrNO2" << ": ";
    std::cout << std::setw(9);
    std::cout << BrNO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "IEPOX" << ": ";
    std::cout << std::setw(9);
    std::cout << IEPOX[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PMNN" << ": ";
    std::cout << std::setw(9);
    std::cout << PMNN[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "N2O" << ": ";
    std::cout << std::setw(9);
    std::cout << N2O[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "N" << ": ";
    std::cout << std::setw(9);
    std::cout << N[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PAN" << ": ";
    std::cout << std::setw(9);
    std::cout << PAN[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ALK4" << ": ";
    std::cout << std::setw(9);
    std::cout << ALK4[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MAP" << ": ";
    std::cout << std::setw(9);
    std::cout << MAP[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MPN" << ": ";
    std::cout << std::setw(9);
    std::cout << MPN[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "Cl2O2" << ": ";
    std::cout << std::setw(9);
    std::cout << Cl2O2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ETP" << ": ";
    std::cout << std::setw(9);
    std::cout << ETP[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HNO2" << ": ";
    std::cout << std::setw(9);
    std::cout << HNO2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "C3H8" << ": ";
    std::cout << std::setw(9);
    std::cout << C3H8[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RA3P" << ": ";
    std::cout << std::setw(9);
    std::cout << RA3P[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RB3P" << ": ";
    std::cout << std::setw(9);
    std::cout << RB3P[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "OClO" << ": ";
    std::cout << std::setw(9);
    std::cout << OClO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ClNO2" << ": ";
    std::cout << std::setw(9);
    std::cout << ClNO2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISOP" << ": ";
    std::cout << std::setw(9);
    std::cout << ISOP[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HNO4" << ": ";
    std::cout << std::setw(9);
    std::cout << HNO4[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MAOP" << ": ";
    std::cout << std::setw(9);
    std::cout << MAOP[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MP" << ": ";
    std::cout << std::setw(9);
    std::cout << MP[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ClOO" << ": ";
    std::cout << std::setw(9);
    std::cout << ClOO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RP" << ": ";
    std::cout << std::setw(9);
    std::cout << RP[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "BrCl" << ": ";
    std::cout << std::setw(9);
    std::cout << BrCl[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PP" << ": ";
    std::cout << std::setw(9);
    std::cout << PP[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PRPN" << ": ";
    std::cout << std::setw(9);
    std::cout << PRPN[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "SO4" << ": ";
    std::cout << std::setw(9);
    std::cout << SO4[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "Br2" << ": ";
    std::cout << std::setw(9);
    std::cout << Br2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ETHLN" << ": ";
    std::cout << std::setw(9);
    std::cout << ETHLN[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MVKN" << ": ";
    std::cout << std::setw(9);
    std::cout << MVKN[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "R4P" << ": ";
    std::cout << std::setw(9);
    std::cout << R4P[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "C2H6" << ": ";
    std::cout << std::setw(9);
    std::cout << C2H6[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RIP" << ": ";
    std::cout << std::setw(9);
    std::cout << RIP[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "VRP" << ": ";
    std::cout << std::setw(9);
    std::cout << VRP[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ATOOH" << ": ";
    std::cout << std::setw(9);
    std::cout << ATOOH[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "IAP" << ": ";
    std::cout << std::setw(9);
    std::cout << IAP[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "DHMOB" << ": ";
    std::cout << std::setw(9);
    std::cout << DHMOB[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MOBA" << ": ";
    std::cout << std::setw(9);
    std::cout << MOBA[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MRP" << ": ";
    std::cout << std::setw(9);
    std::cout << MRP[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "N2O5" << ": ";
    std::cout << std::setw(9);
    std::cout << N2O5[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISNOHOO" << ": ";
    std::cout << std::setw(9);
    std::cout << ISNOHOO[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISNP" << ": ";
    std::cout << std::setw(9);
    std::cout << ISNP[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISOPNB" << ": ";
    std::cout << std::setw(9);
    std::cout << ISOPNB[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "IEPOXOO" << ": ";
    std::cout << std::setw(9);
    std::cout << IEPOXOO[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MACRNO2" << ": ";
    std::cout << std::setw(9);
    std::cout << MACRNO2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ROH" << ": ";
    std::cout << std::setw(9);
    std::cout << ROH[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MOBAOO" << ": ";
    std::cout << std::setw(9);
    std::cout << MOBAOO[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "DIBOO" << ": ";
    std::cout << std::setw(9);
    std::cout << DIBOO[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PMN" << ": ";
    std::cout << std::setw(9);
    std::cout << PMN[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISNOOB" << ": ";
    std::cout << std::setw(9);
    std::cout << ISNOOB[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "INPN" << ": ";
    std::cout << std::setw(9);
    std::cout << INPN[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "H" << ": ";
    std::cout << std::setw(9);
    std::cout << H[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "BrNO3" << ": ";
    std::cout << std::setw(9);
    std::cout << BrNO3[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PRPE" << ": ";
    std::cout << std::setw(9);
    std::cout << PRPE[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MVKOO" << ": ";
    std::cout << std::setw(9);
    std::cout << MVKOO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "Cl2" << ": ";
    std::cout << std::setw(9);
    std::cout << Cl2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISOPND" << ": ";
    std::cout << std::setw(9);
    std::cout << ISOPND[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HOBr" << ": ";
    std::cout << std::setw(9);
    std::cout << HOBr[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "A3O2" << ": ";
    std::cout << std::setw(9);
    std::cout << A3O2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PROPNN" << ": ";
    std::cout << std::setw(9);
    std::cout << PROPNN[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "GLYX" << ": ";
    std::cout << std::setw(9);
    std::cout << GLYX[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MAOPO2" << ": ";
    std::cout << std::setw(9);
    std::cout << MAOPO2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "CH4" << ": ";
    std::cout << std::setw(9);
    std::cout << CH4[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "GAOO" << ": ";
    std::cout << std::setw(9);
    std::cout << GAOO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "B3O2" << ": ";
    std::cout << std::setw(9);
    std::cout << B3O2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ACET" << ": ";
    std::cout << std::setw(9);
    std::cout << ACET[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MACRN" << ": ";
    std::cout << std::setw(9);
    std::cout << MACRN[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "CH2OO" << ": ";
    std::cout << std::setw(9);
    std::cout << CH2OO[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MGLYOO" << ": ";
    std::cout << std::setw(9);
    std::cout << MGLYOO[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "VRO2" << ": ";
    std::cout << std::setw(9);
    std::cout << VRO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MGLOO" << ": ";
    std::cout << std::setw(9);
    std::cout << MGLOO[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MACROO" << ": ";
    std::cout << std::setw(9);
    std::cout << MACROO[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PO2" << ": ";
    std::cout << std::setw(9);
    std::cout << PO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "CH3CHOO" << ": ";
    std::cout << std::setw(9);
    std::cout << CH3CHOO[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MAN2" << ": ";
    std::cout << std::setw(9);
    std::cout << MAN2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISNOOA" << ": ";
    std::cout << std::setw(9);
    std::cout << ISNOOA[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "H2O2" << ": ";
    std::cout << std::setw(9);
    std::cout << H2O2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PRN1" << ": ";
    std::cout << std::setw(9);
    std::cout << PRN1[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ETO2" << ": ";
    std::cout << std::setw(9);
    std::cout << ETO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "KO2" << ": ";
    std::cout << std::setw(9);
    std::cout << KO2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RCO3" << ": ";
    std::cout << std::setw(9);
    std::cout << RCO3[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HC5OO" << ": ";
    std::cout << std::setw(9);
    std::cout << HC5OO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "GLYC" << ": ";
    std::cout << std::setw(9);
    std::cout << GLYC[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ClNO3" << ": ";
    std::cout << std::setw(9);
    std::cout << ClNO3[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RIO2" << ": ";
    std::cout << std::setw(9);
    std::cout << RIO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "R4N1" << ": ";
    std::cout << std::setw(9);
    std::cout << R4N1[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HOCl" << ": ";
    std::cout << std::setw(9);
    std::cout << HOCl[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ATO2" << ": ";
    std::cout << std::setw(9);
    std::cout << ATO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HNO3" << ": ";
    std::cout << std::setw(9);
    std::cout << HNO3[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISN1" << ": ";
    std::cout << std::setw(9);
    std::cout << ISN1[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MAO3" << ": ";
    std::cout << std::setw(9);
    std::cout << MAO3[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MRO2" << ": ";
    std::cout << std::setw(9);
    std::cout << MRO2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "INO2" << ": ";
    std::cout << std::setw(9);
    std::cout << INO2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HAC" << ": ";
    std::cout << std::setw(9);
    std::cout << HAC[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HC5" << ": ";
    std::cout << std::setw(9);
    std::cout << HC5[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MGLY" << ": ";
    std::cout << std::setw(9);
    std::cout << MGLY[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISOPNBO2" << ": ";
    std::cout << std::setw(9);
    std::cout << ISOPNBO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISOPNDO2" << ": ";
    std::cout << std::setw(9);
    std::cout << ISOPNDO2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "R4O2" << ": ";
    std::cout << std::setw(9);
    std::cout << R4O2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "R4N2" << ": ";
    std::cout << std::setw(9);
    std::cout << R4N2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "BrO" << ": ";
    std::cout << std::setw(9);
    std::cout << BrO[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RCHO" << ": ";
    std::cout << std::setw(9);
    std::cout << RCHO[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MEK" << ": ";
    std::cout << std::setw(9);
    std::cout << MEK[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ClO" << ": ";
    std::cout << std::setw(9);
    std::cout << ClO[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MACR" << ": ";
    std::cout << std::setw(9);
    std::cout << MACR[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "SO2" << ": ";
    std::cout << std::setw(9);
    std::cout << SO2[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MVK" << ": ";
    std::cout << std::setw(9);
    std::cout << MVK[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ALD2" << ": ";
    std::cout << std::setw(9);
    std::cout << ALD2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MCO3" << ": ";
    std::cout << std::setw(9);
    std::cout << MCO3[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "CH2O" << ": ";
    std::cout << std::setw(9);
    std::cout << CH2O[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "H2O" << ": ";
    std::cout << std::setw(9);
    std::cout << H2O[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "Br" << ": ";
    std::cout << std::setw(9);
    std::cout << Br[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "NO" << ": ";
    std::cout << std::setw(9);
    std::cout << NO[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "NO3" << ": ";
    std::cout << std::setw(9);
    std::cout << NO3[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "Cl" << ": ";
    std::cout << std::setw(9);
    std::cout << Cl[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "O" << ": ";
    std::cout << std::setw(9);
    std::cout << O[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "O1D" << ": ";
    std::cout << std::setw(9);
    std::cout << O1D[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "O3" << ": ";
    std::cout << std::setw(9);
    std::cout << O3[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HO2" << ": ";
    std::cout << std::setw(9);
    std::cout << HO2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "NO2" << ": ";
    std::cout << std::setw(9);
    std::cout << NO2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "OH" << ": ";
    std::cout << std::setw(9);
    std::cout << OH[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HBr" << ": ";
    std::cout << std::setw(9);
    std::cout << HBr[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HCl" << ": ";
    std::cout << std::setw(9);
    std::cout << HCl[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "CO" << ": ";
    std::cout << std::setw(9);
    std::cout << CO[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MO2" << ": ";
    std::cout << std::setw(9);
    std::cout << MO2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";//" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ACTA" << ": ";
    std::cout << std::setw(9);
    std::cout << ACTA[jNy][iNx]/airDens*1.0E+09 << " [ppb],";// << std::endl;
    std::cout << std::setw(9);
    std::cout << "EOH" << ": ";
    std::cout << std::setw(9);
    std::cout << EOH[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "H2" << ": ";
    std::cout << std::setw(9);
    std::cout << H2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";// << std::endl;
    std::cout << std::setw(9);
    std::cout << "HCOOH" << ": ";
    std::cout << std::setw(9);
    std::cout << HCOOH[jNy][iNx]/airDens*1.0E+09 << " [ppb],";// << std::endl;
    std::cout << std::setw(9);
    std::cout << "MOH" << ": ";
    std::cout << std::setw(9);
    std::cout << MOH[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "N2" << ": ";
    std::cout << std::setw(9);
    std::cout << N2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";// << std::endl;
    std::cout << std::setw(9);
    std::cout << "O2" << ": ";
    std::cout << std::setw(9);
    std::cout << O2[jNy][iNx]/airDens*1.0E+09 << " [ppb],";// << std::endl;
    std::cout << std::setw(9);
    std::cout << "RCOOH" << ": ";
    std::cout << std::setw(9);
    std::cout << RCOOH[jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

} /* End of Solution::Debug */

/* End of Structure.cpp */
