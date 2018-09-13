#include "Structure.hpp"

double pSat_H2Ol( double T );
//double pSat_H2Os( double T );

Solution::Solution( const int nVar, const int n_x, const int n_y ) : \
        nVariables( nVar ), size_x( n_x ), size_y( n_y )
{
    /* Constructor */
}

Solution::~Solution()
{
    /* Destructor */
}

void Solution::Clear ( std::vector<std::vector<double> >& vector_2D )
{
    for ( unsigned int i = 0; i < vector_2D.size(); i++ ) {
        vector_2D[i].clear();
    }
    vector_2D.clear();

}

void Solution::SetShape ( std::vector<std::vector<double> >& vector_2D, unsigned int n_x, unsigned int n_y, double value )
{
    Clear( vector_2D );

    /* Dimensions are transposed! */
    for ( unsigned int i = 0; i < n_y; i++ ) {
        vector_2D.push_back( std::vector<double>( n_x, value ) );
    }

}

void Solution::SetToValue( std::vector<std::vector<double> >& vector_2D, double value )
{
    for ( unsigned int i = 0; i < vector_2D.size(); i++ ) {
        for ( unsigned int j = 0; j < vector_2D[0].size(); j++ ) {
            vector_2D[i][j] = (double) value;
        }
    }

}

void Solution::Print( std::vector<std::vector<double> >& vector_2D, unsigned int i_max, unsigned int j_max )
{
    for ( unsigned int i = 0; i < i_max; i++ ) {
        for ( unsigned int j = 0; j < j_max; j++ ) {
            std::cout << vector_2D[i][j];
        }
        std::cout << "" << std::endl;
    }
}

void Solution::Initialize( char const *fileName, double temperature, double airDens, double relHum )
{
    double amb_Value[N_SPC];
    std::ifstream file;

    file.open( fileName );

    if ( file.is_open() )
    {
        std::cout << "Reading ambient data from file: " << fileName << std::endl;
        std::string line;
        unsigned int i = 0;

        while ( std::getline( file, line ) ) {
            if ( line[0] != '#' ) {
                std::istringstream iss(line);
                iss >> amb_Value[i];
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
    SetShape( CO2  , NX, NY, amb_Value[  0] * airDens );
    SetShape( PPN  , NX, NY, amb_Value[  1] * airDens );
    SetShape( BrNO2, NX, NY, amb_Value[  2] * airDens );
    SetShape( IEPOX, NX, NY, amb_Value[  3] * airDens );
    SetShape( PMNN , NX, NY, amb_Value[  4] * airDens );
    SetShape( N2O  , NX, NY, amb_Value[  5] * airDens );
    SetShape( N    , NX, NY, amb_Value[  6] * airDens );
    SetShape( PAN  , NX, NY, amb_Value[  7] * airDens );
    SetShape( ALK4 , NX, NY, amb_Value[  8] * airDens );
    SetShape( MAP  , NX, NY, amb_Value[  9] * airDens );
    SetShape( MPN  , NX, NY, amb_Value[ 10] * airDens );
    SetShape( Cl2O2, NX, NY, amb_Value[ 11] * airDens );
    SetShape( ETP  , NX, NY, amb_Value[ 12] * airDens );
    SetShape( HNO2 , NX, NY, amb_Value[ 13] * airDens );
    SetShape( C3H8 , NX, NY, amb_Value[ 14] * airDens );
    SetShape( RA3P , NX, NY, amb_Value[ 15] * airDens );
    SetShape( RB3P , NX, NY, amb_Value[ 16] * airDens );
    SetShape( OClO , NX, NY, amb_Value[ 17] * airDens );
    SetShape( ClNO2, NX, NY, amb_Value[ 18] * airDens );
    SetShape( ISOP , NX, NY, amb_Value[ 19] * airDens );
    SetShape( HNO4 , NX, NY, amb_Value[ 20] * airDens );
    SetShape( MAOP , NX, NY, amb_Value[ 21] * airDens );
    SetShape( MP   , NX, NY, amb_Value[ 22] * airDens );
    SetShape( ClOO , NX, NY, amb_Value[ 23] * airDens );
    SetShape( RP   , NX, NY, amb_Value[ 24] * airDens );
    SetShape( BrCl , NX, NY, amb_Value[ 25] * airDens );
    SetShape( PP   , NX, NY, amb_Value[ 26] * airDens );
    SetShape( PRPN , NX, NY, amb_Value[ 27] * airDens );
    SetShape( SO4  , NX, NY, amb_Value[ 28] * airDens );
    SetShape( Br2  , NX, NY, amb_Value[ 29] * airDens );
    SetShape( ETHLN, NX, NY, amb_Value[ 30] * airDens );
    SetShape( MVKN , NX, NY, amb_Value[ 31] * airDens );
    SetShape( R4P  , NX, NY, amb_Value[ 32] * airDens );
    SetShape( C2H6 , NX, NY, amb_Value[ 33] * airDens );
    SetShape( RIP  , NX, NY, amb_Value[ 34] * airDens );
    SetShape( VRP  , NX, NY, amb_Value[ 35] * airDens );
    SetShape( ATOOH, NX, NY, amb_Value[ 36] * airDens );
    SetShape( IAP  , NX, NY, amb_Value[ 37] * airDens );
    SetShape( DHMOB, NX, NY, amb_Value[ 38] * airDens );
    SetShape( MOBA , NX, NY, amb_Value[ 39] * airDens );
    SetShape( MRP  , NX, NY, amb_Value[ 40] * airDens );
    SetShape( N2O5 , NX, NY, amb_Value[ 41] * airDens );
    SetShape( ISNOHOO, NX, NY, amb_Value[ 42] * airDens );
    SetShape( ISNP , NX, NY, amb_Value[ 43] * airDens );
    SetShape( ISOPNB, NX, NY, amb_Value[ 44] * airDens );
    SetShape( IEPOXOO, NX, NY, amb_Value[ 45] * airDens );
    SetShape( MACRNO2, NX, NY, amb_Value[ 46] * airDens );
    SetShape( ROH  , NX, NY, amb_Value[ 47] * airDens );
    SetShape( MOBAOO, NX, NY, amb_Value[ 48] * airDens );
    SetShape( DIBOO, NX, NY, amb_Value[ 49] * airDens );
    SetShape( PMN  , NX, NY, amb_Value[ 50] * airDens );
    SetShape( ISNOOB, NX, NY, amb_Value[ 51] * airDens );
    SetShape( INPN , NX, NY, amb_Value[ 52] * airDens );
    SetShape( H    , NX, NY, amb_Value[ 53] * airDens );
    SetShape( BrNO3, NX, NY, amb_Value[ 54] * airDens );
    SetShape( PRPE , NX, NY, amb_Value[ 55] * airDens );
    SetShape( MVKOO, NX, NY, amb_Value[ 56] * airDens );
    SetShape( Cl2  , NX, NY, amb_Value[ 57] * airDens );
    SetShape( ISOPND, NX, NY, amb_Value[ 58] * airDens );
    SetShape( HOBr , NX, NY, amb_Value[ 59] * airDens );
    SetShape( A3O2 , NX, NY, amb_Value[ 60] * airDens );
    SetShape( PROPNN, NX, NY, amb_Value[ 61] * airDens );
    SetShape( GLYX , NX, NY, amb_Value[ 62] * airDens );
    SetShape( MAOPO2, NX, NY, amb_Value[ 63] * airDens );
    SetShape( CH4  , NX, NY, amb_Value[ 64] * airDens );
    SetShape( GAOO , NX, NY, amb_Value[ 65] * airDens );
    SetShape( B3O2 , NX, NY, amb_Value[ 66] * airDens );
    SetShape( ACET , NX, NY, amb_Value[ 67] * airDens );
    SetShape( MACRN, NX, NY, amb_Value[ 68] * airDens );
    SetShape( CH2OO, NX, NY, amb_Value[ 69] * airDens );
    SetShape( MGLYOO, NX, NY, amb_Value[ 70] * airDens );
    SetShape( VRO2 , NX, NY, amb_Value[ 71] * airDens );
    SetShape( MGLOO, NX, NY, amb_Value[ 72] * airDens );
    SetShape( MACROO, NX, NY, amb_Value[ 73] * airDens );
    SetShape( PO2  , NX, NY, amb_Value[ 74] * airDens );
    SetShape( CH3CHOO, NX, NY, amb_Value[ 75] * airDens );
    SetShape( MAN2 , NX, NY, amb_Value[ 76] * airDens );
    SetShape( ISNOOA, NX, NY, amb_Value[ 77] * airDens );
    SetShape( H2O2 , NX, NY, amb_Value[ 78] * airDens );
    SetShape( PRN1 , NX, NY, amb_Value[ 79] * airDens );
    SetShape( ETO2 , NX, NY, amb_Value[ 80] * airDens );
    SetShape( KO2  , NX, NY, amb_Value[ 81] * airDens );
    SetShape( RCO3 , NX, NY, amb_Value[ 82] * airDens );
    SetShape( HC5OO, NX, NY, amb_Value[ 83] * airDens );
    SetShape( GLYC , NX, NY, amb_Value[ 84] * airDens );
    SetShape( ClNO3, NX, NY, amb_Value[ 85] * airDens );
    SetShape( RIO2 , NX, NY, amb_Value[ 86] * airDens );
    SetShape( R4N1 , NX, NY, amb_Value[ 87] * airDens );
    SetShape( HOCl , NX, NY, amb_Value[ 88] * airDens );
    SetShape( ATO2 , NX, NY, amb_Value[ 89] * airDens );
    SetShape( HNO3 , NX, NY, amb_Value[ 90] * airDens );
    SetShape( ISN1 , NX, NY, amb_Value[ 91] * airDens );
    SetShape( MAO3 , NX, NY, amb_Value[ 92] * airDens );
    SetShape( MRO2 , NX, NY, amb_Value[ 93] * airDens );
    SetShape( INO2 , NX, NY, amb_Value[ 94] * airDens );
    SetShape( HAC  , NX, NY, amb_Value[ 95] * airDens );
    SetShape( HC5  , NX, NY, amb_Value[ 96] * airDens );
    SetShape( MGLY , NX, NY, amb_Value[ 97] * airDens );
    SetShape( ISOPNBO2, NX, NY, amb_Value[ 98] * airDens );
    SetShape( ISOPNDO2, NX, NY, amb_Value[ 99] * airDens );
    SetShape( R4O2 , NX, NY, amb_Value[100] * airDens );
    SetShape( R4N2 , NX, NY, amb_Value[101] * airDens );
    SetShape( BrO  , NX, NY, amb_Value[102] * airDens );
    SetShape( RCHO , NX, NY, amb_Value[103] * airDens );
    SetShape( MEK  , NX, NY, amb_Value[104] * airDens );
    SetShape( ClO  , NX, NY, amb_Value[105] * airDens );
    SetShape( MACR , NX, NY, amb_Value[106] * airDens );
    SetShape( SO2  , NX, NY, amb_Value[107] * airDens );
    SetShape( MVK  , NX, NY, amb_Value[108] * airDens );
    SetShape( ALD2 , NX, NY, amb_Value[109] * airDens );
    SetShape( MCO3 , NX, NY, amb_Value[110] * airDens );
    SetShape( CH2O , NX, NY, amb_Value[111] * airDens );
    SetShape( H2O  , NX, NY, amb_Value[112] * airDens );
    SetShape( Br   , NX, NY, amb_Value[113] * airDens );
    SetShape( NO   , NX, NY, amb_Value[114] * airDens );
    SetShape( NO3  , NX, NY, amb_Value[115] * airDens );
    SetShape( Cl   , NX, NY, amb_Value[116] * airDens );
    SetShape( O    , NX, NY, amb_Value[117] * airDens );
    SetShape( O1D  , NX, NY, amb_Value[118] * airDens );
    SetShape( O3   , NX, NY, amb_Value[119] * airDens );
    SetShape( HO2  , NX, NY, amb_Value[120] * airDens );
    SetShape( NO2  , NX, NY, amb_Value[121] * airDens );
    SetShape( OH   , NX, NY, amb_Value[122] * airDens );
    SetShape( HBr  , NX, NY, amb_Value[123] * airDens );
    SetShape( HCl  , NX, NY, amb_Value[124] * airDens );
    SetShape( CO   , NX, NY, amb_Value[125] * airDens );
    SetShape( MO2  , NX, NY, amb_Value[126] * airDens );
    SetShape( ACTA , NX, NY, amb_Value[127] * airDens );
    SetShape( EOH  , NX, NY, amb_Value[128] * airDens );
    SetShape( H2   , NX, NY, amb_Value[129] * airDens );
    SetShape( HCOOH, NX, NY, amb_Value[130] * airDens );
    SetShape( MOH  , NX, NY, amb_Value[131] * airDens );
    SetShape( N2   , NX, NY, amb_Value[132] * airDens );
    SetShape( O2   , NX, NY, amb_Value[133] * airDens );
    SetShape( RCOOH, NX, NY, amb_Value[134] * airDens );

    /* Setting up ambient water vapor */
    for ( unsigned int i = 0; i < NX; i++ ) {
        for ( unsigned int j = 0; j < NY; j++ ) {
            H2O[j][i] = (relHum/((double) 100.0) * \
                             pSat_H2Ol( temperature ) / ( kB * temperature )) / 1.00E+06;
            /* RH_w = x_H2O*P/Psat_H2Ol = [H2O]*1E6 * kB*T/Psat_H2Ol */

        }
    }

    /* Liquid/solid species */
    SetShape( SO4L , NX, NY, (double) 0.0 );
    SetShape( H2OL , NX, NY, (double) 0.0 );
    SetShape( H2OS , NX, NY, (double) 0.0 );
    SetShape( HNO3L, NX, NY, (double) 0.0 );
    SetShape( HNO3S, NX, NY, (double) 0.0 );
    SetShape( HClL , NX, NY, (double) 0.0 );
    SetShape( HOClL, NX, NY, (double) 0.0 );
    SetShape( HBrL , NX, NY, (double) 0.0 );
    SetShape( HOBrL, NX, NY, (double) 0.0 );

    /* Aerosols */
    SetShape( Soot, NX, NY, (double) 0.0 );
    SetShape( SLA , NX, NY, (double) 0.0 );
    SetShape( SPA , NX, NY, (double) 0.0 );

}

void Solution::GetData( double varArray[], double fixArray[], unsigned int i, unsigned int j )
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

}

/* End of Structure.cpp */
