/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Engine and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Engine Header File                                               */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Engine.cpp                                */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Engine.hpp"

Engine::Engine( const char *fileName, const char *engineName, double tempe_K, double pres_Pa, double relHum_w, double vFlight_ms )
{
    /* Constructor */

    std::ifstream enginefile;
    OpenFile( fileName, enginefile );
    
    std::string search( engineName );
    std::string delimiter = ",";

    std::string endch(1,search.back());
    if ( delimiter.compare(endch) != 0)
        search += delimiter;

    size_t pos;
    bool found = 0;
    std::string line;
    std::string idle, approach, climbout, takeoff;
    while ((enginefile.eof() == 0) && (found == 0)) {
        getline( enginefile, line );
        pos = line.find( search );
        if ( pos != std::string::npos ) {
            approach = line;
            getline( enginefile, line );
            climbout = line;
            getline( enginefile, line );
            takeoff = line;
            getline( enginefile, line );
            idle = line;
            found = 1;
        }
    }

    CloseFile( enginefile );

    ratedThrust.push_back(4);
    LTO_fuelflow.push_back(4);
    LTO_NOx.push_back(4);
    LTO_CO.push_back(4);
    LTO_HC.push_back(4);


    /* BOEING FUEL FLOW METHOD 2 (BFFM2) */
    /* See:
     * SAGE ( System for assessing Aviation's Global Emissions ),
     * FAA, 
     * Version 1.5, Technical Manual (2005)
     */

    /* Conversion of uninstalled conditions to installed conditions */
    /* Adjustment/correction factor for installation effects (engine air bleed) */
    std::vector<double> fuelAdjustmentFactor = {1.100, 1.020, 1.013, 1.010};

    std::vector<std::string> tokens;
    std::string token;
    std::string::size_type sz;
    std::istringstream tokenStream;

    /* Idle */
    tokenStream.str(idle);
    while ( std::getline( tokenStream, token, delimiter.c_str()[0] ) )
        tokens.push_back( token );
  
    ratedThrust[0]  = 0.07;
    LTO_fuelflow[0] = std::stod(tokens[7]) * fuelAdjustmentFactor[0];
    LTO_NOx[0]      = std::stod(tokens[4]);
    LTO_CO[0]       = std::stod(tokens[2]);
    LTO_HC[0]       = std::stod(tokens[3]);
    tokens.clear();
    tokenStream.clear();

    /* Approach */
    tokenStream.str(approach);
    while ( std::getline( tokenStream, token, delimiter.c_str()[0] ) )
        tokens.push_back( token );
  
    ratedThrust[1]  = 0.30;
    LTO_fuelflow[1] = std::stod(tokens[7]) * fuelAdjustmentFactor[1];
    LTO_NOx[1]      = std::stod(tokens[4]);
    LTO_CO[1]       = std::stod(tokens[2]);
    LTO_HC[1]       = std::stod(tokens[3]);
    tokens.clear();
    tokenStream.clear();

    /* Climb out */
    tokenStream.str(climbout);
    while ( std::getline( tokenStream, token, delimiter.c_str()[0] ) )
        tokens.push_back( token );
  
    ratedThrust[2]  = 0.70;
    LTO_fuelflow[2] = std::stod(tokens[7]) * fuelAdjustmentFactor[2];
    LTO_NOx[2]      = std::stod(tokens[4]);
    LTO_CO[2]       = std::stod(tokens[2]);
    LTO_HC[2]       = std::stod(tokens[3]);
    tokens.clear();
    tokenStream.clear();

    /* Take off */
    tokenStream.str(takeoff);
    while ( std::getline( tokenStream, token, delimiter.c_str()[0] ) )
        tokens.push_back( token );
  
    ratedThrust[3]  = 1.00;
    LTO_fuelflow[3] = std::stod(tokens[7]) * fuelAdjustmentFactor[3];
    LTO_NOx[3]      = std::stod(tokens[4]);
    LTO_CO[3]       = std::stod(tokens[2]);
    LTO_HC[3]       = std::stod(tokens[3]);
    tokens.clear();
    tokenStream.clear();

    /* Check that all values are strictly positives */
    for ( unsigned int i = 0; i < 4; i++ ) {
        if ( LTO_fuelflow[i] <= 0.0 ) {
            std::cout << "LTO_fuelflow is negative for engine: " << engineName << " LTO index: " << i << ", fuelflow: " << LTO_fuelflow[i] << std::endl;
            LTO_fuelflow[i] = 0.1;
        }
        if ( LTO_NOx[i] <= 0.0 ) {
            std::cout << "LTO_NOx is negative for engine: " << engineName << " LTO index: " << i << ", EI_NOx: " << LTO_NOx[i] << std::endl;
            LTO_NOx[i] = 0.01;
        }
        if ( LTO_CO[i] <= 0.0 ) {
            std::cout << "LTO_CO is negative for engine: " << engineName << " LTO index: " << i << ", EI_CO: " << LTO_CO[i] << std::endl;
            LTO_CO[i] = 0.1;
        }
        if ( LTO_HC[i] <= 0.0 ) {
            std::cout << "LTO_HC is negative for engine: " << engineName << " LTO index: " << i << ", EI_HC: " << LTO_HC[i] << std::endl;
            LTO_CO[i] = 0.1;
        }
    }

    /* Finding log-log relationship between fuelflow and EI_NOx */
    std::vector<double> log_LTO_fuelflow(4);
    std::vector<double> log_LTO_NOx(4);

    std::vector<double> p(2); /* Parameters y = p(1)*x + p(2) */

    /* (1) Compute logs */
    for ( unsigned int i = 0; i < 4; i++ ) {
        log_LTO_fuelflow[i] = log10( LTO_fuelflow[i] );
        log_LTO_NOx[i]      = log10( LTO_NOx[i] );
    }

    /* (2) Build Vandermonde matrix */
    std::vector<std::vector<double> > V{{log_LTO_fuelflow[0], 1.0},
                                        {log_LTO_fuelflow[1], 1.0},
                                        {log_LTO_fuelflow[2], 1.0},
                                        {log_LTO_fuelflow[3], 1.0}};

    /* (3) Some linear algebra magic */
    std::vector<std::vector<double> > VV, invVV;
    std::vector<double> vect(2);

    for ( unsigned int i = 0; i < 2; i++ ) {
        VV.push_back( std::vector<double>( 2 ) );
        invVV.push_back( std::vector<double>( 2 ) );
        for ( unsigned int j = 0; j < 2; j++ ) {
            for ( unsigned int k = 0; k < 4; k++ )
                VV[i][j] += V[k][i] * V[k][j];
        }
    }

    double determinant = VV[0][0] * VV[1][1] - VV[0][1] * VV[1][0];
    if ( abs(determinant) < 1E-20 )
        std::cout << "Matrix is badly-scaled or singular" << std::endl;

    invVV[0][0] =  VV[1][1] / determinant;
    invVV[0][1] = -VV[0][1] / determinant;
    invVV[1][0] = -VV[1][0] / determinant;
    invVV[1][1] =  VV[0][0] / determinant;

    for ( unsigned int i = 0; i < 2; i++ ) {
        for ( unsigned int k = 0; k < 4; k++ )
            vect[i] = vect[i] + V[k][i] * log_LTO_NOx[k];
    }

    for ( unsigned int i = 0; i < 2; i++ ) {
        for ( unsigned int k = 0; k < 2; k++ )
            p[i] += invVV[i][k] * vect[k];
    }


    double fuelflow = 0.8;

    double delta = pres_Pa / PRES_SL;
    double theta = tempe_K / TEMP_SL;
    double mach  = vFlight_ms / sqrt( gamma * R_Air * tempe_K );

    double fuelflow_factor;
    
    /* Fuel flow rate converted to SLS-ISA conditions (installed engine) */
    fuelflow_factor = fuelflow / delta * pow( theta, 3.8 ) * exp( 0.2 * mach );
    EI_NOx = pow( 10.0, log10( fuelflow_factor ) * p[0] + p[1] );

    /* Cruise correction for NOx */
    double beta, Pv, H;
    beta = 7.90298 * ( 1.0 - ( 373.16 ) / ( tempe_K + 0.01 ) ) + 3.00571 + \
           5.02808 * log10(( 373.16 ) / ( tempe_K + 0.01 )) + \
           1.3816E-07 * ( 1.0 - ( pow( 10.0, 11.344 * ( 1.0 - (( tempe_K + 0.01 ) / ( 373.16 ))) ) )) + \
           8.1328E-03 * (( pow( 10.0, 3.49149 * ( 1.0 -  ( 373.16 ) / ( tempe_K + 0.01 )))) - 1.0 );
    Pv = 0.014504 * pow( 10.0, beta );
    H = -19.0 * ( 0.37318 * relHum_w/100.0 * Pv ) / ( 14.696 * delta - relHum_w/100.0 * Pv );

    EI_NOx *= exp( H ) * pow( pow( delta, 1.02 ) / pow( theta, 3.3 ), 0.5 );

    std::cout << "NOx Emission Index: " << EI_NOx << std::endl;

    EI_NO   = NOxtoNO   * EI_NOx;
    EI_NO2  = NOxtoNO2  * EI_NOx;
    EI_HNO2 = NOxtoHNO2 * EI_NOx;


    double line1, line2, line3;
    line1 = (log10( LTO_CO[1] ) - log10( LTO_CO[0] )) / ( log10( LTO_fuelflow[1] ) - log10( LTO_fuelflow[0]) );
    line2 = log10( LTO_fuelflow[0] );
    line3 = log10( LTO_CO[0] );

    double horzline;
    horzline = log10( LTO_CO[2] ) + log10( LTO_CO[3] ) / 2.0;

    double intercept;
    intercept = ( 2 * log10( LTO_fuelflow[0] ) * line1 + log10( LTO_CO[2] ) + log10( LTO_CO[3] ) - 2 * log10( LTO_CO[0]) ) / ( 2 * line1 );

    std::cout << intercept << std::endl;
    if ( intercept > log10( LTO_fuelflow[2] ) )
        intercept = log10( LTO_fuelflow[2] );
    else if ( intercept < log10( LTO_fuelflow[1] ) && ( line1 < 0 ) ) {
        horzline = log10( LTO_CO[1] );
        intercept = log10( LTO_fuelflow[1] );
    }
    else if ( line1 >= 0 ) {
        line1 = 0;
        line2 = 0;
        line3 = horzline;
        intercept = log10( LTO_fuelflow[1] );
    }




    /* Cruise correction for CO */
    EI_CO *= pow( theta, 3.3 ) / pow(delta, 1.02 );

    /* Cruise correction for HC */
    EI_HC *= pow( theta, 3.3 ) / pow(delta, 1.02 );


}

Engine::~Engine( )
{
    /* Destructor */
}

bool Engine::CheckFile( const char *fileName )
{
    std::ifstream f( fileName );
    return f.is_open();
}

void Engine::OpenFile( const char *fileName, std::ifstream &file )
{
    const bool file_exists = CheckFile( fileName );

    if (!file_exists) {
        std::string const currFunc("Engine::OpenFile");
        std::cout << "ERROR: In " << currFunc << ": Can't read (" << fileName << ")" << std::endl;
        return;
    }
    std::cout << "Reading emission data from file: " << fileName << std::endl;

    file.open( fileName );

}

void Engine::CloseFile( std::ifstream &file )
{
    file.close( );
}

/* End of Engine.cpp */

