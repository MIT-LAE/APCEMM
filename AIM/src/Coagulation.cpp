/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*                        AIrcraft Microphysics                     */
/*                              (AIM)                               */
/*                                                                  */
/* Coagulation Program File                                         */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 9/27/2018                                 */
/* File                 : Coagulation.cpp                           */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Coagulation.hpp"

AIM::Coagulation::Coagulation( )
{

    /* Default Constructor */

} /* End of Coagulation::Coagulation */

AIM::Coagulation::Coagulation( const char* phase, Vector_1D const &bin_Centers_1, RealDouble rho_1, Vector_1D const &bin_Centers_2, RealDouble rho_2, RealDouble temperature_K_, RealDouble pressure_Pa_ ): Kernel( bin_Centers_1.size() )
{

    /* Constructor */

    RealDouble temperature_K = temperature_K_;
    RealDouble pressure_Pa = pressure_Pa_;

    /* Set shape */
    for ( unsigned int iBin_1 = 0; iBin_1 < bin_Centers_1.size(); iBin_1++ ) {
        Kernel.push_back( Vector_1D( bin_Centers_2.size() ) );
    }

    if ( ( strcmp ( phase, "liq" ) == 0) || ( strcmp( phase, "liquid" ) == 0) || ( strcmp( phase, "sulfate" ) == 0) ) {
        /* Liquid-liquid coagulation */

        /* Declare and assign coagulation kernels for different physical processes */
        Vector_2D K_Brow = buildBrownianKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_2, rho_2 );
        Vector_2D K_DE   = buildDEKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_2, rho_2, K_Brow );

        /* Fill in total coagulation kernel */
        for ( unsigned int iBin_1 = 0; iBin_1 < bin_Centers_1.size(); iBin_1++ ) {
            for ( unsigned int iBin_2 = 0; iBin_2 < bin_Centers_2.size(); iBin_2++ ) {
                Kernel[iBin_1][iBin_2] = ( K_Brow[iBin_1][iBin_2] + K_DE[iBin_1][iBin_2] ) * 1.00E+06; /* Conversion from m^3/s to cm^3/s */
            }
        }

    } 
    else if ( ( strcmp ( phase, "ice" ) == 0) ) { 
        /* Ice-Ice coagulation */

        /* Declare and assign coagulation kernels for different physical processes */
        Vector_2D K_Brow = buildBrownianKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_2, rho_2 );
        Vector_2D K_DE   = buildDEKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_2, rho_2, K_Brow );
        Vector_2D K_GC   = buildGCKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_2, rho_2 );
        Vector_2D K_TI   = buildTIKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_2, rho_2 );
        Vector_2D K_TS   = buildTSKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_2, rho_2 );
        
        /* Fill in total coagulation kernel */
        for ( unsigned int iBin_1 = 0; iBin_1 < bin_Centers_1.size(); iBin_1++ ) {
            for ( unsigned int iBin_2 = 0; iBin_2 < bin_Centers_2.size(); iBin_2++ ) {
                Kernel[iBin_1][iBin_2] = ( K_Brow[iBin_1][iBin_2] + K_DE[iBin_1][iBin_2] + K_GC[iBin_1][iBin_2] + K_TI[iBin_1][iBin_2] + K_TS[iBin_1][iBin_2] ) * 1.00E+06; /* Conversion from m^3/s to cm^3/s */
            }
        }

    }
    else if ( ( strcmp ( phase, "soot" ) == 0) || ( strcmp( phase, "bc" ) == 0) ) {
        /* Soot-soot coagulation */

        /* Declare and assign coagulation kernels for different physical processes */
        Vector_2D K_Brow = buildBrownianKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_2, rho_2 );
        Vector_2D K_DE   = buildDEKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_2, rho_2, K_Brow );

        /* Fill in total coagulation kernel */
        for ( unsigned int iBin_1 = 0; iBin_1 < bin_Centers_1.size(); iBin_1++ ) {
            for ( unsigned int iBin_2 = 0; iBin_2 < bin_Centers_2.size(); iBin_2++ ) {
                Kernel[iBin_1][iBin_2] = ( K_Brow[iBin_1][iBin_2] + K_DE[iBin_1][iBin_2] ) * 1.00E+06; /* Conversion from m^3/s to cm^3/s */
            }
        }

    }
    else {
        std::cout << "\nIn AIM::Coagulation::Coagulation: phase " << phase << " is not defined.";
        std::cout << "\nOptions are: liquid, ice or soot.";
    }


} /* End of Coagulation::Coagulation */
        
AIM::Coagulation::Coagulation( const char* phase, Vector_1D const &bin_Centers_1, RealDouble rho_1, RealDouble bin_Centers_2, RealDouble rho_2, RealDouble temperature_K_, RealDouble pressure_Pa_ ):
    Kernel_1D( bin_Centers_1.size() )
{

    /* Constructor */

    RealDouble temperature_K = temperature_K_;
    RealDouble pressure_Pa = pressure_Pa_;

    /* Declare coagulation kernel and set shape */

    if ( ( strcmp ( phase, "liq" ) == 0) || ( strcmp( phase, "liquid" ) == 0) ) {
        /* Liquid-liquid coagulation */

        /* Declare and assign coagulation kernels for different physical processes */
        Vector_1D K_Brow = buildBrownianKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_2, rho_2 );
        Vector_1D K_DE   = buildDEKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_2, rho_2, K_Brow );

        /* Fill in total coagulation kernel */
        for ( unsigned int iBin_1 = 0; iBin_1 < bin_Centers_1.size(); iBin_1++ ) {
            Kernel_1D[iBin_1] = ( K_Brow[iBin_1] + K_DE[iBin_1] ) * 1.00E+06; /* Conversion from m^3/s to cm^3/s */
        }

    } 
    else if ( ( strcmp ( phase, "ice" ) == 0) ) { 
        /* Ice-Ice coagulation */

        /* Declare and assign coagulation kernels for different physical processes */
        Vector_1D K_Brow = buildBrownianKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_2, rho_2 );
        Vector_1D K_DE   = buildDEKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_2, rho_2, K_Brow );
        Vector_1D K_GC   = buildGCKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_2, rho_2 );
        Vector_1D K_TI   = buildTIKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_2, rho_2 );
        Vector_1D K_TS   = buildTSKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_2, rho_2 );
        
        /* Fill in total coagulation kernel */
        for ( unsigned int iBin_1 = 0; iBin_1 < bin_Centers_1.size(); iBin_1++ ) {
            Kernel_1D[iBin_1] = ( K_Brow[iBin_1] + K_DE[iBin_1] + K_GC[iBin_1] + K_TI[iBin_1] + K_TS[iBin_1] ) * 1.00E+06; /* Conversion from m^3/s to cm^3/s */
        }

    }
    else if ( ( strcmp ( phase, "soot" ) == 0) || ( strcmp( phase, "bc" ) == 0) ) {
        /* Soot-soot coagulation */

        /* Declare and assign coagulation kernels for different physical processes */
        Vector_1D K_Brow = buildBrownianKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_2, rho_2 );
        Vector_1D K_DE   = buildDEKernel( temperature_K, pressure_Pa, bin_Centers_1, rho_1, bin_Centers_2, rho_2, K_Brow );

        /* Fill in total coagulation kernel */
        for ( unsigned int iBin_1 = 0; iBin_1 < bin_Centers_1.size(); iBin_1++ ) {
            Kernel_1D[iBin_1] = ( K_Brow[iBin_1] + K_DE[iBin_1] ) * 1.00E+06; /* Conversion from m^3/s to cm^3/s */
        }

    }
    else {
        std::cout << "\nIn AIM::Coagulation::Coagulation: phase " << phase << " is not defined.";
        std::cout << "\nOptions are: liquid, ice or soot.";
    }

} /* End of Coagulation::Coagulation */

AIM::Coagulation::~Coagulation( )
{

    /* Destructor */

} /* End of Coagulation::~Coagulation */

AIM::Coagulation::Coagulation( const Coagulation &k )
{

    Kernel = k.Kernel;
    Kernel_1D = k.Kernel_1D;

} /* End of Coagulation::Coagulation */

AIM::Coagulation& AIM::Coagulation::operator=( const Coagulation &k )
{

    if ( &k == this )
        return *this;

    Kernel = k.Kernel;
    Kernel_1D = k.Kernel_1D;
    return *this;

} /* End of Coagulation::operator= */

Vector_2D AIM::Coagulation::getKernel() const
{

    return Kernel;

} /* End of Coagulation::getKernel */
        
Vector_1D AIM::Coagulation::getKernel_1D() const
{

    return Kernel_1D;

} /* End of Coagulation::getKernel_1D */

void AIM::Coagulation::printKernel_1D( const char* fileName ) const
{

    printKernel2File( Kernel_1D, fileName ); 

    return;

} /* End of Coagulation::printKernel_1D */

void AIM::Coagulation::printKernel_2D( const char* fileName) const
{

    printKernel2File( Kernel, fileName ); 

    return;

} /* End of Coagulation::printKernel_2D */

/* End of Coagulation.cpp */
