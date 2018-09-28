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

namespace AIM
{

    Coagulation::Coagulation( )
    {

        /* Default Constructor */

    } /* End of Coagulation::Coagulation */

    Coagulation::Coagulation( const char* phase, Vector_1D binCenters, RealDouble temperature_K_, RealDouble pressure_Pa_ )
    {

        /* Constructor */

        RealDouble temperature_K = temperature_K_;
        RealDouble pressure_Pa = pressure_Pa_;

        if ( (strcmp ( phase, "liq") == 0) || (strcmp( phase, "liquid") == 0) ) {
        } 
        else if ( (strcmp ( phase, "ice") == 0) ) { 
        }
        else if ( (strcmp ( phase, "soot") == 0) || (strcmp( phase, "bc") == 0) ) {
        }
        else {
            std::cout << "\nIn AIM::Coagulation::Coagulation: phase " << phase << " is not defined.";
            std::cout << "\nOptions are: liquid, ice or soot.";
        }


    } /* End of Coagulation::Coagulation */

    Coagulation::~Coagulation( )
    {

        /* Destructor */

    } /* End of Coagulation::~Coagulation */

    Coagulation::Coagulation( const Coagulation &k )
    {

        Kernel = k.Kernel;

    } /* End of Coagulation::Coagulation */

    Coagulation& Coagulation::operator=( const Coagulation &k )
    {

        if ( &k == this )
            return *this;

        Kernel = k.Kernel;
        return *this;

    } /* End of Coagulation::operator= */

    Vector_2D Coagulation::getKernel() const
    {

        return Kernel;

    } /* End of Coagulation::getKernel */

}

/* End of Coagulation.cpp */
