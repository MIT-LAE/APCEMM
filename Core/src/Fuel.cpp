/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Fuel Program File                                                */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Fuel.cpp                                  */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Fuel.hpp"

Fuel::Fuel( const char *fuelChem )
{
    /* Constructor */
    GetAtoms( fuelChem );

    FSC = 1600; /* [ppm] */
}

Fuel::~Fuel( )
{
    /* Destructor */
}

void Fuel::GetAtoms( const char *fuelChem )
{
    atomC = 0;
    atomH = 0;
    atomN = 0;
    atomS = 0;

    std::string fuelChem_low;

    char *temp = strdup(fuelChem);
    unsigned char *tptr = (unsigned char *)temp;
    while (*tptr) {
        if (((*tptr >= 'A') && (*tptr <= 'Z')) || ((*tptr >= 'a') && (*tptr <= 'z')) || ((*tptr >= '0') && (*tptr <= '9')) || (*tptr == '.')) {
            *tptr = tolower(*tptr);
            fuelChem_low += *tptr;
        }
        tptr++;
    }

    const char *atoms = {"hcns"}; 
    std::size_t found, next_found;

    found = fuelChem_low.find_first_of( atoms );

    std::string substr;
    std::string::size_type sz;
    while (found != std::string::npos) {
        next_found = fuelChem_low.find_first_of( atoms, found+1 );
        substr = fuelChem_low.substr(found+1, next_found-1);
        if (fuelChem_low[found] == 'c') {
            atomC = std::stod(substr, &sz);
        }
        if (fuelChem_low[found] == 'h') {
            atomH = std::stod(substr, &sz);
        }
        if (fuelChem_low[found] == 'n') {
            atomN = std::stod(substr, &sz);
        }
        if (fuelChem_low[found] == 's') {
            atomS = std::stod(substr, &sz);
        }
        found = fuelChem_low.find_first_of( atoms, found+1 );
    }

    if ( atomC == 0 ) {
        std::cout << "Fuel doesn't contain any carbon: " << fuelChem << std::endl;
        return;
    }
    if ( atomH == 0 ) {
        std::cout << "Fuel doesn't contain any hydrogen: " << fuelChem << std::endl;
        return;
    }

}

