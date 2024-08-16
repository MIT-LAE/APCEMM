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
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include "Core/Fuel.hpp"

Fuel::Fuel( )
{

    /* Default Constructor */

} /* End of Fuel::Fuel */

Fuel::Fuel( const char *fuelChem )
{

    /* Constructor */
    getAtoms( fuelChem );

    FSC = 1600; /* [ppm] */

    std::string ChemFormula( fuelChem );

} /* End of Fuel::Fuel */

Fuel::Fuel( const Fuel &f )
{

    /* Constructor */
    
    atomC = f.getAtomC();
    atomH = f.getAtomH();
    atomN = f.getAtomN();
    atomS = f.getAtomS();
    FSC = f.getFSC();
    ChemFormula = f.getChemFormula();

} /* End of Fuel::Fuel */

Fuel& Fuel::operator=( const Fuel &f )
{

    if ( &f == this )
        return *this;

    atomC = f.getAtomC();
    atomH = f.getAtomH();
    atomN = f.getAtomN();
    atomS = f.getAtomS();
    FSC = f.getFSC();
    ChemFormula = f.getChemFormula();

    return *this;

} /* End of Fuel::operator= */

Fuel::~Fuel( )
{

    /* Destructor */

} /* End of Fuel::~Fuel */

void Fuel::getAtoms( const char *fuelChem )
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

    free((char*) temp); temp = NULL;
//    free((unsigned char*) tptr); tptr = NULL;

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

} /* End of Fuel::getAtoms */

double Fuel::getAtomC() const
{

    return atomC;

} /* End of Fuel::getAtomC */

double Fuel::getAtomH() const
{

    return atomH;

} /* End of Fuel::getAtomH */

double Fuel::getAtomN() const
{

    return atomN;

} /* End of Fuel::getAtomN */

double Fuel::getAtomS() const
{

    return atomS;

} /* End of Fuel::getAtomS */

double Fuel::getFSC() const
{

    return FSC;

} /* End of Fuel::getFSC */

void Fuel::setFSC( const double FSC_ )
{

    if ( FSC_ > 0.0E+00 )
        FSC = FSC_;

} /* End of Fuel::setFSC */

std::string Fuel::getChemFormula() const
{

    return ChemFormula;

} /* End of Fuel::getChemFormula */


/* End of Fuel.cpp */

