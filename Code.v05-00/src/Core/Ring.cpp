/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Ring Program File                                                */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 8/12/2018                                 */
/* File                 : Ring.cpp                                  */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include "Core/Ring.hpp"

Ring::Ring( )
{
    /* Default Constructor */

    horizontalAxis = 0.0;
    verticalAxis = 0.0;

} /* End of Ring::Ring */

Ring::Ring( double a, double b )
{

    /* Constructor */

    if ( a < 0 || b < 0 ) {
        std::cout << "Horizontal and/or vertical axis are negative!" << std::endl;
        a = 0.0;
        b = 0.0;
    }

    horizontalAxis = a;
    verticalAxis = b;
    
} /* End of Ring::Ring */

Ring::Ring( const Ring &r )
{

    /* Constructor */

    horizontalAxis = r.getHAxis();
    verticalAxis = r.getVAxis();
    
} /* End of Ring::Ring */

Ring& Ring::operator=( const Ring &r )
{

    if ( &r == this )
        return *this;

    horizontalAxis = r.getHAxis();
    verticalAxis = r.getVAxis();

    return *this;

} /* End of Ring::operator= */

Ring::~Ring( )
{
    /* Destructor */

} /* End of Ring::~Ring */

void Ring::Assign( double a, double b )
{

    if ( a < 0 || b < 0 ) {
        std::cout << "Horizontal and/or vertical axis are negative!" << std::endl;
        a = 0.0;
        b = 0.0;
    }

    horizontalAxis = a;
    verticalAxis = b;
    

} /* End of Ring::Assign */

void Ring::Print( ) const
{
    std::cout << "Ring's horizontal and vertical axis: " << horizontalAxis << ", " << verticalAxis << " [m]" << std::endl;

} /* End of Ring::Print */

double Ring::getHAxis( ) const
{
    
    return horizontalAxis;

} /* End of Ring::getHAxis */

double Ring::getVAxis( ) const
{
    
    return verticalAxis;

} /* End of Ring::getVAxis */

/* End of Ring.cpp */
