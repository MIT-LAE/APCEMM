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
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Ring.hpp"

Ring::Ring( )
{
    /* Constructor */

} /* End of Ring::Ring */

void Ring::Create( double a, double b )
{

    if ( a < 0 || b < 0 ) {
        std::cout << "Horizontal and/or vertical axis are negative!" << std::endl;
        return;
    }

    horizontalAxis = a;
    verticalAxis = b;
    
} /* End of Ring::Create */

Ring::~Ring( )
{
    /* Destructor */

} /* End of Ring::~Ring */

void Ring::Print( ) const
{
    std::cout << "Ring's horizontal and vertical axis: " << horizontalAxis << ", " << verticalAxis << " [m]" << std::endl;

} /* End of Ring::Print */

double Ring::GetHAxis( ) const
{
    
    return horizontalAxis;

} /* End of Ring::GetHAxis */

double Ring::GetVAxis( ) const
{
    
    return verticalAxis;

} /* End of Ring::GetVAxis */

Ring Ring::Copy( )
{

    Ring copy;

    copy.Create( horizontalAxis, verticalAxis );

    return copy;

} /* End of Ring::Copy */

/* End of Ring.cpp */
