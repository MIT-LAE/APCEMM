/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* PrintMessage Program File                                        */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : PrintMessage.cpp                          */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <string>

void PrintMessage( bool doPrint )
{
    std::string welcomeMessage, authorsMessage;

    welcomeMessage = "\n\n\
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\
            ~~~~~~~~                                  ~~~~~~~~\n\
            ~~~~~~~~             APCEMM               ~~~~~~~~\n\
            ~~~~~~~~                                  ~~~~~~~~\n\
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\0";

/*
                            _\
                           | \ \
                          _|  \____________________________\
                         /    o  o  o  o  o  o  o  o  o  |_\ \
                         \_________________________________/ \
                                         /    /\
                                        /    /\
                                       /    /\
                                      /____/\0
*/

    authorsMessage = "\
            \n   Version: 5.0\n\
            \n   Author: Thibaud M. Fritz (fritzt@mit.edu),\
            \n     with contributions from:\
            \n         - Sebastian D. Eastham (seastham@mit.edu),\
            \n         - Raymond L. Speth (speth@mit.edu).\n\
            \n   Corresponding author: Sebastian D. Eastham (seastham@mit.edu)\n\
            \n   This project was funded by NASA and developed at \
            \n   the Laboratory for Aviation and the Environment,\
            \n   Department of Aeronautics and Astronautics\
            \n   Massachusetts Institute of Technology\
            \n   Cambridge, MA, USA\n\0";

    std::cout << welcomeMessage << std::endl;

    if ( doPrint )
        std::cout << authorsMessage << std::endl;

}
