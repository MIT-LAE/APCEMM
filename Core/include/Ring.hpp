/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Ring Header File                                                 */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 8/12/2018                                 */
/* File                 : Ring.hpp                                  */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef RING_H_INCLUDED
#define RING_H_INCLUDED

#include <iostream>

class Ring
{

    public:

        Ring( );
        void Create( double a, double b );
        ~Ring( );
        void Print( ) const;
        double GetHAxis( ) const;
        double GetVAxis( ) const;
        Ring Copy( );

        double horizontalAxis;
        double verticalAxis;

    private:

};

#endif /* RING_H_INCLUDED */
