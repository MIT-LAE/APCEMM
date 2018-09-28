/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*                        AIrcraft Microphysics                     */
/*                              (AIM)                               */
/*                                                                  */
/* Coagulation Header File                                          */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 9/27/2018                                 */
/* File                 : Coagulation.hpp                           */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef COAGULATION_H_INCLUDED
#define COAGULATION_H_INCLUDED

#include <iostream>
#include <vector>
#include <cstring>

#include "../../Headers/PhysConstant.hpp"
#include "../../Headers/PhysFunction.hpp"

typedef double RealDouble;
typedef std::vector<RealDouble> Vector_1D;
typedef std::vector<Vector_1D> Vector_2D;

namespace AIM
{

    class Coagulation
    {

        public:

            Coagulation( );
            Coagulation( const char* phase, Vector_1D binCenters, RealDouble temperature_K = 220.0, RealDouble pressure_Pa = 220.0E+02 );
            ~Coagulation( );
            Coagulation( const Coagulation& k );
            Coagulation& operator=( const Coagulation& k );
            Vector_2D getKernel() const;

        protected:

            Vector_2D Kernel;

        private:

    };

}

#endif /* COAGULATION_H_INCLUDED */

