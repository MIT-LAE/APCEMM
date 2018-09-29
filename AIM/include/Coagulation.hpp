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

#include "buildKernel.hpp"
#include "../../Headers/PhysConstant.hpp"
#include "../../Headers/PhysFunction.hpp"

typedef double RealDouble;
typedef std::vector<RealDouble> Vector_1D;
typedef std::vector<Vector_1D> Vector_2D;

namespace AIM
{

    class Coagulation;

}

class AIM::Coagulation
{

    public:

        Coagulation( );
        Coagulation( const char* phase, Vector_1D const &bin_Centers_1, RealDouble rho_1, Vector_1D const &bin_Centers_2, RealDouble rho_2, RealDouble temperature_K_, RealDouble pressure_Pa_ );
        Coagulation( const char* phase, Vector_1D const &bin_Centers_1, RealDouble rho_1, RealDouble bin_Centers_2, RealDouble rho_2, RealDouble temperature_K_, RealDouble pressure_Pa_ );
            
        ~Coagulation( );
        Coagulation( const Coagulation& k );
        Coagulation& operator=( const Coagulation& k );
        Vector_2D getKernel() const;
        Vector_1D getKernel_1D() const;
        void printKernel_1D( const char* fileName ) const;
        void printKernel_2D( const char* fileName ) const;

    protected:

        Vector_2D Kernel;
        Vector_1D Kernel_1D;

    private:

};

#endif /* COAGULATION_H_INCLUDED */

