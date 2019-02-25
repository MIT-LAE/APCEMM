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
#include <algorithm>
#include <vector>
#include <cstring>

#include "Util/ForwardDecl.hpp"
#include "Util/PhysConstant.hpp"
#include "Util/PhysFunction.hpp"
#include "AIM/buildKernel.hpp"

namespace AIM
{

    /* We need forward declaration here */
    class Aerosol;
    class Grid_Aerosol;

    class Coagulation;

}

class AIM::Coagulation
{

    friend class Aerosol;
    friend class Grid_Aerosol;

    public:

        Coagulation( );
        Coagulation( const char* phase, Vector_1D const &bin_Centers_1, Vector_1D &bin_VCenters_1, RealDouble rho_1, Vector_1D const &bin_Centers_2, RealDouble rho_2, RealDouble temperature_K_, RealDouble pressure_Pa_ );
        Coagulation( const char* phase, Vector_1D const &bin_Centers_1, Vector_1D &bin_VCenters_1, RealDouble rho_1, RealDouble temperature_K_, RealDouble pressure_Pa_ );
        Coagulation( const char* phase, Vector_1D const &bin_Centers_1, RealDouble rho_1, RealDouble bin_Centers_2, RealDouble rho_2, RealDouble temperature_K_, RealDouble pressure_Pa_ );
            
        ~Coagulation( );
        Coagulation( const Coagulation& k );
        Coagulation& operator=( const Coagulation& k );
        void buildBeta( const Vector_1D &bin_Centers );
        void buildF( Vector_1D &bin_VCenters );
        void buildF( Vector_3D &bin_VCenters, const UInt jNy, const UInt iNx );
        Vector_2D getKernel() const;
        Vector_1D getKernel_1D() const;
        Vector_2D getBeta() const;
        Vector_3D getF() const;
        void printKernel_1D( const char* fileName ) const;
        void printKernel_2D( const char* fileName ) const;
        
        Vector_3D f;
        std::vector<std::vector<std::vector<UInt> > > indices;

    protected:

        Vector_2D Kernel;
        Vector_1D Kernel_1D;
        Vector_2D beta;

    private:

        const RealDouble A0 = 5.07;
        const RealDouble A1 = -5.94;
        const RealDouble A2 = 7.27;
        const RealDouble A3 = -5.29;


};

#endif /* COAGULATION_H_INCLUDED */

