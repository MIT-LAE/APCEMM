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

#include <vector>
#include <cstring>
#include "Util/ForwardDecl.hpp"

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
        Coagulation( const char* phase, Vector_1D const &bin_Centers_1, Vector_1D &bin_VCenters_1, double rho_1, Vector_1D const &bin_Centers_2, double rho_2, double temperature_K_, double pressure_Pa_ );
        Coagulation( const char* phase, Vector_1D const &bin_Centers_1, Vector_1D &bin_VCenters_1, double rho_1, double temperature_K_, double pressure_Pa_ );
        Coagulation( const char* phase, Vector_1D const &bin_Centers_1, double rho_1, double bin_Centers_2, double rho_2, double temperature_K_, double pressure_Pa_ );
            
        ~Coagulation( );
        Coagulation( const Coagulation& k );
        Coagulation& operator=( const Coagulation& k );
        void buildBeta( const Vector_1D &bin_Centers );
        void buildF( const Vector_1D &bin_VCenters );
        void buildF( const Vector_3D &bin_VCenters, const UInt jNy, const UInt iNx );
        Vector_2D getKernel() const;
        Vector_1D getKernel_1D() const;
        Vector_2D getBeta() const;
        Vector_3D getF() const;
        
        Vector_3D f;
        std::vector<std::vector<std::vector<UInt> > > indices;

    protected:

        Vector_2D Kernel;
        Vector_1D Kernel_1D;
        Vector_2D beta;

    private:

        const double A0 = 5.07;
        const double A1 = -5.94;
        const double A2 = 7.27;
        const double A3 = -5.29;


};

#endif /* COAGULATION_H_INCLUDED */

