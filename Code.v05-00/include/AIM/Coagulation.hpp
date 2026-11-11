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
#include "Util/ForwardDecl.hpp"

namespace AIM
{

    /* We need forward declaration here */
    class Aerosol;
    class Grid_Aerosol;

    class Coagulation;

    /* Colliding particle pair. It selects the physical processes in the kernel. */
    enum class CoagulationPair {
        LiquidLiquid, /* Brownian + convective diffusion enhancement */
        LiquidSoot,   /* Brownian + convective diffusion enhancement */
        IceIce        /* Same as above, plus gravitational collection, turbulent inertial motion and turbulent shear */
    };

}

class AIM::Coagulation
{

    friend class Aerosol;
    friend class Grid_Aerosol;

    public:

        Coagulation( CoagulationPair pair, Vector_1D const &bin_Centers_1, Vector_1D const &bin_VCenters_1, double rho_1, double temperature_K_, double pressure_Pa_ );
        Coagulation( CoagulationPair pair, Vector_1D const &bin_Centers_1, double rho_1, double radius_2, double rho_2, double temperature_K_, double pressure_Pa_ );

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

};

#endif /* COAGULATION_H_INCLUDED */

