/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Meteorology Header File                                          */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 11/2/2018                                 */
/* File                 : Meteorology.hpp                           */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef METEOROLOGY_H_INCLUDED
#define METEOROLOGY_H_INCLUDED

#include <iostream>
#include "APCEMM.h"
#ifdef OMP
    #include "omp.h"
#endif /* OMP */

/* Include Parameters.hpp for multithreading option */
#include "Core/Parameters.hpp"
#include "Util/ForwardDecl.hpp"
#include "Core/Mesh.hpp"
#include "Core/Input_Mod.hpp"
#include "Util/PhysConstant.hpp"
#include "Util/MetFunction.hpp"
/*#include <netcdfcpp.h>*/
#include <netcdf>
#include <limits>

using namespace netCDF;
using namespace netCDF::exceptions;

class Meteorology
{

    public:

        Meteorology( );
        Meteorology( const OptInput &USERINPUT,      \
                     const RealDouble solarTime_h,   \
                     const Mesh &m,                  \
                     const RealDouble temperature_K,  \
                     const RealDouble relHumidity_i, \
                     const RealDouble pressure_Pa,   \
                     const RealDouble shear_pm,   \
                     const bool DBG = 0 );
        Meteorology( const Meteorology &met );
        ~Meteorology( );

        void Update( const OptInput &USERINPUT, const RealDouble solarTime_h, \
                     const RealDouble simTime_h, \
                     const Mesh &m, const RealDouble dTrav_x, const RealDouble dTrav_y, \
		     const bool DBG = 0 );
        
        void loadMet(const OptInput &USERINPUT,      \
                     const RealDouble solarTime_h,   \
                     const Mesh &m,                  \
                     const RealDouble temperature_K,  \
                     const RealDouble relHumidity_i, \
                     const RealDouble pressure_Pa,   \
                     const RealDouble shear_pm,   \
                     const bool DBG = 0 );
        
        void applyTempPerturb(const OptInput& USERINPUT, const Mesh& m);
        RealDouble alt( UInt j ) const { return alt_[j]; }
	RealDouble press( UInt j ) const { return press_[j]; }
        RealDouble shear( UInt j ) const { return shear_[j]; }

        RealDouble temp( UInt j, UInt i ) const { return temp_[j][i]; }
        RealDouble airDens( UInt j, UInt i ) const { return airDens_[j][i]; }
        RealDouble H2O( UInt j, UInt i ) const { return H2O_[j][i]; }

        const Vector_2D& Temp() const { return temp_; }
        const Vector_1D& Press() const { return press_; }
        const Vector_1D& Shear() const { return shear_; }

        friend class Solution;

        /* Temperature, pressure and humidity fields if input from user-defined file */
        RealDouble alt_user;
        RealDouble temp_user;
        RealDouble pres_user;
        RealDouble RHw_user;
        RealDouble satdepth_user;
        RealDouble S_user;

    protected:

        /* Met input type */
        UInt TYPE;

        /* Ambient parameters */
        const RealDouble TEMPERATURE;
        const RealDouble PRESSURE;
        const RealDouble RHI;
        const RealDouble SHEAR;
        RealDouble ALTITUDE;

        /* Temperature lapse rate */
        RealDouble LAPSERATE;

        /* Diurnal temperature variations */
        RealDouble DIURNAL_AMPL;
        RealDouble DIURNAL_PHASE;
        RealDouble diurnalPert;
        
        /* Non-uniformity length scales and temperature perturbations */
        RealDouble DELTAT;
        RealDouble TOP, BOT, LEFT, RIGHT;

        RealDouble RH;
        RealDouble RH_star;
        RealDouble RH_far;

        /* Assume that pressure only depends on the vertical coordinate */
        Vector_1D alt_;
	Vector_1D press_;
        Vector_1D shear_;

        /* Temperature, air density and humidity fields can potentially be
         * 2D fields */
        Vector_2D temp_;
        Vector_2D airDens_;
        Vector_2D H2O_;


        /* Temperature input from user can be 2D [T(z,t)] */
	Vector_2D temperature_store_;
        Vector_1D altitude_store_;
        Vector_2D shear_store_;

};


#endif /* METEOROLOGY_H_INCLUDED */
