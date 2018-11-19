/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Input Header File                                                */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 11/19/2018                                */
/* File                 : Input.hpp                                 */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef INPUT_H_INCLUDED
#define INPUT_H_INCLUDED

#include <iostream>
#include <vector>
#include "Util/ForwardDecl.hpp"

class Input
{
        
     RealDouble temperature_K_;
     RealDouble pressure_Pa_;
     RealDouble relHumidity_w_;

     RealDouble longitude_deg_;
     RealDouble latitude_deg_;

     UInt dayGMT_;
     RealDouble emissionTime_;

     RealDouble EI_NOx_;
     RealDouble EI_CO_;
     RealDouble EI_HC_;
     RealDouble EI_Soot_;
     RealDouble sootRad_;

     RealDouble fuelFlow_;


    public:

        Input( const RealDouble temperature_K, \
               const RealDouble pressure_Pa,   \
               const RealDouble relHumidity_w, \
               const RealDouble long_deg,      \
               const RealDouble lat_deg,       \
               const unsigned int dayGMT,      \
               const RealDouble emissionTime,  \
               const Vector_1D emissionInput);
        Input( unsigned int iCase, \
               const Vector_2D &parameters );

        ~Input();

        RealDouble temperature_K() const { return temperature_K_; }
        RealDouble pressure_Pa() const { return pressure_Pa_; }
        RealDouble relHumidity_w() const { return relHumidity_w_; }
        
        RealDouble longitude_deg() const { return longitude_deg_; }
        RealDouble latitude_deg() const { return latitude_deg_; }

        UInt dayGMT() const { return dayGMT_; }
        RealDouble emissionTime() const { return emissionTime_; }

        RealDouble EI_NOx() const { return EI_NOx_; }
        RealDouble EI_CO() const { return EI_CO_; }
        RealDouble EI_HC() const { return EI_HC_; }
        RealDouble EI_Soot() const { return EI_Soot_; }
        RealDouble sootRad() const { return sootRad_; }
        
        RealDouble fuelFlow() const { return fuelFlow_; }


};

#endif /* INPUT_H_INCLUDED */
