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

    RealDouble backgNOx_;
    RealDouble backgHNO3_;
    RealDouble backgO3_;
    RealDouble backgCO_;
    RealDouble backgCH4_;

    std::string fileName_;
    std::string fileName_ADJ_;

    public:

        Input( const RealDouble temperature_K, \
               const RealDouble pressure_Pa,   \
               const RealDouble relHumidity_w, \
               const RealDouble long_deg,      \
               const RealDouble lat_deg,       \
               const unsigned int dayGMT,      \
               const RealDouble emissionTime,  \
               const Vector_1D emissionInput,  \
               const Vector_1D backgMixRatio,  \
               const std::string fileName,     \
               const std::string fileName_ADJ );
        Input( unsigned int iCase,          \
               const Vector_2D &parameters, \
               const std::string fileName,  \
               const std::string fileName_ADJ );

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

        RealDouble backgNOx() const { return backgNOx_; }
        RealDouble backgHNO3() const { return backgHNO3_; }
        RealDouble backgO3() const { return backgO3_; }
        RealDouble backgCO() const { return backgCO_; }
        RealDouble backgCH4() const { return backgCH4_; }

        std::string fileName() const { return fileName_; }
        std::string fileName_ADJ() const { return fileName_ADJ_; }
        const char* fileName2char() const { return fileName_.c_str(); }
        const char* fileName_ADJ2char() const { return fileName_ADJ_.c_str(); }

};

#endif /* INPUT_H_INCLUDED */
