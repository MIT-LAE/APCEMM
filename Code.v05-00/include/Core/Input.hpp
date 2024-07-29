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
#include <unordered_map>
#include "Util/ForwardDecl.hpp"

class Input
{

    UInt Case_;

    double simulationTime_;

    double temperature_K_;
    double relHumidity_w_;
    double horizDiff_;
    double vertiDiff_;
    double shear_;

    double longitude_deg_;
    double latitude_deg_;
    double pressure_Pa_;

    UInt emissionDOY_;
    double emissionTime_;

    double EI_NOx_;
    double EI_CO_;
    double EI_HC_;
    double EI_SO2_;
    double EI_SO2TOSO4_;
    double EI_Soot_;
    double sootRad_;

    double fuelFlow_;
    double aircraftMass_;

    double backgNOx_;
    double backgHNO3_;
    double backgO3_;
    double backgCO_;
    double backgCH4_;
    double backgSO2_;
    
    double flightSpeed_;
    double numEngines_;
    double wingspan_;
    double coreExitTemp_;
    double bypassArea_;

    std::string fileName_;
    std::string fileName_ADJ_;
    std::string fileName_BOX_;
    std::string fileName_micro_;
    std::string author_;

    double nBV_;

    UInt emissionDay_;
    UInt emissionMonth_;

    private:
        void checkInputValidity();
        void adjustLatLong();
        void calculate_emissionMonth();

    public:

        Input( unsigned int iCase,               \
               const Vector_2D &parameters,      \
               const std::string fileName,       \
               const std::string fileName_ADJ,   \
               const std::string fileName_BOX,   \
               const std::string fileName_micro, \
               const std::string author          );
        Input( unsigned int iCase,               \
                const std::vector<std::unordered_map<std::string, double>> &parameters,      \
                const std::string fileName,       \
                const std::string fileName_ADJ,   \
                const std::string fileName_BOX,   \
                const std::string fileName_micro, \
                const std::string author          );

        ~Input();
        UInt Case() const { return Case_; }

        double simulationTime() const { return simulationTime_; }
        double temperature_K() const { return temperature_K_; }
        double pressure_Pa() const { return pressure_Pa_; }
        double relHumidity_w() const { return relHumidity_w_; }
        double horizDiff() const { return horizDiff_; }
        double vertiDiff() const { return vertiDiff_; }
        double shear() const { return shear_; }
        inline double nBV() const { return nBV_; }
        
        double longitude_deg() const { return longitude_deg_; }
        double latitude_deg() const { return latitude_deg_; }

        UInt emissionDOY() const { return emissionDOY_; }
        UInt emissionDay() const { return emissionDay_; }
        UInt emissionMonth() const { return emissionMonth_; }
        double emissionTime() const { return emissionTime_; }

        double EI_NOx() const { return EI_NOx_; }
        double EI_CO() const { return EI_CO_; }
        double EI_HC() const { return EI_HC_; }
        double EI_SO2() const { return EI_SO2_; }
        double EI_SO2TOSO4() const { return EI_SO2TOSO4_; }
        double EI_Soot() const { return EI_Soot_; }
        double sootRad() const { return sootRad_; }
        
        double fuelFlow() const { return fuelFlow_; }

        double aircraftMass() const { return aircraftMass_; }
        double flightSpeed() const { return flightSpeed_; }
        double numEngines() const { return numEngines_; }
        double wingspan() const { return wingspan_; }
        double coreExitTemp() const { return coreExitTemp_; }
        double bypassArea() const { return bypassArea_; }
        
        double backgNOx() const { return backgNOx_; }
        double backgHNO3() const { return backgHNO3_; }
        double backgO3() const { return backgO3_; }
        double backgCO() const { return backgCO_; }
        double backgCH4() const { return backgCH4_; }
        double backgSO2() const { return backgSO2_; }

        std::string fileName() const { return fileName_; }
        std::string fileName_ADJ() const { return fileName_ADJ_; }
        std::string fileName_BOX() const { return fileName_BOX_; }
        std::string fileName_micro() const { return fileName_micro_; }
        std::string author() const { return author_; }
        const char* fileName2char() const { return fileName_.c_str(); }
        const char* fileName_ADJ2char() const { return fileName_ADJ_.c_str(); }
        const char* fileName_BOX2char() const { return fileName_BOX_.c_str(); }
        const char* fileName_micro2char() const { return fileName_micro_.c_str(); }

};

#endif /* INPUT_H_INCLUDED */
