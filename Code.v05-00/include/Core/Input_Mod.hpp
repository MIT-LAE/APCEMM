/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Input_Mod Header File                                            */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 12/12/2018                                */
/* File                 : Input_Mod.hpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef INPUT_MOD_H_INCLUDED
#define INPUT_MOD_H_INCLUDED

#include <iostream>
#include <vector>
#include "Util/ForwardDecl.hpp"

class OptInput
{

    public:

        OptInput();
        ~OptInput();

        /* ========================================== */
        /* ---- SIMULATION MENU --------------------- */
        /* ========================================== */

        bool        SIMULATION_PARAMETER_SWEEP;
        std::string SIMULATION_OUTPUT_FOLDER;
        std::string SIMULATION_RUN_DIRECTORY;
        std::string SIMULATION_INPUT_BACKG_COND;
        bool        SIMULATION_SAVE_FORWARD;
        std::string SIMULATION_FORWARD_FILENAME;
        bool        SIMULATION_ADJOINT;
        std::string SIMULATION_ADJOINT_FILENAME;

        /* ========================================== */
        /* ---- PARAMETER MENU ---------------------- */
        /* ========================================== */

        bool        PARAMETER_TEMPERATURE_RANGE;
        const char* PARAMETER_TEMPERATURE_UNIT;
        Vector_1D   PARAMETER_TEMPERATURE;
        bool        PARAMETER_RHW_RANGE;
        const char* PARAMETER_RHW_UNIT;
        Vector_1D   PARAMETER_RHW;
        bool        PARAMETER_LATITUDE_RANGE;
        const char* PARAMETER_LATITUDE_UNIT;
        Vector_1D   PARAMETER_LATITUDE;
        bool        PARAMETER_LONGITUDE_RANGE;
        const char* PARAMETER_LONGITUDE_UNIT;
        Vector_1D   PARAMETER_LONGITUDE;
        bool        PARAMETER_PRESSURE_RANGE;
        const char* PARAMETER_PRESSURE_UNIT;
        Vector_1D   PARAMETER_PRESSURE;
        bool        PARAMETER_EDAY_RANGE;
        const char* PARAMETER_EDAY_UNIT;
        Vector_1D   PARAMETER_EDAY;
        bool        PARAMETER_ETIME_RANGE;
        const char* PARAMETER_ETIME_UNIT;
        Vector_1D   PARAMETER_ETIME;
        bool        PARAMETER_BACKG_NOX_RANGE;
        const char* PARAMETER_BACKG_NOX_UNIT;
        Vector_1D   PARAMETER_BACKG_NOX;
        bool        PARAMETER_BACKG_HNO3_RANGE;
        const char* PARAMETER_BACKG_HNO3_UNIT;
        Vector_1D   PARAMETER_BACKG_HNO3;
        bool        PARAMETER_BACKG_O3_RANGE;
        const char* PARAMETER_BACKG_O3_UNIT;
        Vector_1D   PARAMETER_BACKG_O3;
        bool        PARAMETER_BACKG_CO_RANGE;
        const char* PARAMETER_BACKG_CO_UNIT;
        Vector_1D   PARAMETER_BACKG_CO;
        bool        PARAMETER_BACKG_CH4_RANGE;
        const char* PARAMETER_BACKG_CH4_UNIT;
        Vector_1D   PARAMETER_BACKG_CH4;
        bool        PARAMETER_BACKG_SO2_RANGE;
        const char* PARAMETER_BACKG_SO2_UNIT;
        Vector_1D   PARAMETER_BACKG_SO2;
        bool        PARAMETER_EI_NOX_RANGE;
        const char* PARAMETER_EI_NOX_UNIT;
        Vector_1D   PARAMETER_EI_NOX;
        bool        PARAMETER_EI_CO_RANGE;
        const char* PARAMETER_EI_CO_UNIT;
        Vector_1D   PARAMETER_EI_CO;
        bool        PARAMETER_EI_UHC_RANGE;
        const char* PARAMETER_EI_UHC_UNIT;
        Vector_1D   PARAMETER_EI_UHC;
        bool        PARAMETER_EI_SO2_RANGE;
        const char* PARAMETER_EI_SO2_UNIT;
        Vector_1D   PARAMETER_EI_SO2;
        bool        PARAMETER_EI_SO2TOSO4_RANGE;
        const char* PARAMETER_EI_SO2TOSO4_UNIT;
        Vector_1D   PARAMETER_EI_SO2TOSO4;
        bool        PARAMETER_EI_SOOT_RANGE;
        const char* PARAMETER_EI_SOOT_UNIT;
        Vector_1D   PARAMETER_EI_SOOT;
        bool        PARAMETER_EI_SOOTRAD_RANGE;
        const char* PARAMETER_EI_SOOTRAD_UNIT;
        Vector_1D   PARAMETER_EI_SOOTRAD;
        


};


#endif /* INPUT_MOD_H_INCLUDED */


