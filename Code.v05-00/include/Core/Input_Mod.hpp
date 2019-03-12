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
        bool        SIMULATION_MONTECARLO;
        int         SIMULATION_MCRUNS;
        std::string SIMULATION_OUTPUT_FOLDER;
        bool        SIMULATION_OVERWRITE;
        std::string SIMULATION_RUN_DIRECTORY;
        std::string SIMULATION_INPUT_BACKG_COND;
        bool        SIMULATION_SAVE_FORWARD;
        std::string SIMULATION_FORWARD_FILENAME;
        bool        SIMULATION_ADJOINT;
        std::string SIMULATION_ADJOINT_FILENAME;

        /* ========================================== */
        /* ---- PARAMETER MENU ---------------------- */
        /* ========================================== */

        bool        PARAMETER_PLUMEPROCESS_RANGE;
        std::string PARAMETER_PLUMEPROCESS_UNIT;
        Vector_1D   PARAMETER_PLUMEPROCESS;
        bool        PARAMETER_TEMPERATURE_RANGE;
        std::string PARAMETER_TEMPERATURE_UNIT;
        Vector_1D   PARAMETER_TEMPERATURE;
        bool        PARAMETER_RHW_RANGE;
        std::string PARAMETER_RHW_UNIT;
        Vector_1D   PARAMETER_RHW;
        bool        PARAMETER_SHEAR_RANGE;
        std::string PARAMETER_SHEAR_UNIT;
        Vector_1D   PARAMETER_SHEAR;
        bool        PARAMETER_LATITUDE_RANGE;
        std::string PARAMETER_LATITUDE_UNIT;
        Vector_1D   PARAMETER_LATITUDE;
        bool        PARAMETER_LONGITUDE_RANGE;
        std::string PARAMETER_LONGITUDE_UNIT;
        Vector_1D   PARAMETER_LONGITUDE;
        bool        PARAMETER_PRESSURE_RANGE;
        std::string PARAMETER_PRESSURE_UNIT;
        Vector_1D   PARAMETER_PRESSURE;
        bool        PARAMETER_EDAY_RANGE;
        std::string PARAMETER_EDAY_UNIT;
        Vector_1D   PARAMETER_EDAY;
        bool        PARAMETER_ETIME_RANGE;
        std::string PARAMETER_ETIME_UNIT;
        Vector_1D   PARAMETER_ETIME;
        bool        PARAMETER_BACKG_NOX_RANGE;
        std::string PARAMETER_BACKG_NOX_UNIT;
        Vector_1D   PARAMETER_BACKG_NOX;
        bool        PARAMETER_BACKG_HNO3_RANGE;
        std::string PARAMETER_BACKG_HNO3_UNIT;
        Vector_1D   PARAMETER_BACKG_HNO3;
        bool        PARAMETER_BACKG_O3_RANGE;
        std::string PARAMETER_BACKG_O3_UNIT;
        Vector_1D   PARAMETER_BACKG_O3;
        bool        PARAMETER_BACKG_CO_RANGE;
        std::string PARAMETER_BACKG_CO_UNIT;
        Vector_1D   PARAMETER_BACKG_CO;
        bool        PARAMETER_BACKG_CH4_RANGE;
        std::string PARAMETER_BACKG_CH4_UNIT;
        Vector_1D   PARAMETER_BACKG_CH4;
        bool        PARAMETER_BACKG_SO2_RANGE;
        std::string PARAMETER_BACKG_SO2_UNIT;
        Vector_1D   PARAMETER_BACKG_SO2;
        bool        PARAMETER_EI_NOX_RANGE;
        std::string PARAMETER_EI_NOX_UNIT;
        Vector_1D   PARAMETER_EI_NOX;
        bool        PARAMETER_EI_CO_RANGE;
        std::string PARAMETER_EI_CO_UNIT;
        Vector_1D   PARAMETER_EI_CO;
        bool        PARAMETER_EI_UHC_RANGE;
        std::string PARAMETER_EI_UHC_UNIT;
        Vector_1D   PARAMETER_EI_UHC;
        bool        PARAMETER_EI_SO2_RANGE;
        std::string PARAMETER_EI_SO2_UNIT;
        Vector_1D   PARAMETER_EI_SO2;
        bool        PARAMETER_EI_SO2TOSO4_RANGE;
        std::string PARAMETER_EI_SO2TOSO4_UNIT;
        Vector_1D   PARAMETER_EI_SO2TOSO4;
        bool        PARAMETER_EI_SOOT_RANGE;
        std::string PARAMETER_EI_SOOT_UNIT;
        Vector_1D   PARAMETER_EI_SOOT;
        bool        PARAMETER_EI_SOOTRAD_RANGE;
        std::string PARAMETER_EI_SOOTRAD_UNIT;
        Vector_1D   PARAMETER_EI_SOOTRAD;
        bool        PARAMETER_FF_RANGE;
        std::string PARAMETER_FF_UNIT;
        Vector_1D   PARAMETER_FF;
        
        /* ========================================== */
        /* ---- TRANSPORT MENU ---------------------- */
        /* ========================================== */

        bool        TRANSPORT_TRANSPORT;
        bool        TRANSPORT_FILL;
        RealDouble  TRANSPORT_TIMESTEP;
        bool        TRANSPORT_UPDRAFT;
        RealDouble  TRANSPORT_UPDRAFT_TIMESCALE;
        RealDouble  TRANSPORT_UPDRAFT_VELOCITY;

        /* ========================================== */
        /* ---- CHEMISTRY MENU ---------------------- */
        /* ========================================== */

        bool        CHEMISTRY_CHEMISTRY;
        bool        CHEMISTRY_HETCHEM;
        std::string CHEMISTRY_JRATE_FOLDER;
        RealDouble  CHEMISTRY_TIMESTEP;

        /* ========================================== */
        /* ---- AEROSOL MENU ------------------------ */
        /* ========================================== */

        bool        AEROSOL_GRAVSETTLING;
        bool        AEROSOL_COAGULATION;
        RealDouble  AEROSOL_COAGULATION_TIMESTEP;
        bool        AEROSOL_ICE_GROWTH;
        
        /* ========================================== */
        /* ---- METEOROLOGY MENU -------------------- */
        /* ========================================== */

        bool        MET_LOADMET;
        std::string MET_FILENAME;
        bool        MET_LOADTEMP;
        bool        MET_LOADH2O;
        bool        MET_FIXDEPTH;
        RealDouble  MET_DEPTH;
        bool        MET_FIXLAPSERATE;
        RealDouble  MET_LAPSERATE;
        bool        MET_DIURNAL;
        
        /* ========================================== */
        /* ---- DIAGNOSTIC MENU --------------------- */
        /* ========================================== */

        std::string DIAG_FILENAME;

        /* ========================================== */
        /* ---- TIMESERIES MENU --------------------- */
        /* ========================================== */

        bool             TS_SPEC;
        std::string      TS_FILENAME;
        std::vector<int> TS_SPECIES;
        RealDouble       TS_FREQ;
        bool             TS_AERO;
        std::string      TS_AERO_FILENAME;
        std::vector<int> TS_AEROSOL;
        RealDouble       TS_AERO_FREQ;

        /* ========================================== */
        /* ---- PROD & LOSS MENU -------------------- */
        /* ========================================== */

        bool PL_PL;
        bool PL_O3;


};


#endif /* INPUT_MOD_H_INCLUDED */


