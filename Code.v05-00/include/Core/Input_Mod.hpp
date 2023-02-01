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
#include <unordered_map>
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
        bool        SIMULATION_THREADED_FFT;
        bool        SIMULATION_USE_FFTW_WISDOM;
        std::string SIMULATION_DIRECTORY_W_WRITE_PERMISSION;
        std::string SIMULATION_INPUT_BACKG_COND;
        std::string SIMULATION_INPUT_ENG_EI;
        bool        SIMULATION_SAVE_FORWARD;
        std::string SIMULATION_FORWARD_FILENAME;
        bool        SIMULATION_ADJOINT;
        std::string SIMULATION_ADJOINT_FILENAME;
        bool        SIMULATION_BOXMODEL;
        std::string SIMULATION_BOX_FILENAME;

        /* ========================================== */
        /* ---- PARAMETER MENU ---------------------- */
        /* ========================================== */

        std::unordered_map<std::string, Vector_1D> PARAMETER_PARAM_MAP;    
            
        /* ========================================== */
        /* ---- TRANSPORT MENU ---------------------- */
        /* ========================================== */

        bool        TRANSPORT_TRANSPORT;
        bool        TRANSPORT_FILL;
        RealDouble  TRANSPORT_TIMESTEP;
        bool        TRANSPORT_PART_FLUX;
        bool        TRANSPORT_UPDRAFT;
        RealDouble  TRANSPORT_UPDRAFT_TIMESCALE;
        RealDouble  TRANSPORT_UPDRAFT_VELOCITY;

        /* ========================================== */
        /* ---- CHEMISTRY MENU ---------------------- */
        /* ========================================== */

        bool        CHEMISTRY_CHEMISTRY;
        bool        CHEMISTRY_HETCHEM;
        RealDouble  CHEMISTRY_TIMESTEP;
        std::string CHEMISTRY_JRATE_FOLDER;

        /* ========================================== */
        /* ---- AEROSOL MENU ------------------------ */
        /* ========================================== */

        bool        AEROSOL_GRAVSETTLING;
        bool        AEROSOL_COAGULATION_SOLID;
        bool        AEROSOL_COAGULATION_LIQUID;
        RealDouble  AEROSOL_COAGULATION_TIMESTEP;
        bool        AEROSOL_ICE_GROWTH;
        
        /* ========================================== */
        /* ---- METEOROLOGY MENU -------------------- */
        /* ========================================== */

        bool        MET_LOADMET;
        std::string MET_FILENAME;
        RealDouble  MET_DT;
        bool        MET_LOADTEMP;
        bool        MET_TEMPTIMESERIES;
        bool        MET_LOADRH;
        bool        MET_RHTIMESERIES;
        bool        MET_LOADSHEAR;
        bool        MET_SHEARTIMESERIES;
        bool        MET_FIXDEPTH;
        RealDouble  MET_DEPTH;
        bool        MET_FIXLAPSERATE;
        RealDouble  MET_LAPSERATE;
        bool        MET_DIURNAL;
        bool        MET_ENABLE_TEMP_PERTURB;
        double        MET_TEMP_PERTURB_AMPLITUDE;
        double        MET_TEMP_PERTURB_TIMESCALE;
        
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


