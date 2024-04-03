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

struct OptInput
{

    OptInput();
    /* ========================================== */
    /* ---- SIMULATION MENU --------------------- */
    /* ========================================== */

    int         SIMULATION_OMP_NUM_THREADS;
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
    double      TRANSPORT_TIMESTEP;
    bool        TRANSPORT_UPDRAFT;
    double      TRANSPORT_UPDRAFT_TIMESCALE;
    double      TRANSPORT_UPDRAFT_VELOCITY;

    /* ========================================== */
    /* ---- CHEMISTRY MENU ---------------------- */
    /* ========================================== */

    bool        CHEMISTRY_CHEMISTRY;
    bool        CHEMISTRY_HETCHEM;
    double      CHEMISTRY_TIMESTEP;
    std::string CHEMISTRY_JRATE_FOLDER;

    /* ========================================== */
    /* ---- AEROSOL MENU ------------------------ */
    /* ========================================== */

    bool        AEROSOL_GRAVSETTLING;
    bool        AEROSOL_COAGULATION_SOLID;
    bool        AEROSOL_COAGULATION_LIQUID;
    double      AEROSOL_COAGULATION_TIMESTEP;
    bool        AEROSOL_ICE_GROWTH;
    double      AEROSOL_ICE_GROWTH_TIMESTEP;
    
    /* ========================================== */
    /* ---- METEOROLOGY MENU -------------------- */
    /* ========================================== */

    bool        MET_LOADMET;
    std::string MET_FILENAME;
    double      MET_DT;
    bool        MET_LOADTEMP;
    bool        MET_TEMPTIMESERIES;
    bool        MET_INTERPTEMPDATA;
    bool        MET_LOADRH;
    bool        MET_RHTIMESERIES;
    bool        MET_INTERPRHDATA;
    bool        MET_LOADSHEAR;
    bool        MET_SHEARTIMESERIES;
    bool        MET_INTERPSHEARDATA;
    bool        MET_LOADVERTVELOC;
    bool        MET_VERTVELOCTIMESERIES;
    bool        MET_INTERPVERTVELOC;
    bool        MET_FIXDEPTH;
    double      MET_DEPTH;
    double      MET_SUBSAT_RHI;
    bool        MET_FIXLAPSERATE;
    double      MET_LAPSERATE;
    bool        MET_DIURNAL;
    bool        MET_ENABLE_TEMP_PERTURB;
    double      MET_TEMP_PERTURB_AMPLITUDE;
    double      MET_TEMP_PERTURB_TIMESCALE;
    std::string MET_HUMIDSCAL_MODIFICATION_SCHEME;
    double      MET_HUMIDSCAL_CONST_RHI;
    double      MET_HUMIDSCAL_SCALING_A;
    double      MET_HUMIDSCAL_SCALING_B;
    
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
    double           TS_FREQ;
    bool             TS_AERO;
    std::string      TS_AERO_FILENAME;
    std::vector<int> TS_AEROSOL;
    double           TS_AERO_FREQ;

    /* ========================================== */
    /* ---- PROD & LOSS MENU -------------------- */
    /* ========================================== */

    bool PL_PL;
    bool PL_O3;

    /* =============================================== */
    /* ---- ADVANCED OPTIONS MENU -------------------- */
    /* =============================================== */

    int ADV_GRID_NX;
    int ADV_GRID_NY;
    double ADV_GRID_XLIM_RIGHT;
    double ADV_GRID_XLIM_LEFT;
    double ADV_GRID_YLIM_UP;
    double ADV_GRID_YLIM_DOWN;
    double ADV_CSIZE_DEPTH_BASE;
    double ADV_CSIZE_DEPTH_SCALING_FACTOR;
    double ADV_CSIZE_WIDTH_BASE;
    double ADV_CSIZE_WIDTH_SCALING_FACTOR;
    double ADV_AMBIENT_LAPSERATE;
    double ADV_TROPOPAUSE_PRESSURE;
        

};


#endif /* INPUT_MOD_H_INCLUDED */


