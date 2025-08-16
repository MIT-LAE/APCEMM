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

#include <string>
#include <vector>
#include <unordered_map>
#include "Util/ForwardDecl.hpp"

enum class epm_type { EPM_ORIGINAL, EPM_EXTERNAL, EPM_NEW_PHYSICS };

struct OptInput
{

    OptInput();
    /* ========================================== */
    /* ---- SIMULATION MENU --------------------- */
    /* ========================================== */

    int          SIMULATION_OMP_NUM_THREADS;
    bool         SIMULATION_PARAMETER_SWEEP;
    bool         SIMULATION_MONTECARLO;
    unsigned int SIMULATION_MCRUNS;
    std::string  SIMULATION_OUTPUT_FOLDER;
    bool         SIMULATION_OVERWRITE;
    bool         SIMULATION_THREADED_FFT;
    bool         SIMULATION_USE_FFTW_WISDOM;
    std::string  SIMULATION_DIRECTORY_W_WRITE_PERMISSION;
    std::string  SIMULATION_INPUT_BACKG_COND;
    std::string  SIMULATION_INPUT_ENG_EI;
    bool         SIMULATION_SAVE_FORWARD;
    std::string  SIMULATION_FORWARD_FILENAME;
    bool         SIMULATION_ADJOINT;
    std::string  SIMULATION_ADJOINT_FILENAME;
    bool         SIMULATION_BOXMODEL;
    std::string  SIMULATION_BOX_FILENAME;
    bool         SIMULATION_FORCE_SEED;
    int          SIMULATION_SEED_VALUE;
    epm_type     SIMULATION_EPM_TYPE;
    std::string  SIMULATION_EXTERNAL_EPM_NETCDF_FILENAME;

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
    bool        MET_ENABLE_TEMP_PERTURB;
    double      MET_TEMP_PERTURB_AMPLITUDE;
    double      MET_TEMP_PERTURB_TIMESCALE;

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

    unsigned int ADV_GRID_NX;
    unsigned int ADV_GRID_NY;
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
    double ADV_EP_N_REF;
    double ADV_EP_WINGSPAN_REF;
    bool ADV_EP_N_POSTJET_OVERRIDE;
    double ADV_EP_N_POSTJET;

};

#endif /* INPUT_MOD_H_INCLUDED */


