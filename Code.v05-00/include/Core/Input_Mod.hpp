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
#include "Util/ForwardDecl.hpp"

enum class epm_type { EPM_ORIGINAL, EPM_EXTERNAL, EPM_NEW_PHYSICS };

struct OptInput
{

    OptInput();
    /* ========================================== */
    /* ---- SIMULATION MENU --------------------- */
    /* ========================================== */

    int          SIMULATION_OMP_NUM_THREADS;
    std::string  SIMULATION_OUTPUT_FOLDER;
    bool         SIMULATION_OVERWRITE;
    std::string  SIMULATION_INPUT_BACKG_COND;
    std::string  SIMULATION_INPUT_ENG_EI;
    bool         SIMULATION_FORCE_SEED;
    unsigned int SIMULATION_SEED_VALUE;
    epm_type     SIMULATION_EPM_TYPE;
    std::string  SIMULATION_EXTERNAL_EPM_NETCDF_FILENAME;

    /* PARAMETER MENU is read into an Input object
     * because it describes the simulation scenario, whereas
     * OptInput describes the model configuration. */

    /* ========================================== */
    /* ---- TRANSPORT MENU ---------------------- */
    /* ========================================== */

    bool        TRANSPORT_TRANSPORT;
    double      TRANSPORT_TIMESTEP;

    /* ========================================== */
    /* ---- AEROSOL MENU ------------------------ */
    /* ========================================== */

    bool        AEROSOL_GRAVSETTLING;
    bool        AEROSOL_COAGULATION_SOLID;
    bool        AEROSOL_ICE_GROWTH;
    double      AEROSOL_ICE_GROWTH_TIMESTEP;

    /* ========================================== */
    /* ---- METEOROLOGY MENU -------------------- */
    /* ========================================== */

    std::string MET_FILENAME;
    double      MET_DT;
    bool        MET_TEMPTIMESERIES;
    bool        MET_INTERPTEMPDATA;
    bool        MET_RHTIMESERIES;
    bool        MET_INTERPRHDATA;
    bool        MET_SHEARTIMESERIES;
    bool        MET_INTERPSHEARDATA;
    bool        MET_VERTVELOCTIMESERIES;
    bool        MET_INTERPVERTVELOC;
    bool        MET_ENABLE_TEMP_PERTURB;
    double      MET_TEMP_PERTURB_AMPLITUDE;

    /* ========================================== */
    /* ---- TIMESERIES MENU --------------------- */
    /* ========================================== */

    bool             TS_AERO;
    std::string      TS_AERO_FILENAME;
    double           TS_AERO_FREQ;

    /* =============================================== */
    /* ---- ADVANCED OPTIONS MENU -------------------- */
    /* =============================================== */

    unsigned int ADV_GRID_NX;
    unsigned int ADV_GRID_NY;
    double ADV_GRID_XLIM_RIGHT;
    double ADV_GRID_XLIM_LEFT;
    double ADV_GRID_YLIM_UP;
    double ADV_GRID_YLIM_DOWN;
    double ADV_AMBIENT_LAPSERATE;
    double ADV_TROPOPAUSE_PRESSURE;
    double ADV_EP_N_REF;
    double ADV_EP_WINGSPAN_REF;
    bool ADV_EP_N_POSTJET_OVERRIDE;
    double ADV_EP_N_POSTJET;
    bool ADV_SAVE_PSD_GRID;

};

#endif /* INPUT_MOD_H_INCLUDED */


