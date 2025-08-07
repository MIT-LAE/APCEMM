#include "Core/Input_Mod.hpp"
#include "Core/Input.hpp"

#ifndef SIMVARSWRAPPER_H_INCLUDED
#define SIMVARSWRAPPER_H_INCLUDED

struct MPMSimVarsWrapper{ 
    /* ======================================================================= */
    /* --- Input options from the SIMULATION MENU ---------------------------- */
    /* ======================================================================= */

    bool RUN_BOXMODEL;
    bool BUILD_LUT;
    bool SAVE_FORWARD;
    bool ADJOINT;
    std::string BACKG_FILENAME;
    bool THREADED_FFT;
    bool USE_WISDOM;
    std::string FFTW_DIR;

    /* ======================================================================= */
    /* ---- Input options from the TRANSPORT MENU ---------------------------- */
    /* ======================================================================= */

    bool TRANSPORT;
    bool FILLNEG;

    bool UPDRAFT;
    double UPDRAFT_TIME;
    double UPDRAFT_VEL;

    /* ======================================================================= */
    /* ---- Input options from the CHEMISTRY MENU ---------------------------- */
    /* ======================================================================= */

    bool CHEMISTRY;
    bool HETCHEM;
    std::string JRATE_FOLDER;

    /* ======================================================================= */
    /* ---- Input options from the AEROSOL MENU ------------------------------ */
    /* ======================================================================= */

    bool GRAVSETTLING;
    bool ICE_COAG;
    bool LIQ_COAG;
    bool ICE_GROWTH;

    /* ======================================================================= */
    /* ---- Input options from the METEOROLOGY MENU -------------------------- */
    /* ======================================================================= */

    /* Input options from the METEOROLOGY MENU are read in the Meteorology
     * subroutine */

    bool TEMP_PERTURB;
    double metDepth;

    /* ======================================================================= */
    /* ---- Input options from the DIAGNOSTIC MENU --------------------------- */
    /* ======================================================================= */
    
    std::string DIAG_FILENAME;
    
    /* ======================================================================= */
    /* ---- Input options from the TIMESERIES MENU --------------------------- */
    /* ======================================================================= */

    std::string TS_FOLDER;
    bool TS_SPEC;
    std::string TS_SPEC_FILEPATH;
    std::vector<int> TS_SPEC_LIST;
    double TS_FREQ;

    bool TS_AERO;
    std::string TS_AERO_FILEPATH;
    std::vector<int> TS_AERO_LIST;
    double TS_AERO_FREQ;

    /* ======================================================================= */
    /* ---- Input options from the PROD & LOSS MENU -------------------------- */
    /* ======================================================================= */

    /* TODO: Implement PL rates */
    bool SAVE_PL;
    bool SAVE_O3PL;

    /* If DIAG_OUTPUT is turned on, make sure that output timestep is a
     * multiple of the dynamic timestep */

    double temperature_K;
    double pressure_Pa;
    double relHumidity_w;

    /* Compute relative humidity w.r.t ice */
    double relHumidity_i;

    MPMSimVarsWrapper() = default;
    MPMSimVarsWrapper(const Input& input, const OptInput& Input_Opt, const double depth_estimate);
    inline int LA_MICROPHYSICS() const { return LA_MICROPHYSICS_; }
    inline int PA_MICROPHYSICS() const { return PA_MICROPHYSICS_; }
    inline bool TRANSPORT_LA() const { return TRANSPORT_LA_; }
    inline bool TRANSPORT_PA() const { return TRANSPORT_PA_; }

    // Making these private since these were const in the original PlumeModel
    private:
        int LA_MICROPHYSICS_;
        int PA_MICROPHYSICS_;
        bool TRANSPORT_LA_;
        bool TRANSPORT_PA_;


};
#endif
