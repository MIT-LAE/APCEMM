#include "Core/Input_Mod.hpp"
#include "Core/Input.hpp"
#include "Core/Structure.hpp"

#ifndef SIMVARSWRAPPER_H_INCLUDED
#define SIMVARSWRAPPER_H_INCLUDED

struct MPMSimVarsWrapper{ 
    /* ======================================================================= */
    /* --- Input options from the SIMULATION MENU ---------------------------- */
    /* ======================================================================= */

    const bool RUN_BOXMODEL;
    const bool BUILD_LUT;
    const bool SAVE_FORWARD;
    const bool ADJOINT;
    const std::string BACKG_FILENAME;
    const bool THREADED_FFT;
    const bool USE_WISDOM;
    const std::string FFTW_DIR;

    /* ======================================================================= */
    /* ---- Input options from the TRANSPORT MENU ---------------------------- */
    /* ======================================================================= */

    const bool TRANSPORT;
    const bool FILLNEG;

    const bool UPDRAFT;
    const double UPDRAFT_TIME;
    const double UPDRAFT_VEL;

    /* ======================================================================= */
    /* ---- Input options from the CHEMISTRY MENU ---------------------------- */
    /* ======================================================================= */

    const bool CHEMISTRY;
    const bool HETCHEM;
    const std::string JRATE_FOLDER;

    /* ======================================================================= */
    /* ---- Input options from the AEROSOL MENU ------------------------------ */
    /* ======================================================================= */

    const bool GRAVSETTLING;
    const bool ICE_COAG;
    const bool LIQ_COAG;
    const bool ICE_GROWTH;

    /* ======================================================================= */
    /* ---- Input options from the METEOROLOGY MENU -------------------------- */
    /* ======================================================================= */

    /* Input options from the METEOROLOGY MENU are read in the Meteorology
     * subroutine */

    const bool TEMP_PERTURB;
    double metDepth;

    /* ======================================================================= */
    /* ---- Input options from the DIAGNOSTIC MENU --------------------------- */
    /* ======================================================================= */
    
    const std::string DIAG_FILENAME;
    
    /* ======================================================================= */
    /* ---- Input options from the TIMESERIES MENU --------------------------- */
    /* ======================================================================= */

    const std::string TS_FOLDER;
    const bool TS_SPEC;
    const std::string TS_SPEC_FILEPATH;
    const std::vector<int> TS_SPEC_LIST;
    const double TS_FREQ;

    const bool TS_AERO;
    const std::string TS_AERO_FILEPATH;
    const std::vector<int> TS_AERO_LIST;
    const double TS_AERO_FREQ;

    /* ======================================================================= */
    /* ---- Input options from the PROD & LOSS MENU -------------------------- */
    /* ======================================================================= */

    /* TODO: Implement PL rates */
    const bool SAVE_PL;
    const bool SAVE_O3PL;

    /* If DIAG_OUTPUT is turned on, make sure that output timestep is a
     * multiple of the dynamic timestep */

    double temperature_K;
    double pressure_Pa;
    double relHumidity_w;

    /* Compute relative humidity w.r.t ice */
    double relHumidity_i;

    MPMSimVarsWrapper() = default;
    MPMSimVarsWrapper(const Input& input, const OptInput& Input_mod);
    void initMicrophysicsVars(const Solution& data, const AIM::Aerosol& epm_LA, double iceDen);
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
