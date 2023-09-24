#include <Core/Input.hpp>
#include <Core/Input_Mod.hpp>
#ifndef TIMESTEPVARSWRAPPER_H_INCLUDED
#define TIMESTEPVARSWRAPPER_H_INCLUDED
struct TimestepVarsWrapper
{

    // DO NOT change the order of the variables declared here.
    // It may cause undefined behavior in the constructor.
    // - MX
    /*
     *  - tEmission is the local emission time expressed in hours
     *  (between 0.0 and 24.0)
     *  - tInitial is the local time at which the simulation starts in hours
     *  - simulationTime represents the simulation time (in hours) (now read from
     *    input file)
     *  - tFinal corresponds to the final time of the simulation expressed in hours
     */

    /* Define emission and simulation time */
    const double tEmission_h;
    const double tInitial_h;
    const double tFinal_h;
    const double tInitial_s;
    const double tFinal_s;

    /* Current time in [s] */
    double curr_Time_s; /* [s] */
    /* Time step in [s] */
    double dt; /* [s] */

    /* Create time array */

    /* Vector of time in [s] */
    std::vector<double> timeArray;

    /* Time counter [-] */
    int nTime;

    bool LAST_STEP;
    bool ITS_TIME_FOR_TRANSPORT;
    bool ITS_TIME_FOR_CHEM;
    bool ITS_TIME_FOR_LIQ_COAGULATION;
    bool ITS_TIME_FOR_ICE_COAGULATION;
    bool ITS_TIME_FOR_ICE_GROWTH;
    bool ITS_TIME_FOR_TEMPPERTURB;

    const double TRANSPORT_DT;    // [s]
    const double CHEMISTRY_DT;    // [s]
    const double COAG_DT;         //[s]
    const double TEMP_PERTURB_DT; //[s]
    const double ICE_GROWTH_DT;   // [s]

    double lastTimeTransport;
    double lastTimeChem;
    double lastTimeLiqCoag;
    double lastTimeIceCoag;
    double lastTimeIceGrowth;
    double lastTimeTempPerturb;

    double totalIceParticles_before;
    double totalIceMass_before;
    double totalIceParticles_initial;
    double totalIceMass_initial;
    double totalIceParticles_now;
    double totalIceMass_now;
    double totalIceParticles_last;
    double totalIceMass_last;
    double totalIceParticles_after;
    double totalIceMass_after;
    double totPart_lost;
    double totIce_lost;
    const double ABORT_THRESHOLD;

    TimestepVarsWrapper() = default;
    TimestepVarsWrapper(const Input &input, const OptInput &Input_Opt);
    inline void setTimeArray(const Vector_1D &vec)
    {
        timeArray = vec;
    }
    inline bool checkLastStep()
    {
        LAST_STEP = (curr_Time_s + dt >= tFinal_s);
        return LAST_STEP;
    }
    inline bool checkTimeForTransport()
    {
        ITS_TIME_FOR_TRANSPORT = (((curr_Time_s + dt - lastTimeTransport) >= TRANSPORT_DT) || LAST_STEP);
        return ITS_TIME_FOR_TRANSPORT;
    }
    inline bool checkTimeForChem()
    {
        ITS_TIME_FOR_CHEM = (((curr_Time_s + dt - lastTimeChem) >= CHEMISTRY_DT) || LAST_STEP);
        return ITS_TIME_FOR_CHEM;
    }
    inline bool checkTimeForTempPerturb()
    {
        ITS_TIME_FOR_TEMPPERTURB = (((curr_Time_s + dt - lastTimeTempPerturb) >= TEMP_PERTURB_DT) || LAST_STEP);
        return ITS_TIME_FOR_TEMPPERTURB;
    }
    inline bool checkTimeForLiqCoag()
    {
        ITS_TIME_FOR_LIQ_COAGULATION = (((curr_Time_s + dt - lastTimeLiqCoag) >= COAG_DT) || LAST_STEP);
        return ITS_TIME_FOR_LIQ_COAGULATION;
    }
    inline bool checkTimeForIceCoag()
    {
        ITS_TIME_FOR_ICE_COAGULATION = (((curr_Time_s + dt - lastTimeIceCoag) >= COAG_DT) || LAST_STEP);
        return ITS_TIME_FOR_ICE_COAGULATION;
    }
    inline bool checkTimeForIceGrowth()
    {
        /* TODO: For now perform growth at every time step */
        ITS_TIME_FOR_ICE_GROWTH = (((curr_Time_s + dt - lastTimeIceGrowth) >= ICE_GROWTH_DT) || LAST_STEP);
        return ITS_TIME_FOR_ICE_GROWTH;
    }
};

#endif