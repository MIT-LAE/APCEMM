// Enum to monitor the return status of APCEMM
#ifndef STATUS_H_INCLUDED
#define STATUS_H_INCLUDED

enum class SimStatus : int {
    // Exit status: Contrail formed and disappeared before APCEMM reached max simulation timestep
    Complete,
    // Exit status: Contrail formed but did not disappear before APCEMM reached max simulation timestep
    Incomplete,
    // Exit status: No contrail formation due to not satisfying SAC
    NoWaterSaturation,
    // Exit status: Contrail forms but does not persist because RHi < 1
    // This status may need more thought if APCEMM is used for plume chemistry
    // regardless of contrail formation
    NoPersistence,
    // Exit status: Ice crystals do not survive vortex sinking
    NoSurvivalVortex,
    // Exit status: Simulation failed, placeholder to replace iERR < 0. In practice, iERR was never 
    // assigned a negative (error) value in Main.cpp
    Failed,

    // Intermediary status: EPM completed and contrail forms, simulation can proceed to vortex parametrization
    // and then continue to transport model
    EPMSuccess,
};

#endif // STATUS_H_INCLUDED
