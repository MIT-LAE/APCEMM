# Source files that need to be compiled
set(SRCS
    Aircraft.cpp
    #BoxModel.cpp
    Cluster.cpp
    Diag_Mod.cpp
    Emission.cpp
    Engine.cpp
    Fuel.cpp
    Input_Mod.cpp
    Input.cpp
    LAGRIDPlumeModel.cpp
    LiquidAer.cpp
    Meteorology.cpp
    Mesh.cpp
    MPMSimVarsWrapper.cpp
    PlumeModel.cpp
    ReadJRates.cpp
    Ring.cpp
    #Save.cpp
    Species.cpp
    Structure.cpp
    SZA.cpp
    Util.cpp
    TimestepVarsWrapper.cpp
    Vortex.cpp
    )

# This command ensures the static library gets build
add_library(Core STATIC ${SRCS})

# This command defines the dependencies of libCore.a
target_link_libraries(Core ${netCDFCxx_LIBRARIES} ${NETCDF_LIBRARIES} ${YAML_CPP_LIBRARIES}
                     FVM_ANDS AIM Util EPM KPP YamlInputReader)
