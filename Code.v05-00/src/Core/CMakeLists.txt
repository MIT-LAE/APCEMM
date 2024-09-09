# Source files that need to be compiled
set(SRCS
    Aircraft.cpp
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
    Ring.cpp
    Species.cpp
    Structure.cpp
    SZA.cpp
    Util.cpp
    TimestepVarsWrapper.cpp
    Vortex.cpp
    )

# This command ensures the static library gets build
add_library(Core STATIC ${SRCS})

target_link_libraries(Core PRIVATE netCDF::netcdf netCDF::netcdf-cxx4)
target_link_libraries(Core PRIVATE yaml-cpp::yaml-cpp)

# This command defines the dependencies of libCore.a
target_link_libraries(Core PRIVATE FVM_ANDS AIM Util EPM KPP YamlInputReader)
