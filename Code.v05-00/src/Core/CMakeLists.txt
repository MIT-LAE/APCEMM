# Source files that need to be compiled
set(SRCS
    Advection.cpp
    Aircraft.cpp
    Ambient.cpp
    BoxModel.cpp
    BuildTime.cpp
    Cluster.cpp
    Diag_Mod.cpp
    Diffusion.cpp
    Emission.cpp
    Engine.cpp
    Fuel.cpp
    Input_Mod.cpp
    Input.cpp
    LiquidAer.cpp
    Meteorology.cpp
    Mesh.cpp
    PlumeModel.cpp
    ReadJRates.cpp
    Ring.cpp
    Save.cpp
    Species.cpp
    Structure.cpp
    SZA.cpp
    Timer.cpp
    Util.cpp
    Vortex.cpp
    )

# This command ensures the static library gets build
add_library(Core STATIC ${SRCS})

# This command defines the dependencies of libCore.a
target_link_libraries(Core ${netCDFCxx_LIBRARIES} ${NETCDF_LIBRARIES} ${YAML_CPP_LIBRARIES}
                     AIM Util EPM KPP SANDS YamlInputReader)