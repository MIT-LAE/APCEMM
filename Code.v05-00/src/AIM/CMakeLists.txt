# Source files that need to be compiled
set(SRCS
    Aerosol.cpp
    buildKernel.cpp
    Coagulation.cpp
    Nucleation.cpp
    Settling.cpp)

# This command ensures the static library gets build
add_library(AIM STATIC ${SRCS})

# This command defines the dependencies of libAIM.a
target_link_libraries(AIM Util OpenMP::OpenMP_CXX)
