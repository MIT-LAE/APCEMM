# Source files that need to be compiled
set(SRCS
    odeSolver.cpp
    Integrate.cpp
    RHS.cpp
    )

# This command ensures the static library gets build
add_library(EPM STATIC ${SRCS})

# This command defines the dependencies of libEPM.a
target_link_libraries(EPM AIM Core Util)
