# Source files that need to be compiled
set(SRCS
    Solver.cpp
    FFT.cpp
    )

# This command ensures the static library gets build
add_library(SANDS STATIC ${SRCS})

# This command defines the dependencies of libSANDS.a
target_link_libraries(SANDS ${FFTW_LIBRARIES} OpenMP::OpenMP_CXX Core Util)