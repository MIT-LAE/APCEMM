# Source files that need to be compiled
set(SRCS
    Error.cpp
    MC_Rand.cpp
    PhysFunction.cpp
    MetFunction.cpp
    PlumeModelUtils.cpp
    VectorUtils.cpp
)

# This command ensures the static library gets build
add_library(Util STATIC ${SRCS})