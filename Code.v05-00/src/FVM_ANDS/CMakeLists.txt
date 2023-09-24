# Source files that need to be compiled
set(SRCS
    AdvDiffSystem.cpp
    BoundaryCondition.cpp
    FVM_ANDS_HelperFunctions.cpp
    FVM_Solver.cpp
    )

# This command ensures the static library gets build
add_library(FVM_ANDS STATIC ${SRCS})
# This command defines the dependencies of libFDM_ANDS.a
target_link_libraries(FVM_ANDS Util Eigen3::Eigen OpenMP::OpenMP_CXX)