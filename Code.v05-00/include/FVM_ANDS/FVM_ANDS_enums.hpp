#include <limits>
#ifndef FVM_ANDS_ENUMS_H
#define FVM_ANDS_ENUMS_H
namespace FVM_ANDS{
    enum class BoundaryConditionFlag: unsigned char {
        INTERIOR = 0,
        DIRICHLET_INT_BPOINT = 1,
        PERIODIC_INT_BPOINT = 2,
        DIRICHLET_GHOSTPOINT = 101,
        PERIODIC_GHOSTPOINT = 102,
        ERROR = std::numeric_limits<unsigned char>::max()
    };
    enum class FaceDirection: unsigned char {
        NORTH,
        SOUTH,
        EAST,
        WEST,
        ERROR
    };
    enum class SolutionMethod : unsigned char {
        ForwardEuler,
        BackwardEuler,
        OperatorSplitting
    };
    enum class AdvectionScheme : unsigned char {
        FirstOrderUpwind,
        CentralDifference,
        MinMod
    };
    enum class vecFormat: unsigned char {
        ROWMAJOR,
        COLMAJOR
    };
}
#endif