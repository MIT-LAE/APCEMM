#include "FVM_ANDS/BoundaryCondition.hpp"

namespace FVM_ANDS {
    Point::Point(BoundaryConditionFlag bcType){
        bc_ = BoundaryCondDescription(bcType);
    }
    Point::Point(BoundaryCondDescription bcDesc){
        bc_ = bcDesc;
    }

    GhostPoint::GhostPoint(BoundaryCondDescription bc) {
        const bool is_ghost_flag = static_cast<int>(bc.bcType) > 100;
        if(!is_ghost_flag){
            bc.bcType = static_cast<BoundaryConditionFlag>(static_cast<int>(bc.bcType) + 100);
        }
        setBC(bc);
    }
    IntBoundPoint::IntBoundPoint(BoundaryCondDescription bc1)
     : Point(bc1) {
        secondary_bc_= std::nullopt;
    }
    IntBoundPoint::IntBoundPoint(BoundaryCondDescription bc1, BoundaryCondDescription bc2)
     : Point(bc1) {
        secondary_bc_ = bc2;
    }

    BoundaryCondDescription::BoundaryCondDescription()
    :   bcType(BoundaryConditionFlag::ERROR),
        direction(FaceDirection::ERROR),
        bcVal(-1),
        corrPoint(-1){

    }
    BoundaryCondDescription::BoundaryCondDescription(BoundaryConditionFlag bcType)
    :   bcType(bcType),
        direction(FaceDirection::ERROR),
        bcVal(-1),
        corrPoint(-1){

    }
    BoundaryCondDescription::BoundaryCondDescription(BoundaryConditionFlag bcType, FaceDirection direction, double bcVal, int corrPoint)
    :   bcType(bcType),
        direction (direction),
        bcVal (bcVal),
        corrPoint(corrPoint){
        
    }
}