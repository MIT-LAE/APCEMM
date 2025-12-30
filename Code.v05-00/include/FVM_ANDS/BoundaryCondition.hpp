#ifndef FVM_ANDS_BOUNDARYCONDITION_H
#define FVM_ANDS_BOUNDARYCONDITION_H

#include <Util/ForwardDecl.hpp>
#include <optional>
#include "FVM_ANDS_enums.hpp"
namespace FVM_ANDS{
    struct BoundaryCondDescription {
        BoundaryCondDescription();
        BoundaryCondDescription(BoundaryConditionFlag bcType);
        BoundaryCondDescription(BoundaryConditionFlag bcType, FaceDirection direction, double bcVal, int corrPoint);
        BoundaryConditionFlag bcType;
        FaceDirection direction;
        double bcVal;
        int corrPoint;
    };
    class Point {
        public:
            Point() = default;
            virtual ~Point() = default;
            Point(BoundaryConditionFlag bcType);
            Point(BoundaryCondDescription bcDesc);
            inline void setBC(BoundaryCondDescription bc){
                bc_ = bc;
            }
            inline BoundaryCondDescription boundaryConds() const {
                return bc_;
            }
            inline BoundaryConditionFlag bcType() const noexcept {
                return bc_.bcType;
            }
            inline double bcVal() const {
                return bc_.bcVal;
            }
            inline FaceDirection bcDirection() const {
                return bc_.direction;
            }
            inline int corrPoint() const {
                return bc_.corrPoint;
            }
            inline void setBCVal(double bcVal){
                bc_.bcVal = bcVal;
            }
            inline virtual void setBCType(BoundaryConditionFlag bcType){
                bc_.bcType = bcType;
            }
            constexpr virtual bool isGhost() const {
                return false;
            }
            inline virtual void setSecondaryBC(BoundaryCondDescription bc)  {
                return;
            }
            virtual inline std::optional<BoundaryCondDescription> secondBoundaryConds() const {
                return std::nullopt;
            }
        protected:
            BoundaryCondDescription bc_;
    };
    struct GhostPoint : public Point {
        public:
            GhostPoint() = delete;
            GhostPoint(BoundaryCondDescription bc);
            constexpr bool isGhost() const override {
                return true;
            }
            inline void setBCType(BoundaryConditionFlag bcType) override {
                bc_.bcType = static_cast<BoundaryConditionFlag>((static_cast<int>(bcType) < 100) * 100 + static_cast<int>(bcType));
            }
    };
    struct IntBoundPoint : public Point {
        public:
            IntBoundPoint() = delete;
            IntBoundPoint(BoundaryCondDescription bc1);
            IntBoundPoint(BoundaryCondDescription bc1, BoundaryCondDescription bc2);
            inline void setSecondaryBC(BoundaryCondDescription bc) override {
                secondary_bc_ = bc;
            }
            inline std::optional<BoundaryCondDescription> secondBoundaryConds() const override{
                return secondary_bc_;
            }

        private:
            std::optional<BoundaryCondDescription> secondary_bc_;
    };
    struct BoundaryConditions {
        BoundaryConditionFlag bcType_top;
        BoundaryConditionFlag bcType_left;
        BoundaryConditionFlag bcType_right;
        BoundaryConditionFlag bcType_bot;
        Vector_1D bcVals_top;
        Vector_1D bcVals_left;
        Vector_1D bcVals_right; 
        Vector_1D bcVals_bot;
    };
}
#endif