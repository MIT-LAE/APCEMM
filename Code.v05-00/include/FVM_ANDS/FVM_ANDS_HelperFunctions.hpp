#ifndef FVM_ANDS_HELPERFUNCTIONS_H
#define FVM_ANDS_HELPERFUNCTIONS_H

#include <Eigen/Sparse>
#include <Util/ForwardDecl.hpp>
#include <Util/PhysConstant.hpp>
#include "FVM_ANDS/BoundaryCondition.hpp"

namespace FVM_ANDS{
    int twoDIdx_to_vecIdx(int idx_x, int idx_y, int nx, int ny, vecFormat format = vecFormat::COLMAJOR);
    Eigen::VectorXd std2dVec_to_eigenVec(const Vector_2D& phi, vecFormat format = vecFormat::COLMAJOR);
    BoundaryConditions bcFrom2DVector(const Vector_2D& initialVec, bool zeroBC = false);
    Vector_2D eigenVec_to_std2dVec(Eigen::VectorXd eig_vec, int nx, int ny);
} 
#endif