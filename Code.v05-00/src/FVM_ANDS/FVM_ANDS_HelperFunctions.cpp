#include <FVM_ANDS/FVM_ANDS_HelperFunctions.hpp>
namespace FVM_ANDS{
    int twoDIdx_to_vecIdx(int idx_x, int idx_y, int nx, int ny, vecFormat format){
        return (format == vecFormat::ROWMAJOR) ?
                idx_y * nx + idx_x :
                idx_x * ny + idx_y;
    }
    Eigen::VectorXd std2dVec_to_eigenVec(const Vector_2D& phi, vecFormat format){
        //Doesn't create ghost nodes but its fine because AdvDiffSys constructor resizes.
        int ny = phi.size();
        int nx = ny > 0 ? phi[0].size() : 0;
        Eigen::VectorXd vec(nx*ny);
        for(int i = 0; i < nx; i++){
            for(int j = 0; j < ny; j++){
                int vecIdx = twoDIdx_to_vecIdx(i, j, nx, ny, format);
                vec[vecIdx] = phi[j][i];
            }
        }
        return vec;
    }
    Vector_2D eigenVec_to_std2dVec(Eigen::VectorXd eig_vec, int nx, int ny){
        Vector_2D std2dVec(ny, Vector_1D(nx, 0));
        for(int i = 0; i < nx; i++){
            for(int j = 0; j < ny; j++){
                std2dVec[j][i] = eig_vec(ny*i + j);
            }
        }
        return std2dVec;    
    }
    BoundaryConditions bcFrom2DVector(const Vector_2D& initialVec, bool zeroBC){
        int ny = initialVec.size();
        int nx = initialVec[0].size();
        BoundaryConditions bc;
        bc.bcType_top = BoundaryConditionFlag::DIRICHLET_INT_BPOINT;
        bc.bcType_left = BoundaryConditionFlag::DIRICHLET_INT_BPOINT;
        bc.bcType_right = BoundaryConditionFlag::DIRICHLET_INT_BPOINT;
        bc.bcType_bot = BoundaryConditionFlag::DIRICHLET_INT_BPOINT;
        bc.bcVals_top = Vector_1D(nx, 0);
        bc.bcVals_bot = Vector_1D(nx, 0);
        bc.bcVals_left = Vector_1D(ny, 0);
        bc.bcVals_right = Vector_1D(ny, 0);

        if(!zeroBC){
            for(int i = 0; i < nx; i++){
                bc.bcVals_top[i] = initialVec[ny-1][i];
                bc.bcVals_bot[i] = initialVec[0][i];
            }
            for(int j = 0; j < ny; j++){
                bc.bcVals_left[j] = initialVec[j][0];
                bc.bcVals_right[j] = initialVec[j][nx-1];
            }
        }
        return bc;
    }

}