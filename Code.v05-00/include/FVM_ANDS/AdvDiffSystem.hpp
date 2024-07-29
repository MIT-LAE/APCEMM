#ifndef FVM_ANDS_ADVDIFFSYSTEM_H
#define FVM_ANDS_ADVDIFFSYSTEM_H

#include <Eigen/Sparse>
#include <unordered_map>
#include "Core/Mesh.hpp"
#include <chrono>
#include "FVM_ANDS/BoundaryCondition.hpp"
#include "FVM_ANDS_HelperFunctions.hpp"
#include <memory>

namespace FVM_ANDS{

    // Separate the SOR solver for testing without having to build an AdvDiffSystem object
    void sor_solve(const Eigen::SparseMatrix<double, Eigen::RowMajor> &A, const Eigen::VectorXd &rhs, Eigen::VectorXd &phi, double omega = 1.0, double threshold = 1e-3, int n_iters = 3);

    struct AdvDiffParams {
        AdvDiffParams(double u, double v, double shear, double Dh, double Dv, double dt){
            this->u = u;
            this->v = v;
            this->shear = shear;
            this->Dh = Dh;
            this->Dv = Dv;
            this->dt = dt;
        }
        double u;
        double v;
        double shear;
        double Dh;
        double Dv;
        double dt;
    };
    class AdvDiffSystem{
        public:
            AdvDiffSystem() = delete;
            AdvDiffSystem(const AdvDiffParams& params, const Vector_1D xCoords, const Vector_1D yCoords, const BoundaryConditions& bc, const Eigen::VectorXd& phi_init, vecFormat format = vecFormat::COLMAJOR);
            void buildCoeffMatrix(bool operatorSplit = false);
            const Eigen::VectorXd& calcRHS();
            void applyBoundaryCondition();
            void updateBoundaryCondition(const BoundaryConditions& bc);
            Eigen::VectorXd forwardEulerAdvection(bool operatorSplit = false, bool parallelAdvection = false) const noexcept;
            // Breakup the implementation of sor_solve to allow for easy testing by inputing an arbitrary linear system to solve:
            // Implementation is moved outside of the class, and make class method to be used in code
            void sor_solve(double omega = 1.0, double threshold = 1e-3, int n_iters = 3){ FVM_ANDS::sor_solve(totalCoefMatrix_, rhs_, phi_, omega, threshold, n_iters); };
            inline const Eigen::VectorXd& getRHS() const { return rhs_; }
            inline const Eigen::VectorXd& phi() const { return phi_; }
            inline const std::vector<std::unique_ptr<Point>>& points() const { return points_; }
            inline const Eigen::SparseMatrix<double, Eigen::RowMajor>& getCoefMatrix() const { return totalCoefMatrix_; }
            inline void updatePhi(const Eigen::VectorXd& phi_new){ 
                //Need to resize to account for grid changing in size.
                phi_.resize(nx_ * ny_ + 2*nx_ + 2*ny_);
                phi_(Eigen::seq(0, nx_ * ny_ - 1)) = phi_new(Eigen::seq(0, nx_ * ny_ - 1));
            }
            inline void addSource(const Eigen::VectorXd& source){ source_ = source; }
            inline void updateDiffusion(double Dh, double Dv){
                for(int i = 0; i < nx_; i++){
                    for(int j = 0; j < ny_; j++){
                        int vector_idx = twoDIdx_to_vecIdx(i, j, nx_, ny_, format_);
                        Dh_vec_[vector_idx] = Dh;
                        Dv_vec_[vector_idx] = Dv;
                    }
                }
            }
            inline void updateDiffusion(const Vector_2D& Dh, const Vector_2D& Dv){
                Dh_vec_ = std2dVec_to_eigenVec(Dh);
                Dv_vec_ = std2dVec_to_eigenVec(Dv);
            }
            inline void updateAdvection(double u, double v, double shear){
                u_double_ = u;
                v_double_ = v;
                shear_ = shear;
                initVelocVecs();
            }
            inline double timestep() const { return dt_; }
            inline void updateDy(double dy_new) { 
                dy_ = dy_new;
                invdy_ = dy_new;  
            }
            inline void updateDx(double dx_new) { 
                dx_ = dx_new;
                invdx_ = dx_new;  
            }
            inline void updateYCoord(const Vector_1D& yCoord_new) { 
                yCoord_ = yCoord_new;
            }
            inline void updateNx(int nx_new) { nx_ = nx_new; }
            inline void updateNy(int ny_new) { ny_ = ny_new; }

            inline void updateSpacing(const Vector_1D& yCoord_new, double dx_new, int nx_new) {
                updateYCoord(yCoord_new);
                updateDy(yCoord_new[1] - yCoord_new[0]);
                updateDy(dx_new);
                updateNy(yCoord_new.size());
                updateNx(nx_new);
            }
            inline void updateTimestep(double dt){ dt_ = dt; }
            inline double courant() const{
                auto maxCoeffAbsolute = [] (const Eigen::VectorXd& vec) -> double {
                    return std::max(vec.maxCoeff(), std::abs(vec.minCoeff()));
                };
                return maxCoeffAbsolute(u_vec_) * dt_ / dx_ + maxCoeffAbsolute(v_vec_) * dt_ / dy_;
            }
            inline void scaleRHS(double scalingFactor) { rhs_ = rhs_ * scalingFactor; }
            inline void scalePhi(double scalingFactor) { phi_ = phi_ * scalingFactor; }
            void scaleBC(double scalingFactor) {
                for (int i = 0; i < nx_; i++){
                    bcVals_top_[i] *= scalingFactor;
                    bcVals_bot_[i] *= scalingFactor;
                }
                for (int i = 0; i < ny_; i++){
                    bcVals_left_[i] *= scalingFactor;
                    bcVals_right_[i] *= scalingFactor;
                }
            }
            

        private:

            vecFormat format_;
            Eigen::VectorXd u_vec_;
            Eigen::VectorXd v_vec_;
            double u_double_;
            double v_double_;
            double shear_;
            Eigen::VectorXd Dh_vec_;
            Eigen::VectorXd Dv_vec_;
            double dt_;
            double dx_;
            double dy_;
            double invdx_;
            double invdy_;
            int nx_;
            int ny_;
            int nInteriorPoints_;
            int nGhostPoints_;
            int nTotalPoints_;
            Vector_1D yCoord_;
            BoundaryConditionFlag bcType_top_;
            BoundaryConditionFlag bcType_left_;
            BoundaryConditionFlag bcType_right_;
            BoundaryConditionFlag bcType_bot_;
            Vector_1D bcVals_top_;
            Vector_1D bcVals_left_;
            Vector_1D bcVals_right_; 
            Vector_1D bcVals_bot_;
            std::vector<std::unique_ptr<Point>> points_;
            Eigen::SparseMatrix<double, Eigen::RowMajor> totalCoefMatrix_;
            Eigen::VectorXd rhs_;
            Eigen::VectorXd phi_;
            Eigen::VectorXd source_;
            Eigen::VectorXd deferredCorr_;

            void initVelocVecs();
            void buildPointList();
            void buildAdvectionCoeffs(int i, double& coeff_C, double& coeff_N, double& coeff_S, double& coeff_E, double& coeff_W);
            void updateGhostNodes();

            inline bool isValidPointID(int idx) const {
                return (idx >= 0 && idx < phi_.rows());
            }

            inline double minmod(int pointID, FaceDirection face, double faceVelocity) const noexcept {
                //INDEXING ASSUMES ROW-MAJOR 2DVEC -> 1DVEC MAPPING
                //Min-mod flux limiter function to ensure TVD / monoticity
                //r = upwind gradient ratio
                //psi = flux limiter function
                //See Versteeg & Malasekera (2007) Ch. 5.10 for specifics on notation

                double phi_P = phi_[pointID];
                double r = 0; //shut up valgrind
                switch(face){
                    case FaceDirection::NORTH:{
                        double phi_N = phi_[neighbor_point(FaceDirection::NORTH, pointID)];
                        if(phi_N - phi_P == 0) return 0;
                        if(faceVelocity >= 0){
                            double phi_S = phi_[neighbor_point(FaceDirection::SOUTH, pointID)];
                            r = (phi_P - phi_S) / (phi_N - phi_P);
                            break;
                        }
                        else {
                            int pointID_N = neighbor_point(FaceDirection::NORTH, pointID);
                            double phi_NN = phi_[neighbor_point(FaceDirection::NORTH, pointID_N)];
                            r = (phi_NN - phi_N) / (phi_N - phi_P);
                            break;
                        }
                    }

                    case FaceDirection::SOUTH:{
                        double phi_S = phi_[neighbor_point(FaceDirection::SOUTH, pointID)];
                        if(phi_P - phi_S == 0) return 0;
                        if(faceVelocity >= 0){
                            int pointID_S = neighbor_point(FaceDirection::SOUTH, pointID);
                            double phi_SS = phi_[neighbor_point(FaceDirection::SOUTH, pointID_S)];
                            r = (phi_S - phi_SS) / (phi_P - phi_S);
                            break;
                        }
                        else {
                            double phi_N = phi_[neighbor_point(FaceDirection::NORTH, pointID)];
                            r = (phi_N - phi_P) / (phi_P - phi_S);
                            break;
                        }
                    }

                    case FaceDirection::EAST:{
                        double phi_E = phi_[neighbor_point(FaceDirection::EAST, pointID)];
                        if (phi_E - phi_P == 0) return 0;                    
                        if(faceVelocity >= 0){
                            double phi_W = phi_[neighbor_point(FaceDirection::WEST, pointID)];
                            r = (phi_P - phi_W) / (phi_E - phi_P);
                            break;
                        }
                        else {
                            int pointID_E = neighbor_point(FaceDirection::EAST, pointID);
                            double phi_EE = phi_[neighbor_point(FaceDirection::EAST, pointID_E)];
                            r = (phi_EE - phi_E) / (phi_E - phi_P);
                            break;
                        }
                    }
                    
                    case FaceDirection::WEST:{
                        double phi_W = phi_[neighbor_point(FaceDirection::WEST, pointID)];
                        if(phi_P - phi_W == 0) return 0;
                        if(faceVelocity >= 0){
                            int pointID_W = neighbor_point(FaceDirection::WEST, pointID);
                            double phi_WW = phi_[neighbor_point(FaceDirection::WEST, pointID_W)];
                            r = (phi_W - phi_WW) / (phi_P - phi_W);
                            break;
                        }
                        else {
                            double phi_E = phi_[neighbor_point(FaceDirection::EAST, pointID)];
                            r = (phi_E - phi_P) / (phi_P - phi_W);
                            break;
                        }
                    }
                    case FaceDirection::ERROR:{
                        throw std::runtime_error("Received invalid face direction");
                    }
                }
                return std::max(0.0, std::min(r, 1.0));
            }
            inline double minmod_N_vPos(int pointID) const noexcept{
                if(!isValidPointID(pointID + 1) || !isValidPointID(pointID - 1)) return 0;
                double phi_P = phi_[pointID];
                double phi_N = phi_[pointID + 1];
                double phi_S = phi_[pointID - 1];
                double r = (phi_N - phi_P == 0) ? 0 : (phi_P - phi_S) / (phi_N - phi_P);
                return std::max(0.0, std::min(r, 1.0));
            }
            inline double minmod_N_vNeg(int pointID) const noexcept{
                if(!isValidPointID(pointID + 2)) return 0;
                double phi_P = phi_[pointID];
                double phi_N = phi_[pointID + 1];
                double phi_NN = phi_[pointID + 2];
                double r = (phi_N - phi_P == 0) ? 0 : (phi_NN - phi_N) / (phi_N - phi_P);
                return std::max(0.0, std::min(r, 1.0));
            }
            inline double minmod_S_vPos(int pointID) const noexcept{
                if(!isValidPointID(pointID - 2)) return 0;
                double phi_P = phi_[pointID];
                double phi_S = phi_[pointID - 1];
                double phi_SS = phi_[pointID - 2];
                double r = (phi_P - phi_S == 0) ? 0 : (phi_S - phi_SS) / (phi_P - phi_S);
                return std::max(0.0, std::min(r, 1.0));
            }
            inline double minmod_S_vNeg(int pointID) const noexcept{
                if(!isValidPointID(pointID - 1) || neighbor_point(FaceDirection::NORTH, pointID)) return 0;
                double phi_P = phi_[pointID];
                double phi_S = phi_[pointID - 1];
                double phi_N = phi_[neighbor_point(FaceDirection::NORTH, pointID)];
                double r = (phi_P - phi_S == 0) ? 0 : (phi_N - phi_P) / (phi_P - phi_S);
                return std::max(0.0, std::min(r, 1.0));
            }
            inline double minmod_E_vPos(int pointID) const noexcept{
                if(!isValidPointID(pointID + ny_) || !isValidPointID(pointID - ny_)) return 0;
                double phi_P = phi_[pointID];
                double phi_E = phi_[pointID + ny_];
                double phi_W = phi_[pointID - ny_];
                double r = (phi_E - phi_P == 0) ? 0 : (phi_P - phi_W) / (phi_E - phi_P);
                return std::max(0.0, std::min(r, 1.0));
            }

            inline double minmod_E_vNeg(int pointID) const noexcept{
                if(!isValidPointID(pointID + 2*ny_)) return 0;
                double phi_P = phi_[pointID];
                double phi_E = phi_[pointID + ny_];
                if (phi_E - phi_P == 0) return 0;                    
                double phi_EE = phi_[pointID + ny_ + ny_];
                double r = (phi_EE - phi_E) / (phi_E - phi_P);
                return std::max(0.0, std::min(r, 1.0));
            }
            inline double minmod_W_vPos(int pointID) const noexcept{
                if(!isValidPointID(pointID - 2*ny_)) return 0;
                double phi_P = phi_[pointID];
                double phi_W = phi_[pointID - ny_];
                if(phi_P - phi_W == 0) return 0;    
                double phi_WW = phi_[pointID - ny_ - ny_];
                double r = (phi_W - phi_WW) / (phi_P - phi_W);
                return std::max(0.0, std::min(r, 1.0));
            }
            inline double minmod_W_vNeg(int pointID) const noexcept{
                if(!isValidPointID(pointID - ny_) || !isValidPointID(pointID + ny_)) return 0;
                double phi_P = phi_[pointID];
                double phi_W = phi_[pointID - ny_];
                if(phi_P - phi_W == 0) return 0;    
                double phi_E = phi_[pointID + ny_];
                double r = (phi_E - phi_P) / (phi_P - phi_W);
                return std::max(0.0, std::min(r, 1.0));
            }
            inline int neighbor_point(FaceDirection direction, int pointID) const noexcept{
                Point* point = points_[pointID].get();
                if(point->bcType() == BoundaryConditionFlag::INTERIOR) return neighbor_point_interior(direction, pointID);
                
                if(point->bcDirection() == direction){
                    return point->corrPoint();
                }
                else if (point->secondBoundaryConds() && point->secondBoundaryConds().value().direction == direction){
                    return point->secondBoundaryConds().value().corrPoint;
                }
                return neighbor_point_interior(direction, pointID);        
            }
            inline int neighbor_point_interior(FaceDirection direction, int pointID) const noexcept{
                switch(direction){
                    case FaceDirection::NORTH:
                        return pointID + 1;

                    case FaceDirection::SOUTH:
                        return pointID - 1;
                    
                    case FaceDirection::EAST:
                        return pointID + ny_;
                    
                    case FaceDirection::WEST:
                        return pointID - ny_;
                    
                    // For consistency with previous implementation
                    default:
                        return -1;
                }
            }
    };
}
#endif
