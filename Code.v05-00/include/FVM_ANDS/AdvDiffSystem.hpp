#ifndef FVM_ANDS_ADVDIFFSYSTEM_H
#define FVM_ANDS_ADVDIFFSYSTEM_H

#include <Eigen/Sparse>
#include "FVM_ANDS/BoundaryCondition.hpp"
#include "FVM_ANDS_HelperFunctions.hpp"
#include <memory>

namespace FVM_ANDS{

    // Separate the SOR solver for testing without having to build an AdvDiffSystem object
    void sor_solve(const Eigen::SparseMatrix<double, Eigen::RowMajor> &A, const Eigen::VectorXd &rhs, Eigen::VectorXd &phi, double omega = 1.0, double threshold = 1e-3, int n_iters = 3);

    /**
     * @brief 1D Flux-Form Semi-Lagrangian (FFSL) advection with Lax-Wendroff TVD subgrid reconstruction.
     *
     * @details
     * Solves the 1D linear advection equation for a scalar field \f$\phi(s, t)\f$:
     * \f[
     *   \frac{\partial \phi}{\partial t} + v \frac{\partial \phi}{\partial s} = 0
     * \f]
     * across a uniform grid with spacing \f$\Delta s\f$ and arbitrary timestep \f$\Delta t\f$ without CFL restrictions.
     *
     * ### Method Overview
     *
     * The method combines two key concepts from atmospheric transport and hyperbolic conservation laws:
     *
     * 1. **Integer-Shift Trajectory Decomposition (Ritchie, 1986)**:
     *    The advective displacement \f$\Delta s_{\text{total}} = v \Delta t\f$ is decomposed into an integer
     *    grid-cell translation \f$k\f$ and a subgrid residual displacement \f$\delta s\f$:
     *    \f[
     *      k = \left\lfloor \frac{|v| \Delta t}{\Delta s} \right\rfloor \in \mathbb{Z}_{\ge 0}, \quad
     *      \delta s = |v| \Delta t - k \Delta s \in [0, \Delta s)
     *    \f]
     *    The fractional Courant number is \f$c_{\text{frac}} = \frac{\delta s}{\Delta s} \in [0, 1)\f$,
     *    and the residual timestep is \f$\Delta t_{\text{rem}} = \frac{\delta s}{|v|}\f$.
     *    - For \f$v > 0\f$: cells are shifted rightward by \f$k\f$ positions (\f$\phi_m \leftarrow \phi_{m-k}\f$),
     *      with inflow \f$m \in [0, k-1]\f$ padded by \f$\text{bc\_left}\f$.
     *    - For \f$v < 0\f$: cells are shifted leftward by \f$k\f$ positions (\f$\phi_m \leftarrow \phi_{m+k}\f$),
     *      with inflow \f$m \in [N-k, N-1]\f$ padded by \f$\text{bc\_right}\f$.
     *
     * 2. **Flux-Form Subgrid Advection with TVD Limiter (Lin & Rood, 1996; LeVeque, 2002)**:
     *    The remaining subgrid displacement is evolved via a single-step conservative finite-volume update:
     *    \f[
     *      \phi_m^{n+1} = \phi_m^n - c_{\text{frac}} \left( F_{m+1/2} - F_{m-1/2} \right)
     *    \f]
     *    To maintain second-order accuracy in time and total-variation-diminishing (TVD) monotonicity,
     *    numerical interface fluxes \f$F_{m+1/2}\f$ represent the time-average of the characteristic departure
     *    interval \f$[x_{m+1/2} - v \Delta t_{\text{rem}}, \, x_{m+1/2}]\f$. As derived in LeVeque (2002, Ch. 6),
     *    evaluating the linear reconstruction at the centroid of this interval introduces the Lax-Wendroff factor
     *    \f$(1 - c_{\text{frac}})\f$:
     *    \f[
     *      F_{m+1/2} = 
     *      \begin{cases}
     *        \phi_m + \frac{1}{2} (1 - c_{\text{frac}}) \, \text{minmod}(\Delta \phi_{m-1/2}, \Delta \phi_{m+1/2}) & \text{if } v > 0 \\
     *        \phi_{m+1} - \frac{1}{2} (1 - c_{\text{frac}}) \, \text{minmod}(\Delta \phi_{m+1/2}, \Delta \phi_{m+3/2}) & \text{if } v < 0
     *      \end{cases}
     *    \f]
     *    where \f$\text{minmod}(a, b) = \text{sgn}(a) \max\left(0, \min(|a|, b \cdot \text{sgn}(a))\right)\f$.
     *
     * ### Key Properties
     * - **Strict Mass Conservation**: Guaranteed by the conservative flux-differencing formulation (Lin & Rood, 1996).
     * - **Monotonicity (TVD)**: The minmod limiter with the \f$(1 - c_{\text{frac}})\f$ Lax-Wendroff correction
     *   prevents spurious numerical oscillations for all Courant numbers.
     * - **Unconditional Stability**: The integer translation ensures the residual Eulerian step always satisfies
     *   \f$c_{\text{frac}} < 1\f$.
     *
     * ### References
     * - Ritchie, H. (1986). Eliminating the interpolation associated with the semi-Lagrangian scheme.
     *   *Monthly Weather Review*, 114(1), 135–146.
     * - Lin, S.-J., & Rood, R. B. (1996). Multidimensional flux-form semi-Lagrangian transport schemes.
     *   *Monthly Weather Review*, 124(9), 2046–2070.
     * - LeVeque, R. J. (2002). *Finite Volume Methods for Hyperbolic Problems*. Cambridge University Press,
     *   Chapters 6 (TVD Limiters) & 9 (Variable-Coefficient and Large Time-Step Methods).
     *
     * @param[in,out] slice     1D vector of cell-centered scalar values across the slice (modified in place).
     * @param[in]     velocity  Advecting velocity (\f$v\f$) along the coordinate direction [m/s].
     * @param[in]     dt        Advective timestep (\f$\Delta t\f$) [s].
     * @param[in]     ds        Grid cell spacing (\f$\Delta s\f$) [m].
     * @param[in]     bc_left   Dirichlet boundary value at the left (inflow/outflow) face.
     * @param[in]     bc_right  Dirichlet boundary value at the right (inflow/outflow) face.
     */
    void semiLagrangianAdvection1D(std::vector<double>& slice, double velocity, double dt, double ds, double bc_left, double bc_right);

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
            void semiLagrangianAdvection(double dt, bool parallelAdvection = false);
            Eigen::VectorXd forwardEulerAdvection(bool operatorSplit = false, bool parallelAdvection = false) const noexcept;
            // Breakup the implementation of sor_solve to allow for easy testing by inputing an arbitrary linear system to solve:
            // Implementation is moved outside of the class, and make class method to be used in code
            void sor_solve(double omega = 1.0, double threshold = 1e-3, int n_iters = 3){ 
                FVM_ANDS::sor_solve(getCoefMatrix(), rhs_, phi_, omega, threshold, n_iters); 
            };
            inline const Eigen::VectorXd& getRHS() const { return rhs_; }
            inline const Eigen::VectorXd& phi() const { return phi_; }
            inline const std::vector<std::unique_ptr<Point>>& points() const { return points_; }
            inline const Eigen::SparseMatrix<double, Eigen::RowMajor>& getCoefMatrix() const { 
                // If using prebuilt matrix, this instance of AdvDiffSystem does not have its own
                // totalCoefMatrix_, instead it holds a pointer to a shared matrix held in another AdvDiffSystem
                // Getter abstracts that away to have one access pattern
                return use_shared_totalCoefMatrix_ ? *shared_totalCoefMatrixPtr_ : totalCoefMatrix_; 
            }
            inline std::shared_ptr<const Eigen::SparseMatrix<double, Eigen::RowMajor>> getCoefMatrixPtr() {
                // Only used by the AdvDiffSystem instance that actually computed the matrix to make it available to
                // other instances without copying the full matrix. Other instances will just hold a copy of
                // this shared pointer
                return std::make_shared<const Eigen::SparseMatrix<double, Eigen::RowMajor>>(totalCoefMatrix_);
            }
            inline void setCoefMatrix(std::shared_ptr<const Eigen::SparseMatrix<double, Eigen::RowMajor>> matrixPtr) {
                // Only used by AdvDiffSystem instances that did not compute the matrix. This sets the shared pointer
                // to the matrix that is held in another instance
                shared_totalCoefMatrixPtr_ = matrixPtr;
                use_shared_totalCoefMatrix_ = true;
            }
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
            inline void updateDiffusion(const Eigen::VectorXd& Dh, const Eigen::VectorXd& Dv){
                Dh_vec_ = Dh;
                Dv_vec_ = Dv;
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
            std::shared_ptr<const Eigen::SparseMatrix<double, Eigen::RowMajor>> shared_totalCoefMatrixPtr_;
            bool use_shared_totalCoefMatrix_ = false;
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
                    /* Function is never called with FaceDirection::ERROR, see call stack
                    Placeholder value instead of throwing exception because function
                    is marked as noexcept... */
                    default:
                        return -1;
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
