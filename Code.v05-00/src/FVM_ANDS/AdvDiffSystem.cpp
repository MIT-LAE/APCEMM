#include <FVM_ANDS/AdvDiffSystem.hpp>
#include <iostream>
#include <math.h>
using std::cout;
using std::endl;
namespace FVM_ANDS{
    AdvDiffSystem::AdvDiffSystem(const AdvDiffParams& params, const Vector_1D xCoords, const Vector_1D yCoords, const BoundaryConditions& bc, const Eigen::VectorXd& phi_init, vecFormat format) :
        format_(format),
        u_double_ (params.u),
        v_double_ (params.v),
        shear_ (params.shear),
        dt_ (params.dt),
        dx_ (xCoords[1] - xCoords[0]),
        dy_ (yCoords[1] - yCoords[0]),
        nx_ (xCoords.size()),
        ny_ (yCoords.size()),
        yCoord_(yCoords),
        bcType_top_ (bc.bcType_top),
        bcType_left_ (bc.bcType_left),
        bcType_right_ (bc.bcType_right),
        bcType_bot_ (bc.bcType_bot),
        bcVals_top_ (bc.bcVals_top),
        bcVals_left_ (bc.bcVals_left),
        bcVals_right_ (bc.bcVals_right),
        bcVals_bot_ (bc.bcVals_bot),
        phi_(phi_init)
    {
        invdx_ = 1.0/dx_;
        invdy_ = 1.0/dy_;
        nInteriorPoints_ = nx_ * ny_;
        nGhostPoints_ = 2*nx_ + 2*ny_;
        nTotalPoints_ = nInteriorPoints_ + nGhostPoints_;

        u_vec_.resize(nInteriorPoints_);
        v_vec_.resize(nInteriorPoints_);
        Dh_vec_.resize(nInteriorPoints_);
        Dv_vec_.resize(nInteriorPoints_);
        rhs_.resize(nTotalPoints_);
        phi_.resize(nTotalPoints_);
        points_.reserve(nTotalPoints_);
        deferredCorr_.resize(nInteriorPoints_);
        deferredCorr_.setZero();
        source_.resize(nInteriorPoints_);
        source_.setZero();
        std::generate_n(std::back_inserter(points_), nTotalPoints_, [] { return std::make_unique<Point>(); });
        totalCoefMatrix_.resize(nTotalPoints_, nTotalPoints_);

        updateDiffusion(params.Dh, params.Dv);
        initVelocVecs();
        buildPointList();
        applyBoundaryCondition();
    }

    void AdvDiffSystem::initVelocVecs(){
        for(int i = 0; i < nx_; i++){
            for(int j = 0; j < ny_; j++){
                int vector_idx = twoDIdx_to_vecIdx(i, j, nx_, ny_, format_);
                u_vec_[vector_idx] = u_double_ - yCoord_[j] * shear_;
                v_vec_[vector_idx] = v_double_;
            }
        }
    }

    void AdvDiffSystem::buildPointList(){
        //interior points
        for(int i = 0; i < nx_; i++){
            for(int j = 0; j < ny_; j++){
                int vector_idx = twoDIdx_to_vecIdx(i, j, nx_, ny_, format_);
                //top interior bound / ghost points
                if(j == ny_ - 1){
                    BoundaryCondDescription bc_ghost = BoundaryCondDescription(bcType_top_, FaceDirection::NORTH, bcVals_top_[i], vector_idx);
                    int corrGhostPoint = nInteriorPoints_ + i;
                    points_[corrGhostPoint] = std::make_unique<GhostPoint>(bc_ghost);

                    BoundaryCondDescription bc_top = BoundaryCondDescription(bcType_top_, FaceDirection::NORTH, bcVals_top_[i], corrGhostPoint);

                    //check for top corner point edge cases
                    //top left
                    if(i == 0){
                        int corrGhostPoint2 = nInteriorPoints_ + nx_ + j;
                        BoundaryCondDescription bc_left = BoundaryCondDescription(bcType_left_, FaceDirection::WEST, bcVals_left_[j], corrGhostPoint2);
                        BoundaryCondDescription bc_ghost2 = BoundaryCondDescription(bcType_left_, FaceDirection::WEST, bcVals_left_[j], vector_idx);
                        points_[vector_idx] = std::make_unique<IntBoundPoint>(bc_top, bc_left);
                        points_[corrGhostPoint2] = std::make_unique<GhostPoint>(bc_ghost2);
                    }
                    //top right
                    else if(i == nx_ - 1){
                        int corrGhostPoint2 = nInteriorPoints_ + nx_ + ny_ + j;
                        BoundaryCondDescription bc_right = BoundaryCondDescription(bcType_right_, FaceDirection::EAST, bcVals_right_[j], corrGhostPoint2);
                        BoundaryCondDescription bc_ghost2 = BoundaryCondDescription(bcType_right_, FaceDirection::EAST, bcVals_right_[j], vector_idx);
                        points_[vector_idx] = std::make_unique<IntBoundPoint>(bc_top, bc_right);
                        points_[corrGhostPoint2] = std::make_unique<GhostPoint>(bc_ghost2);
                    }
                    //if not corner boundary make normal interior node.
                    else{
                        points_[vector_idx] = std::make_unique<IntBoundPoint>(bc_top);
                    }
                }
                //bottom 
                else if (j == 0){
                    int corrGhostPoint = nInteriorPoints_ + nx_ + 2*ny_ + i;
                    BoundaryCondDescription bc_ghost = BoundaryCondDescription(bcType_bot_, FaceDirection::SOUTH, bcVals_bot_[i], vector_idx);
                    points_[corrGhostPoint] = std::make_unique<GhostPoint>(bc_ghost);

                    BoundaryCondDescription bc_bot = BoundaryCondDescription(bcType_bot_, FaceDirection::SOUTH, bcVals_bot_[i], corrGhostPoint);

                    //check for bottom edge cases
                    //bot left
                    if(i == 0){
                        int corrGhostPoint2 = nInteriorPoints_ + nx_ + j;

                        BoundaryCondDescription bc_left = BoundaryCondDescription(bcType_left_, FaceDirection::WEST, bcVals_left_[j], corrGhostPoint2);
                        BoundaryCondDescription bc_ghost2 = BoundaryCondDescription(bcType_left_, FaceDirection::WEST, bcVals_left_[j], vector_idx);
                        
                        points_[vector_idx] = std::make_unique<IntBoundPoint>(bc_bot, bc_left);
                        points_[corrGhostPoint2] = std::make_unique<GhostPoint>(bc_ghost2);
                    }
                    //bot right
                    else if(i == nx_ - 1){
                        int corrGhostPoint2 = nInteriorPoints_ + nx_ + ny_ + j;
                        BoundaryCondDescription bc_right = BoundaryCondDescription(bcType_right_, FaceDirection::EAST, bcVals_right_[j], corrGhostPoint2);
                        BoundaryCondDescription bc_ghost2 = BoundaryCondDescription(bcType_right_, FaceDirection::EAST, bcVals_right_[j], vector_idx);
                        
                        points_[vector_idx] = std::make_unique<IntBoundPoint>(bc_bot, bc_right);
                        points_[corrGhostPoint2] = std::make_unique<GhostPoint>(bc_ghost2);

                    }
                    //if not corner boundary make normal interior node
                    else {
                        points_[vector_idx] = std::make_unique<IntBoundPoint>(bc_bot);
                    }

                }
                // left
                else if (i == 0){
                    int corrGhostPoint = nInteriorPoints_ + nx_ + j;
                    BoundaryCondDescription bc_left = BoundaryCondDescription(bcType_left_, FaceDirection::WEST, bcVals_left_[j], corrGhostPoint);
                    BoundaryCondDescription bc_ghost = BoundaryCondDescription(bcType_left_, FaceDirection::WEST, bcVals_left_[j], vector_idx);

                    points_[vector_idx] = std::make_unique<IntBoundPoint>(bc_left);
                    points_[corrGhostPoint] = std::make_unique<GhostPoint>(bc_ghost);
                } 
                //right 
                else if (i == nx_ - 1){
                    int corrGhostPoint = nInteriorPoints_ + nx_ + ny_ + j;
                    BoundaryCondDescription bc_right = BoundaryCondDescription(bcType_right_, FaceDirection::EAST, bcVals_right_[j], corrGhostPoint);
                    BoundaryCondDescription bc_ghost = BoundaryCondDescription(bcType_right_, FaceDirection::EAST, bcVals_right_[j], vector_idx);

                    points_[vector_idx] = std::make_unique<IntBoundPoint>(bc_right);
                    points_[corrGhostPoint] = std::make_unique<GhostPoint>(bc_ghost);
                }
                else {
                    points_[vector_idx] = std::make_unique<Point>(BoundaryConditionFlag::INTERIOR);
                }
            }
        }
    }
    void AdvDiffSystem::buildCoeffMatrix(bool operatorSplit){
        //Crank-Nicholson Discretization. Builds the Advection terms of the A matrix 
        //in the system A * phi_t+1 = b.

        //Num non-zeros calculation: 
        std::vector<Eigen::Triplet<double>> tripletList;
        auto start = std::chrono::high_resolution_clock::now();
        for(int i = 0; i < nTotalPoints_; i++){

            if(points_[i]->isGhost()){
                switch(points_[i]->bcType()){
                    case BoundaryConditionFlag::DIRICHLET_GHOSTPOINT:{
                        // (phi_int + phi_ghost) / 2 = phi_boundary
                        // if inhomog, the bc value will appear in the rhs.
                        tripletList.emplace_back(i, i, 0.5);
                        tripletList.emplace_back(i, points_[i]->corrPoint(), 0.5);
                        break;
                    }
                    case BoundaryConditionFlag::PERIODIC_GHOSTPOINT:{
                        // Ghost node will exactly equal to the cell on the other end of the domain
                        // effectively makes opposite ends of the domain neighbor points.
                        throw std::runtime_error("Periodic BCs not yet implemented.");
                        break;
                    }
                    default: {
                        throw std::runtime_error("Ghost point doesn't have a bcType associated with being a ghost point!");
                    }
                }
                continue;
            }

            int idx_E = neighbor_point(FaceDirection::EAST, i);
            int idx_W = neighbor_point(FaceDirection::WEST, i);
            int idx_N = neighbor_point(FaceDirection::NORTH, i);
            int idx_S = neighbor_point(FaceDirection::SOUTH, i);

            //Diffusion Terms
            double coeff_C = 1 + 2 * dt_ * (Dh_vec_[i] / (dx_ * dx_) + Dv_vec_[i] / (dy_ * dy_));
            double coeff_E = -Dh_vec_[i] * dt_ / (dx_ * dx_);
            double coeff_W = -Dh_vec_[i] * dt_ / (dx_ * dx_) ;
            double coeff_N = -Dv_vec_[i] * dt_ / (dy_ * dy_);
            double coeff_S = -Dv_vec_[i] * dt_ / (dy_ * dy_);

            //Operator splitting uses implicit only for diffusion
            if(!operatorSplit && (u_double_ > 0 || v_double_ > 0 || shear_ > 0)){
                buildAdvectionCoeffs(i, coeff_C, coeff_N, coeff_S, coeff_E, coeff_W);
            }

            //Triplet Format: row, col, value
            tripletList.emplace_back(i, idx_E, coeff_E);
            tripletList.emplace_back(i, idx_W, coeff_W);
            tripletList.emplace_back(i, idx_N, coeff_N);
            tripletList.emplace_back(i, idx_S, coeff_S);
            tripletList.emplace_back(i, i, coeff_C);
        } 

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
        totalCoefMatrix_.setFromTriplets(tripletList.begin(), tripletList.end());
    }
    void AdvDiffSystem::buildAdvectionCoeffs(int i, double& coeff_C, double& coeff_N, double& coeff_S, double& coeff_E, double& coeff_W){
        //Advection Terms
        //When a boundary condition is in place, phi at the face can be directly calculated using the BC.
        //Therefore, that term goes to the RHS and the contribution of that face to the coeffs goes to 0.

        bool isNorthBoundary = 0, isWestBoundary = 0, isEastBoundary = 0, isSouthBoundary = 0;
        if(points_[i]->bcType() != BoundaryConditionFlag::INTERIOR){
            isNorthBoundary = points_[i]->bcDirection() == FaceDirection::NORTH;
            isSouthBoundary = points_[i]->bcDirection() == FaceDirection::SOUTH;

            //Corner cases...
            bool secondaryWestBound = (points_[i]->secondBoundaryConds() && points_[i]->secondBoundaryConds().value().direction == FaceDirection::WEST);
            bool secondaryEastBound = (points_[i]->secondBoundaryConds() && points_[i]->secondBoundaryConds().value().direction == FaceDirection::EAST);

            isWestBoundary = (points_[i]->bcDirection() == FaceDirection::WEST || secondaryWestBound);
            isEastBoundary = (points_[i]->bcDirection() == FaceDirection::EAST || secondaryEastBound);
        }

        int idx_E = neighbor_point(FaceDirection::EAST, i);
        int idx_W = neighbor_point(FaceDirection::WEST, i);
        int idx_N = neighbor_point(FaceDirection::NORTH, i);
        int idx_S = neighbor_point(FaceDirection::SOUTH, i);

        double u_W = isWestBoundary? u_vec_[i] : 0.5 * (u_vec_[i] + u_vec_[idx_W]);
        double u_E = isEastBoundary? u_vec_[i] : 0.5 * (u_vec_[i] + u_vec_[idx_E]);
        double v_N = isNorthBoundary? v_vec_[i] : 0.5 * (v_vec_[i] + v_vec_[idx_N]);
        double v_S = isSouthBoundary? v_vec_[i] : 0.5 * (v_vec_[i] + v_vec_[idx_S]);

        if (u_E >= 0) coeff_C += u_E * dt_ / dx_;
        if (u_W < 0) coeff_C -= u_W * dt_ / dx_;
        if (v_N >= 0) coeff_C += v_N * dt_ / dy_;
        if (v_S < 0) coeff_C -= v_S * dt_ / dy_;

        if (u_E < 0) coeff_E += u_E * dt_ / dy_;

        if (u_W >= 0) coeff_W -= u_W * dt_ / dx_;

        if (v_N < 0) coeff_N += v_N * dt_ / dy_;

        if (v_S >= 0) coeff_S -= v_S * dt_ / dy_;

        double phi_N_corr = 0, phi_S_corr = 0, phi_W_corr = 0, phi_E_corr = 0;

        auto bool_to_signed = [](bool binary) { return binary ? 1 : -1; };

        if (!isNorthBoundary)
            phi_N_corr -= bool_to_signed(v_N >= 0) * v_N * dt_ / dy_ * 0.5 * minmod(i, FaceDirection::NORTH, 0) * (phi_[idx_N] - phi_[i]);

        if (!isSouthBoundary)
            phi_S_corr += bool_to_signed(v_S >= 0) * v_S * dt_ / dy_ * 0.5 * minmod(i, FaceDirection::SOUTH, 0) * (phi_[i] - phi_[idx_S]);
        
        if (!isWestBoundary)
            phi_W_corr += bool_to_signed(u_W >= 0) * u_W * dt_ / dx_ * 0.5 * minmod(i, FaceDirection::WEST, 0) * (phi_[i] - phi_[idx_W]);

        if (!isEastBoundary)
            phi_E_corr -= bool_to_signed(u_E >= 0) * u_E * dt_ / dx_ * 0.5 * minmod(i, FaceDirection::EAST, 0) * (phi_[idx_E] - phi_[i]);
        
        double TVD_deferred_corr = phi_N_corr + phi_S_corr + phi_W_corr + phi_E_corr;
        deferredCorr_[i] = TVD_deferred_corr;
    }

    const Eigen::VectorXd& AdvDiffSystem::calcRHS(){
        for(int i = 0; i < nTotalPoints_; i++){
            if(points_[i]->isGhost()){
                switch(points_[i]->bcType()){
                    case BoundaryConditionFlag::DIRICHLET_GHOSTPOINT:{
                        // Equation: (phi_int + phi_ghost) / 2 = phi_boundary
                        rhs_[i] = points_[i]->bcVal();
                        break;
                    }
                    case BoundaryConditionFlag::PERIODIC_GHOSTPOINT:{
                        // Ghost node will exactly equal to the cell on the other end of the domain
                        // effectively makes opposite ends of the domain neighbor points.
                        throw std::runtime_error("Periodic BCs not yet implemented.");
                        break;
                    }
                    default: {
                        throw std::runtime_error("Ghost point doesn't have a bcType associated with being a ghost point!");
                    }                
                }
                continue;
            }

            switch(points_[i]->bcType()){
                case BoundaryConditionFlag::INTERIOR:{
                    rhs_[i] = phi_[i] + deferredCorr_[i] + source_[i]*dt_;
                    break;
                }
                case BoundaryConditionFlag::DIRICHLET_INT_BPOINT:{
                    rhs_[i] = phi_[i] + deferredCorr_[i] + source_[i]*dt_;
                    switch(points_[i]->bcDirection()){
                        case FaceDirection::NORTH:
                            rhs_[i] -= v_vec_[i] * dt_ / dy_ * points_[i]->bcVal();
                            break;
                        case FaceDirection::SOUTH:
                            rhs_[i] += v_vec_[i] * dt_ / dy_ * points_[i]->bcVal();
                            break;
                        case FaceDirection::EAST:
                            rhs_[i] -= u_vec_[i] * dt_ / dx_ * points_[i]->bcVal();
                            break;
                        case FaceDirection::WEST:
                            rhs_[i] += u_vec_[i] * dt_ / dx_ * points_[i]->bcVal();
                            break;
                        case FaceDirection::ERROR:
                            throw std::runtime_error("Invalid FaceDirection in Dirichlet boundary condition");
                    }
                    if (!points_[i]->secondBoundaryConds()) break;
                    BoundaryCondDescription bc_2 = points_[i]->secondBoundaryConds().value();
                    switch(bc_2.direction){
                        case FaceDirection::EAST:
                            rhs_[i] -= u_vec_[i] * dt_ / dx_ * bc_2.bcVal;
                            break;
                        case FaceDirection::WEST:
                            rhs_[i] += u_vec_[i] * dt_ / dx_ * bc_2.bcVal;
                            break;
                        default:
                            throw std::runtime_error("Can't have anything but EAST or WEST as secondary BC!");
                    }
                    break;
                }
                case BoundaryConditionFlag::PERIODIC_INT_BPOINT:{
                    throw std::runtime_error("Periodic BCs not yet implemented.");
                    break;
                }
                default: {
                    throw std::runtime_error("Interior boundary point has invalid bcType");
                }                
            }
            continue;
        }

        return rhs_;
    }
    void AdvDiffSystem::applyBoundaryCondition(){
        //top and bottom bc
        for(int i = 0; i < nx_; i++){
            //top
            int bPointID_top = twoDIdx_to_vecIdx(i, ny_ - 1, nx_, ny_, format_);
            switch(bcType_top_){
                case BoundaryConditionFlag::DIRICHLET_INT_BPOINT: {
                    int ghostPointID = points_[bPointID_top]->corrPoint();
                    phi_[ghostPointID] = 2 * points_[bPointID_top]->bcVal() - phi_[bPointID_top];
                    break;
                }
                default: {
                    throw std::runtime_error("Chosen boundary condition not implemented yet");
                }
            }

            //bottom
            int bPointID_bot = twoDIdx_to_vecIdx(i, 0, nx_, ny_, format_);
            switch(bcType_bot_){
                case BoundaryConditionFlag::DIRICHLET_INT_BPOINT: {
                    int ghostPointID = points_[bPointID_bot]->corrPoint();
                    phi_[ghostPointID] = 2 * points_[bPointID_bot]->bcVal() - phi_[bPointID_bot];
                    break;
                }
                default: {
                    throw std::runtime_error("Chosen boundary condition not implemented yet");
                }
            }
        }

        //left and right bc
        for(int j = 0; j < ny_; j++){
            //corner cases
            if(j == 0 || j == ny_ - 1){
                int bPointID_cornerLeft = twoDIdx_to_vecIdx(0, j, nx_, ny_, format_);
                int ghostPointID = points_[bPointID_cornerLeft]->secondBoundaryConds().value().corrPoint;
                double bcVal =  points_[bPointID_cornerLeft]->secondBoundaryConds().value().bcVal;
                switch(bcType_left_){
                    case BoundaryConditionFlag::DIRICHLET_INT_BPOINT: {
                        phi_[ghostPointID] = 2 * bcVal - phi_[bPointID_cornerLeft];
                        break;
                    }
                    default: {
                        throw std::runtime_error("Chosen boundary condition not implemented yet");
                    }
                }

                int bPointID_cornerRight = twoDIdx_to_vecIdx(nx_ - 1, j, nx_, ny_, format_);
                ghostPointID = points_[bPointID_cornerRight]->secondBoundaryConds().value().corrPoint;
                bcVal = points_[bPointID_cornerRight]->secondBoundaryConds().value().bcVal;

                switch(bcType_right_){
                    case BoundaryConditionFlag::DIRICHLET_INT_BPOINT: {
                        phi_[ghostPointID] = 2 * bcVal - phi_[bPointID_cornerRight];
                        break;
                    }
                    default: {
                        throw std::runtime_error("Chosen boundary condition not implemented yet");
                    }
                }
            }

            //left
            int bPointID_left = twoDIdx_to_vecIdx(0, j, nx_, ny_, format_);
            switch(bcType_left_){
                case BoundaryConditionFlag::DIRICHLET_INT_BPOINT: {
                    int ghostPointID = points_[bPointID_left]->corrPoint();
                    phi_[ghostPointID] = 2 * points_[bPointID_left]->bcVal() - phi_[bPointID_left];
                    break;
                }
                default: {
                    throw std::runtime_error("Chosen boundary condition not implemented yet");
                }
            }
            //right
            int bPointID_right = twoDIdx_to_vecIdx(nx_ - 1, j, nx_, ny_, format_);
            switch(bcType_right_){
                case BoundaryConditionFlag::DIRICHLET_INT_BPOINT: {
                    int ghostPointID = points_[bPointID_right]->corrPoint();
                    phi_[ghostPointID] = 2 * points_[bPointID_right]->bcVal() - phi_[bPointID_right];
                    break;
                }
                default: {
                    throw std::runtime_error("Chosen boundary condition not implemented yet");
                }
            }
        }
    }
    void AdvDiffSystem::updateBoundaryCondition(const BoundaryConditions& bc){
        bcType_top_ = bc.bcType_top;
        bcType_left_ = bc.bcType_left;
        bcType_right_ = bc.bcType_right;
        bcType_bot_ = bc.bcType_bot;
        bcVals_top_ = bc.bcVals_top;
        bcVals_left_ = bc.bcVals_left;
        bcVals_right_ = bc.bcVals_right;
        bcVals_bot_ = bc.bcVals_bot;

        //Go through ghost points, and update the bcType and value of them and their corresponding interior nodes
        //As seen in buildPointList(), ghost point order goes top->left->right->bottom
        int currIdx = nInteriorPoints_;
        //top
        for(int i = 0; i < nx_; i++){
            int corrPointID = points_[currIdx]->corrPoint();
            points_[currIdx]->setBCType(bcType_top_); 
            points_[currIdx]->setBCVal(bcVals_top_[i]); 
            points_[corrPointID]->setBCType(bcType_top_); 
            points_[corrPointID]->setBCVal(bcVals_top_[i]); 
            currIdx++;
        }
        //left
        for(int i = 0; i < ny_; i++){
            int corrPointID = points_[currIdx]->corrPoint();
            if(i == 0 || i == ny_ - 1){
                points_[currIdx]->setBCType(bcType_left_); 
                points_[currIdx]->setBCVal(bcVals_left_[i]);
                BoundaryCondDescription bc(bcType_left_, FaceDirection::WEST, bcVals_left_[i], currIdx);
                points_[corrPointID]->setSecondaryBC(bc);
                currIdx++;
                continue;
            }
            points_[currIdx]->setBCType(bcType_left_); 
            points_[currIdx]->setBCVal(bcVals_left_[i]); 
            points_[corrPointID]->setBCType(bcType_left_); 
            points_[corrPointID]->setBCVal(bcVals_left_[i]);
            currIdx++;
        }
        //right
        for(int i = 0; i < ny_; i++){
            int corrPointID = points_[currIdx]->corrPoint();
            if(i == 0 || i == ny_ - 1){
                points_[currIdx]->setBCType(bcType_right_); 
                points_[currIdx]->setBCVal(bcVals_right_[i]);
                BoundaryCondDescription bc(bcType_right_, FaceDirection::EAST, bcVals_right_[i], currIdx);
                points_[corrPointID]->setSecondaryBC(bc);
                currIdx++;
                continue;
            }
            points_[currIdx]->setBCType(bcType_right_); 
            points_[currIdx]->setBCVal(bcVals_right_[i]); 
            points_[corrPointID]->setBCType(bcType_right_); 
            points_[corrPointID]->setBCVal(bcVals_right_[i]);
            currIdx++;
        }
        //bot
        for(int i = 0; i < nx_; i++){
            int corrPointID = points_[currIdx]->corrPoint();
            points_[currIdx]->setBCType(bcType_bot_); 
            points_[currIdx]->setBCVal(bcVals_bot_[i]); 
            points_[corrPointID]->setBCType(bcType_bot_); 
            points_[corrPointID]->setBCVal(bcVals_bot_[i]); 
            currIdx++;
        }
        applyBoundaryCondition(); //need this to calculate minmod function at some timestep.
    }

    Eigen::VectorXd AdvDiffSystem::forwardEulerAdvection(bool operatorSplit, bool parallelAdvection) const noexcept{
        Eigen::VectorXd soln(nTotalPoints_);
        // double avgBackgroundCalcTime = 0;
        //Explicit Time-Stepping
        #pragma omp parallel for    \
        if      ( parallelAdvection ) \
        default ( shared          ) \
        schedule( static, 100      )
        for(int i = 0; i < nInteriorPoints_; i++){
            //When a boundary condition is in place, phi at the face can be directly calculated using the BC.
            //Therefore, that term goes to the RHS and the contribution of that face to the coeffs goes to 0.
            bool isNorthBoundary = 0, isWestBoundary = 0, isEastBoundary = 0, isSouthBoundary = 0, secondaryWestBound = 0, secondaryEastBound = 0;
            int idx_E = i + ny_;
            int idx_W = i - ny_;
            int idx_N = i + 1;
            int idx_S = i - 1;

            //commenting out this results in ~30% speedup
            //The calls involving the optional are maybe 1/3 of the cost. Maybe something to look at later.
            if(points_[i]->bcType() != BoundaryConditionFlag::INTERIOR){
                Point* point = points_[i].get();
                FaceDirection direction = point->bcDirection();
                isNorthBoundary = direction == FaceDirection::NORTH;
                isSouthBoundary = direction == FaceDirection::SOUTH;

                //Corner cases...
                bool secondaryWestBound = (point->secondBoundaryConds() && point->secondBoundaryConds()->direction == FaceDirection::WEST);
                bool secondaryEastBound = (point->secondBoundaryConds() && point->secondBoundaryConds()->direction == FaceDirection::EAST);

                isWestBoundary = (direction == FaceDirection::WEST || secondaryWestBound);
                isEastBoundary = (direction == FaceDirection::EAST || secondaryEastBound);

                //only call this lookup function on boundary nodes which are inconsequential in number
                idx_N = isNorthBoundary? point->corrPoint() : idx_N;
                idx_S = isSouthBoundary? point->corrPoint() : idx_S;
                idx_E = isEastBoundary? (secondaryEastBound ? point->secondBoundaryConds()->corrPoint : point->corrPoint()) : idx_E;
                idx_W = isWestBoundary? (secondaryEastBound ? point->secondBoundaryConds()->corrPoint : point->corrPoint()) : idx_W;
            }
            //When you declare these vars (inside or outside loop) has 0 impact)
            //takes ~ 6 out of 18 ns on background var calcs

            //these cost almost nothing to compute but commenting out anyway for maximum performance
            // double dphi_dx_E = (phi_[idx_E] - phi_[i]) * invdx_;
            // double dphi_dx_W = (phi_[i] - phi_[idx_W]) * invdx_;
            // double dphi_dy_N = (phi_[idx_N] - phi_[i]) * invdy_;
            // double dphi_dy_S = (phi_[i] - phi_[idx_S]) * invdy_;
            
            //ignoreing distinction of faces saves a good amt of time
            // double u_W = isWestBoundary? u_vec_[i] : 0.5 * (u_vec_[i] + u_vec_[idx_W]);
            // double u_E = isEastBoundary? u_vec_[i] : 0.5 * (u_vec_[i] + u_vec_[idx_E]);
            // double v_N = isNorthBoundary? v_vec_[i] : 0.5 * (v_vec_[i] + v_vec_[idx_N]);
            // double v_S = isSouthBoundary? v_vec_[i] : 0.5 * (v_vec_[i] + v_vec_[idx_S]);
            double u_local = u_vec_[i];
            double v_local = v_vec_[i];

            // auto stop = std::chrono::high_resolution_clock::now();
            // auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            // avgBackgroundCalcTime += duration.count();
            //std::cout << "ForwardEuler: Background Variable Calc Time: " << duration.count() << "ns" << std::endl;
            // start = std::chrono::high_resolution_clock::now();
            double phi_N, phi_S, phi_W, phi_E;

            //Unraveling any of these if's into single liners hurts performance
            //Killing the branching completely into 1 statement (not possible) only results in ~10% speedup (not worth it)
            //Using only first order upwind can result in a ~40% speedup of the total advection calc.
            //So... there is significantly more cost from actually doing the calculation than from branching.
            if(isNorthBoundary){
                phi_N = points_[i]->bcVal();
            }
            else if (v_local >= 0){
                phi_N = phi_[i] + 0.5 * minmod_N_vPos(i) * (phi_[idx_N] - phi_[i]);
            }
            else {
                phi_N = phi_[idx_N] + 0.5 * minmod_N_vNeg(i) * (phi_[i] - phi_[idx_N]);
            }
            if(isSouthBoundary){
                phi_S = points_[i]->bcVal();
            }
            else if (v_local >= 0){
                phi_S = phi_[idx_S] +  0.5 * minmod_S_vPos(i) * (phi_[i] - phi_[idx_S]);
            }
            else {
                phi_S = phi_[i] +  0.5 * minmod_S_vNeg(i) * (phi_[idx_S] - phi_[i]);
            }

            if(isWestBoundary){
                phi_W = secondaryWestBound ? points_[i]->secondBoundaryConds().value().bcVal : points_[i]->bcVal();
            }
            else if (u_local >= 0){
                phi_W = phi_[idx_W] + 0.5 * minmod_W_vPos(i) * (phi_[i] - phi_[idx_W]);
            }
            else {
                phi_W = phi_[i] + 0.5 * minmod_W_vNeg(i) * (phi_[idx_W] - phi_[i]);
            }

            if(isEastBoundary){
                phi_E = secondaryEastBound ? points_[i]->secondBoundaryConds().value().bcVal : points_[i]->bcVal();
            }
            else if (u_local >= 0){
                phi_E = phi_[i] + 0.5 * minmod_E_vPos(i) * (phi_[idx_E] - phi_[i]);
            }
            else {
                phi_E = phi_[idx_E] + 0.5 * minmod_E_vNeg(i) * (phi_[i] - phi_[idx_E]);
            }

            //std::cout << "ForwardEuler: Fluxes and Update Time: " << duration.count() << "ns" << std::endl;

            //Even just setting this to 0 is like a 2 ns save out of 12, not sure if worth
            soln[i] = /*(!operatorSplit) * (Dh_ * dt_ * invdx_ * (dphi_dx_E - dphi_dx_W) + Dv_ * dt_ * invdy_ * (dphi_dy_N - dphi_dy_S))\*/
                     dt_ * invdx_ * (u_local * phi_W - u_local * phi_E) + dt_ * invdy_ * (v_local * phi_S - v_local * phi_N)\
                    + source_[i] * dt_ + phi_[i];
            // stop = std::chrono::high_resolution_clock::now();
            // duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            // avgFluxCalcTime += duration.count();
        }
        return soln;
    }
    
void sor_solve(const Eigen::SparseMatrix<double, Eigen::RowMajor> &A, const Eigen::VectorXd &rhs, Eigen::VectorXd &phi, double omega, double threshold, int n_iters) {
    /*
    diagCoeff should always be overwritten in the for loop
    before we get to "x_i *= omega / diagCoeff;"
    Setting it to 0 instead of leaving uninitialized guarantees
    that we get an error if for some reason diagCoeff is not overwritten
    Not sure how to do this better for now...
    */ 
    double diagCoeff = 0;
    double residual = 1;

    while(residual > threshold){
        const double* valuePtr = A.valuePtr();
        const int* innerIdxPtr = A.innerIndexPtr();
        const int* outerIdxPtr = A.outerIndexPtr();

        for(int iteration = 0; iteration < n_iters; iteration++){
            int outerIdx = 0;

            for (int i = 0; i < rhs.size(); i++) {
                double x_i = 0;
                int rowStartIdx = outerIdxPtr[outerIdx];
                int rowEndIdx = outerIdxPtr[outerIdx + 1];
                for (int j = rowStartIdx; j < rowEndIdx; j++) {

                    if (innerIdxPtr[j] == i) {
                        diagCoeff = valuePtr[j];
                        continue;
                    }
                    x_i -= valuePtr[j] * phi[innerIdxPtr[j]];
                }
                x_i += rhs[i];

                x_i *= omega / diagCoeff;
                x_i += (1 - omega) * phi[i];
                phi[i] = x_i;
                outerIdx++;
            } // end inner for loop
        } // end iters for loop
    
        residual = (A * phi - rhs).eval().lpNorm<2>()/ rhs.lpNorm<2>();
        if (isnan(residual)) throw std::runtime_error("NaN residual encountered");
    } // end while loop

}

}