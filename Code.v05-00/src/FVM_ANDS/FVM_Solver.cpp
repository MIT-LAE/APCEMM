#include <chrono>
#include <iostream>
#include <stdexcept>
#include "APCEMM.h"
#include "FVM_ANDS/FVM_Solver.hpp"
namespace FVM_ANDS{
    FVM_Solver::FVM_Solver(const AdvDiffParams& params, const Vector_1D& xCoords, const Vector_1D& yCoords, const BoundaryConditions& bc, const Eigen::VectorXd& phi_init, bool useDiagPreCond, int maxIters, double convergenceThres)
    :   maxIters_(maxIters),
        convergenceThres_(convergenceThres),
        advDiffSys_(AdvDiffSystem(params, xCoords, yCoords, bc, phi_init)),
        useDiagPreCond_(useDiagPreCond){
        solver_.setTolerance(convergenceThres_);
        solver_.setMaxIterations(maxIters_);

        //For SOR solves only, normally Eigen solvers come pre-packed with diagonal preconditioner.
        if(useDiagPreCond_){
            auto diagonalVec = advDiffSys_.getCoefMatrix().diagonal();
            auto diagonalMat = diagonalVec.asDiagonal();
            diagPreCond = diagonalMat.inverse();
            diagPreCond_inv = diagonalMat;
        }
    }
    
    const Eigen::VectorXd& FVM_Solver::solve(){
        //auto start = std::chrono::high_resolution_clock::now();
        advDiffSys_.buildCoeffMatrix();
        advDiffSys_.calcRHS();
        auto mat = advDiffSys_.getCoefMatrix();
        auto b = advDiffSys_.getRHS();
        solver_.compute(mat);
        Eigen::VectorXd solution = solver_.solveWithGuess(b, advDiffSys_.phi());
        advDiffSys_.updatePhi(std::move(solution));
        return advDiffSys_.phi();
    }

    const Eigen::VectorXd& FVM_Solver::solve(const Eigen::VectorXd& source){
        advDiffSys_.buildCoeffMatrix();
        advDiffSys_.addSource(source);
        advDiffSys_.calcRHS();
        auto mat = advDiffSys_.getCoefMatrix();
        auto b = advDiffSys_.getRHS();
        solver_.compute(mat);
        Eigen::VectorXd solution = solver_.solveWithGuess(b, advDiffSys_.phi());
        advDiffSys_.updatePhi(std::move(solution));
        return advDiffSys_.phi();
    }
    const Eigen::VectorXd& FVM_Solver::explicitSolve(){
        //This doesn't work on diffusion. No good reason to use explicit diff anyway.
        Eigen::VectorXd solution = advDiffSys_.forwardEulerAdvection();
        advDiffSys_.updatePhi(std::move(solution));
        advDiffSys_.applyBoundaryCondition();
        return advDiffSys_.phi();
    }
    const Eigen::VectorXd& FVM_Solver::explicitSolve(const Eigen::VectorXd& source){
        advDiffSys_.addSource(source);
        Eigen::VectorXd solution = advDiffSys_.forwardEulerAdvection();
        advDiffSys_.updatePhi(solution); //No point in using std::move because we need to "slice" the vector and as a result copy stuff anyway
        advDiffSys_.applyBoundaryCondition();
        return advDiffSys_.phi();
    }

    const Eigen::VectorXd& FVM_Solver::operatorSplitSolve(bool parallelAdvection, double courant_max) {

        #ifdef ENABLE_TIMING
        auto opsplit_start = std::chrono::high_resolution_clock::now();
        #endif

        //Strang Splitting
        bool operatorSplit = true;
        double dt_max = advDiffSys_.timestep();

        #ifdef ENABLE_TIMING
        auto start = std::chrono::high_resolution_clock::now();
        #endif

        //Step 1: Solve Advection for half timestep via Semi-Lagrangian advection
        advDiffSys_.semiLagrangianAdvection(0.5 * dt_max, parallelAdvection);

        #ifdef ENABLE_TIMING
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
        std::cout << "              First Halfstep Advection Solve Time: " << duration.count() << " ms" << std::endl;

        start = std::chrono::high_resolution_clock::now();
        #endif

        //Step 2: Implicitly solve diffusion (first to help smoothen out potential steep gradients)
        advDiffSys_.updateTimestep(dt_max);

        // Should never happen given using operatorSplit is hard coded into runTransport and above, but serves as
        // a guardrail for future code changes
        if (!operatorSplit && matrix_prebuilt_){
            throw std::runtime_error("When not using operatorSplitSolve, the matrix coefficients cannot be reused");
        }

        // If the matrix has not already been set, compute it
        if (!matrix_prebuilt_){
            advDiffSys_.buildCoeffMatrix(operatorSplit);
        }

        #ifdef ENABLE_TIMING
        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
        std::cout << "              Build diffusion matrix Time: " << duration.count() << " ms" << std::endl;

        start = std::chrono::high_resolution_clock::now();
        #endif

        advDiffSys_.calcRHS();

        #ifdef ENABLE_TIMING
        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
        std::cout << "              Calc RHS matrix Time: " << duration.count() << " ms" << std::endl;

        start = std::chrono::high_resolution_clock::now();
        #endif

        advDiffSys_.sor_solve(1.0, sorThreshold_);

        #ifdef ENABLE_TIMING
        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
        std::cout << "              Diffusion SOR Solve Time: " << duration.count() << " ms" << std::endl;

        start = std::chrono::high_resolution_clock::now();
        #endif

        //Step 3: Solve advection for second half timestep via Semi-Lagrangian advection
        advDiffSys_.semiLagrangianAdvection(0.5 * dt_max, parallelAdvection);

        #ifdef ENABLE_TIMING
        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
        std::cout << "              Advection second step Solve Time: " << duration.count() << " ms" << std::endl;

        duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-opsplit_start);
        std::cout << "          Ran operatorSplitSolve in " << duration.count() << " ms" << std::endl;
        #endif

        return advDiffSys_.phi();
    }

    void FVM_Solver::operatorSplitSolve2DVec(Vector_2D& vec, const BoundaryConditions& bc, bool parallelAdvection, double courant_max ) { 
        Eigen::VectorXd vec_Eigen = std2dVec_to_eigenVec(vec);
        const double VECTORNORM_MIN = 1e-100;
        //Eigen seems to lose way too much precision in calculations with very small numbers.
        //Not sure if that's fixable. For now just ignore these very small numbers before machine precision becomes relevant.
        // std::cout << eigenSqVectorNorm_double(vec)<< std::endl;
        if(eigenSqVectorNorm_double(vec_Eigen) < VECTORNORM_MIN){
            return;
        }
        advDiffSys_.updatePhi(vec_Eigen);
        advDiffSys_.updateBoundaryCondition(bc);
        vec = eigenVec_to_std2dVec(operatorSplitSolve(parallelAdvection, courant_max), vec[0].size(), vec.size());
    }

    void FVM_Solver::operatorSplitSolveTwoMoment2DVec(
        Vector_2D& vecN,
        Vector_2D& vecV,
        const BoundaryConditions& bc,
        double q_min,
        double q_max,
        bool parallelAdvection)
    {
        Eigen::VectorXd vecN_Eigen = std2dVec_to_eigenVec(vecN);
        const double VECTORNORM_MIN = 1e-100;
        if (eigenSqVectorNorm_double(vecN_Eigen) < VECTORNORM_MIN) {
            return;
        }
        Eigen::VectorXd vecV_Eigen = std2dVec_to_eigenVec(vecV);

        const bool operatorSplit = true;
        const double dt_max = advDiffSys_.timestep();
        advDiffSys_.updateBoundaryCondition(bc);

        // Step 1: First half-step coupled two-moment Semi-Lagrangian advection
        advDiffSys_.semiLagrangianAdvectionTwoMoment(vecN_Eigen, vecV_Eigen, 0.5 * dt_max, q_min, q_max, parallelAdvection);

        // Step 2: Implicit diffusion for full timestep on both moments using identical linear operator
        advDiffSys_.updateTimestep(dt_max);
        if (!matrix_prebuilt_) {
            advDiffSys_.buildCoeffMatrix(operatorSplit);
        }

        advDiffSys_.updatePhi(vecN_Eigen);
        advDiffSys_.applyBoundaryCondition();
        advDiffSys_.calcRHS();
        advDiffSys_.sor_solve(1.0, sorThreshold_);
        vecN_Eigen = advDiffSys_.phi();

        advDiffSys_.updatePhi(vecV_Eigen);
        advDiffSys_.applyBoundaryCondition();
        advDiffSys_.calcRHS();
        advDiffSys_.sor_solve(1.0, sorThreshold_);
        vecV_Eigen = advDiffSys_.phi();

        // Step 3: Second half-step coupled two-moment Semi-Lagrangian advection
        advDiffSys_.semiLagrangianAdvectionTwoMoment(vecN_Eigen, vecV_Eigen, 0.5 * dt_max, q_min, q_max, parallelAdvection);

        vecN = eigenVec_to_std2dVec(vecN_Eigen, vecN[0].size(), vecN.size());
        vecV = eigenVec_to_std2dVec(vecV_Eigen, vecV[0].size(), vecV.size());
    }

    void FVM_Solver::advectionHalfTimestepSolve(Vector_2D& vec, const BoundaryConditions& bc, double courant_max){
        Eigen::VectorXd vec_Eigen = std2dVec_to_eigenVec(vec);
        advDiffSys_.updatePhi(vec_Eigen);
        advDiffSys_.updateBoundaryCondition(bc);

        double dt_max = advDiffSys_.timestep();
        advDiffSys_.semiLagrangianAdvection(0.5 * dt_max);
        vec = eigenVec_to_std2dVec(advDiffSys_.phi(), vec[0].size(), vec.size());
    }

}