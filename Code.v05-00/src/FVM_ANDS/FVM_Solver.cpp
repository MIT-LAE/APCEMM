#include "FVM_ANDS/FVM_Solver.hpp"
namespace FVM_ANDS{
    FVM_Solver::FVM_Solver(const AdvDiffParams& params, const Vector_1D xCoords, const Vector_1D yCoords, const BoundaryConditions& bc, const Eigen::VectorXd& phi_init, bool useDiagPreCond, int maxIters, double convergenceThres)
    :   advDiffSys_(AdvDiffSystem(params, xCoords, yCoords, bc, phi_init)),
        maxIters_(maxIters),
        convergenceThres_(convergenceThres),
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
        //Strang Splitting
        //Step 1: Calculate explicit advection timestep based on CFL condition set

        bool operatorSplit = true;
        double courant = advDiffSys_.courant();
        double dt_max = advDiffSys_.timestep();
        double dt_adv = dt_max * (courant_max / courant);

        int n_timesteps_advection_half =  std::ceil((0.5 * dt_max) / dt_adv);
        dt_adv = (0.5 * dt_max) / n_timesteps_advection_half;

        // auto start = std::chrono::high_resolution_clock::now();

        //Step 2: Solve Advection for half timestep
        advDiffSys_.updateTimestep(dt_adv);
        for(int i = 0; i < n_timesteps_advection_half; i++){
            advDiffSys_.updatePhi(advDiffSys_.forwardEulerAdvection(operatorSplit, parallelAdvection));
            advDiffSys_.applyBoundaryCondition();
        }

        // auto stop = std::chrono::high_resolution_clock::now();
        // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
        // std::cout << "First Halfstep Advection Solve Time: " << duration.count() << std::endl;

        // start = std::chrono::high_resolution_clock::now();

        //Step 3: Implicitly solve diffusion (first to help smoothen out potential steep gradients)
        advDiffSys_.updateTimestep(dt_max);
        //Build matrix takes ~40ms atm
        advDiffSys_.buildCoeffMatrix(operatorSplit);
        advDiffSys_.calcRHS();

        // auto mat = advDiffSys_.getCoefMatrix();
        // auto b = advDiffSys_.getRHS();
        // solver_.compute(mat);
        // Eigen::VectorXd solution = solver_.solveWithGuess(b, advDiffSys_.phi());
        // advDiffSys_.updatePhi(std::move(solution));
        
        advDiffSys_.sor_solve();

        // stop = std::chrono::high_resolution_clock::now();
        // duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
        // std::cout << "Diffusion Solve Time: " << duration.count() << std::endl;

        //Step 4: Explicitly solve advection to full timestep

        // start = std::chrono::high_resolution_clock::now();
        advDiffSys_.updateTimestep(dt_adv);
        for(int i = 0; i < n_timesteps_advection_half; i++){
            advDiffSys_.updatePhi(advDiffSys_.forwardEulerAdvection(operatorSplit));
            advDiffSys_.applyBoundaryCondition();
        }

        advDiffSys_.updateTimestep(dt_max);

        // stop = std::chrono::high_resolution_clock::now();
        // duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
        // std::cout << "Advection second step Solve Time: " << duration.count() << std::endl;

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

    void FVM_Solver::advectionHalfTimestepSolve(Vector_2D& vec, const BoundaryConditions& bc, double courant_max){
        Eigen::VectorXd vec_Eigen = std2dVec_to_eigenVec(vec);
        advDiffSys_.updatePhi(vec_Eigen);
        advDiffSys_.updateBoundaryCondition(bc);

        bool operatorSplit = true;
        double courant = advDiffSys_.courant();
        double dt_max = advDiffSys_.timestep();
        double dt_adv = dt_max * (courant_max / courant);

        int n_timesteps_advection_half =  std::ceil((0.5 * dt_max) / dt_adv);
        dt_adv = (0.5 * dt_max) / n_timesteps_advection_half;

        advDiffSys_.updateTimestep(dt_adv);
        for(int i = 0; i < n_timesteps_advection_half; i++){
            advDiffSys_.updatePhi(advDiffSys_.forwardEulerAdvection(operatorSplit));
            advDiffSys_.applyBoundaryCondition();
        } 
        vec = eigenVec_to_std2dVec(advDiffSys_.phi(), vec[0].size(), vec.size());
    }

}