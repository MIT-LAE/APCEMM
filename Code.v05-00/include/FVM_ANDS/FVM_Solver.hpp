#include "FVM_ANDS/AdvDiffSystem.hpp"
#include <unsupported/Eigen/IterativeSolvers>
#include <math.h>
#ifndef FVM_ANDS_SOLVER_H
#define FVM_ANDS_SOLVER_H
namespace FVM_ANDS{
    class FVM_Solver{
        public:
            FVM_Solver(const AdvDiffParams& params, const Vector_1D xCoords, const Vector_1D yCoords, const BoundaryConditions& bc, const Eigen::VectorXd& phi_init, bool useDiagPreCond = false, int maxIters_ = 1000, double convergenceThres_ = 1e-5);
            const Eigen::VectorXd& solve();
            const Eigen::VectorXd& solve(const Eigen::VectorXd& source);
            const Eigen::VectorXd& explicitSolve();
            const Eigen::VectorXd& explicitSolve(const Eigen::VectorXd& source);

            const Eigen::VectorXd& operatorSplitSolve(bool parallelAdvection = false, double courant_max = 0.5);
            void operatorSplitSolve2DVec(Vector_2D& vec, const BoundaryConditions& bc, bool parallelAdvection = false, double courant_max = 0.5);

            void advectionHalfTimestepSolve(Vector_2D& vec, const BoundaryConditions& bc, double courant_max = 0.5);

            void buildCoeffMatrix(bool operatorSplit = false){
                advDiffSys_.buildCoeffMatrix(operatorSplit);
            }
            inline void updateTimestep(double dt){
                advDiffSys_.updateTimestep(dt);
            }
            inline void updateDiffusion(double Dh, double Dv){
                advDiffSys_.updateDiffusion(Dh, Dv);
            }
            inline void updateDiffusion(const Vector_2D& Dh, const Vector_2D& Dv){
                advDiffSys_.updateDiffusion(Dh, Dv);
            }
            inline void updateAdvection(double u, double v, double shear){
                advDiffSys_.updateAdvection(u, v, shear);
            }
            inline void setConvergenceThres(double tol){
                convergenceThres_ = tol;
            }
            inline void setMaxIters(int iters){
                maxIters_ = iters;
            }
            inline void updateBoundaryCondition(const BoundaryConditions& bc){
                advDiffSys_.updateBoundaryCondition(bc);
            }
            inline void updatePhi(const Eigen::VectorXd& phi){
                advDiffSys_.updatePhi(phi);
            }
            inline void updateSpacing(const Vector_1D& yCoords_new, double dx_new, int nx_new) {
                advDiffSys_.updateSpacing(yCoords_new, dx_new, nx_new);
            }
            inline const Eigen::VectorXd& phi(){
                return advDiffSys_.phi();
            }
            inline const Eigen::SparseMatrix<double, Eigen::RowMajor>& coefMatrix(){
                return advDiffSys_.getCoefMatrix();
            }
            inline const std::vector<std::unique_ptr<Point>>& points(){
                return advDiffSys_.points();
            }
            const Eigen::VectorXd& calcRHS(){
                return advDiffSys_.calcRHS();
            }
            inline int numIters(){
                return solver_.iterations();
            }
            inline double error(){
                return solver_.error();
            }
            //For some weird reason, Eigen uses FLOAT accuracy to calcuate norms.
            //Therefore with the extremely small numbers we have in some bins at the start,
            //squaredNorm() or lpNorm<2>() will return zero.
            static double eigenSqVectorNorm_double(const Eigen::VectorXd& vec) {
                double sum = 0;
                for(int i = 0; i < vec.rows(); i++){
                        sum += vec[i] * vec[i];
                }
                sum = sqrt(sum);
                return sum;
            }
            static double eigenSqVectorNorm_double(const Vector_2D& vec) {
                double sum = 0;
                for (std::size_t j = 0; j < vec.size(); j++) {
                    for (std::size_t i = 0; i < vec[0].size(); i++) {
                        if(isnan(vec[j][i])){
                            //std::cout << j << " " << i << std::endl;
                        }
                        sum += vec[j][i] * vec[j][i];
                    }
                }
                sum = sqrt(sum);
                return sum;
            }
        private:
            int maxIters_;
            double convergenceThres_;
            AdvDiffSystem advDiffSys_;
            bool useDiagPreCond_;
            Eigen::DiagonalMatrix<double, -1> diagPreCond;
            Eigen::DiagonalMatrix<double, -1> diagPreCond_inv;
            Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::DiagonalPreconditioner<double> > solver_;
    };
} 
#endif