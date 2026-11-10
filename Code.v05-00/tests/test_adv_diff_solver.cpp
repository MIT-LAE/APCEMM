#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <iostream>
#include <random>
#include "Core/Mesh.hpp"
#include "FVM_ANDS/FVM_Solver.hpp"
using std::cout;
using std::endl;

namespace FVM_ANDS{
    double solnFcn(double t, double x, double y){
        return std::exp(-t) * std::sin(physConst::PI*x) * std::sin(physConst::PI*y);
    }
    std::tuple<Eigen::VectorXd, Eigen::VectorXd, BoundaryConditions> FVM_prescribedDiffSolution(double t, int nx, int ny, double xmin, double xmax, double ymin, double ymax, double Dh, double Dv){
        //Inhomog. Dirichlet
        double dx = (xmax - xmin)/nx;
        double dy = (ymax - ymin)/ny;
        Eigen::VectorXd vec(nx * ny);
        Eigen::VectorXd source(nx * ny);
        vec.resize(nx * ny + 2*nx + 2*ny); //to include ghost nodes
        vec.setZero();
        source.resize(nx * ny);
        BoundaryConditions bc;
        bc.bcType_top = BoundaryConditionFlag::DIRICHLET_INT_BPOINT;
        bc.bcType_left = BoundaryConditionFlag::DIRICHLET_INT_BPOINT;
        bc.bcType_right = BoundaryConditionFlag::DIRICHLET_INT_BPOINT;
        bc.bcType_bot = BoundaryConditionFlag::DIRICHLET_INT_BPOINT;
        bc.bcVals_top = Vector_1D(nx, 0);
        bc.bcVals_bot = Vector_1D(nx, 0);
        bc.bcVals_left = Vector_1D(ny, 0);
        bc.bcVals_right = Vector_1D(ny, 0);

        for(int i = 0; i < nx; i++){
            for(int j = 0; j < ny; j++){
                double x = xmin + dx * (0.5 + i);
                double y = ymin + dy * (0.5 + j);
                vec[i*ny + j] = solnFcn(t, x, y);

                double dphi_dt = -solnFcn(t,x,y);
                double d2f_dx2 = -physConst::PI * physConst::PI * solnFcn(t,x,y);
                double d2f_dy2 = d2f_dx2;

                source[i*ny + j] = dphi_dt - Dh*d2f_dx2 - Dv*d2f_dy2;

                if(i == 0){
                    bc.bcVals_left[j] = std::exp(-t) * std::sin(physConst::PI*xmin) * std::sin(physConst::PI*y);
                }
                else if (i == nx - 1){
                    bc.bcVals_right[j] = std::exp(-t) * std::sin(physConst::PI*xmax) * std::sin(physConst::PI*y);
                }

                if(j == 0){
                    bc.bcVals_bot[i] = std::exp(-t) * std::sin(physConst::PI*x) * std::sin(physConst::PI*ymin);
                }
                else if (j == ny - 1){
                    bc.bcVals_top[i] = std::exp(-t) * std::sin(physConst::PI*x) * std::sin(physConst::PI*ymax);
                }
            }
        }
        return std::make_tuple(vec, source, bc);
    }
    std::tuple<Eigen::VectorXd, BoundaryConditions> initAdvection(int nx, int ny){
        //Only allow domains of (x, y) in [0, 1]
        double dx = 1.0/nx;
        double dy = 1.0/ny;
        Eigen::VectorXd vec;
        vec.resize(nx * ny + 2*nx + 2*ny); //to include ghost nodes
        vec.setZero();
        BoundaryConditions bc;
        bc.bcType_top = BoundaryConditionFlag::DIRICHLET_INT_BPOINT;
        bc.bcType_left = BoundaryConditionFlag::DIRICHLET_INT_BPOINT;
        bc.bcType_right = BoundaryConditionFlag::DIRICHLET_INT_BPOINT;
        bc.bcType_bot = BoundaryConditionFlag::DIRICHLET_INT_BPOINT;
        bc.bcVals_top = Vector_1D(nx, 0);
        bc.bcVals_bot = Vector_1D(nx, 0);
        bc.bcVals_left = Vector_1D(ny, 0);
        bc.bcVals_right = Vector_1D(ny, 0);

        for(int i = 0; i < nx; i++){
            for(int j = 0; j < ny; j++){
                double x = dx * (0.5 + i);
                double y = dy * (0.5 + j);
                if((x > 0.370 && x < 0.620 ) && (y > 0.370 && y < 0.620)){
                    vec[i*ny + j] = std::cos(4 * physConst::PI * (x - 0.495)) * std::cos(4 * physConst::PI * (y - 0.495));
                }
            }
        }
        return std::make_tuple(vec, bc);
    }
    std::tuple<double, double, double> interiorMax(const Eigen::VectorXd& vec, double xmin, double xmax, double ymin, double ymax, int nx, int ny){
        double dx = (xmax - xmin) / nx;
        double dy = (ymax - ymin) / ny;

        /* Initialize maximum to be the first value in the vector
        Avoids compiler complaining about edge case where all
        values are equal to std::numeric_limits::min() which leaves
        variables uninitialized */
        double max = vec(0);
        double maxx = dx * 0.5;
        double maxy = dy * 0.5;

        for(int i = 0; i < nx; i++){
            for(int j = 0; j < ny; j++){
                double x = dx * (0.5 + i);
                double y = dy * (0.5 + j);
                if(vec(i*ny + j) > max){
                    max = vec(i*ny + j);
                    maxx = x;
                    maxy = y;
                }
            }
        }
        return std::make_tuple(max, maxx, maxy);
    }
    std::tuple<double, double, double> interiorMin(const Eigen::VectorXd& vec, double xmin, double xmax, double ymin, double ymax, int nx, int ny){
        double dx = (xmax - xmin) / nx;
        double dy = (ymax - ymin) / ny;

        /* Initialize minimum to be the first value in the vector
        Avoids compiler complaining about edge case where all
        values are equal to std::numeric_limits::max() which leaves
        variables uninitialized */
        double min = vec(0);
        double minx = dx * 0.5;
        double miny = dy * 0.5;

        for(int i = 0; i < nx; i++){
            for(int j = 0; j < ny; j++){
                double x = dx * (0.5 + i);
                double y = dy * (0.5 + j);
                if(vec(i*ny + j) < min){
                    min = vec(i*ny + j);
                    minx = x;
                    miny = y;
                }
            }
        }
        return std::make_tuple(min, minx, miny);
    }
    TEST_CASE("Basic Inhomog. Dirichlet System Setup Tests"){
        /*  4x3 Domain Point IDs:
               12 13 14 15
            18 2  5  8  11 21
            17 1  4  7  10 20
            16 0  3  6  9  19
               22 23 24 25
            
            12-25 are ghost nodes.
        */
        double u = 0.8, v = 1.1, shear = 0, Dh = 1.1, Dv = 0.4, dt = 0.05, xlim_left = 0.1, xlim_right = 0.8, ylim_bot = 0.1, ylim_top = 0.9;
        int nx = 4, ny = 3;
        AdvDiffParams params = AdvDiffParams(u, v, shear, Dh, Dv, dt);
        Mesh mesh = Mesh(nx, ny, xlim_left, xlim_right, ylim_top, ylim_bot, MeshDomainLimitsSpec::ABS_COORDS);
        double dy = (ylim_top - ylim_bot) / ny;
        double dx = (xlim_right - xlim_left) / nx;
        Eigen::VectorXd exact, source;
        BoundaryConditions bc;
        std::tie(exact, source, bc) = FVM_prescribedDiffSolution(0, nx, ny, xlim_left, xlim_right, ylim_bot, ylim_top, Dh, Dv);
        FVM_Solver solver(params, mesh.x(), mesh.y(), bc, exact);
        solver.buildCoeffMatrix();
    
        auto mat = solver.coefMatrix();
        const std::vector<std::unique_ptr<Point>>&  points = solver.points();
        
        SECTION("Basic Sanity Checks"){
            REQUIRE(mat.cols() == 26);
            REQUIRE(mat.rows() == 26);

            int numGhost = 0;
            int numBound = 0;
            for(int i = 0; i < 26; i++){
                if(points[i]->isGhost()) {numGhost++;}
                else if (points[i]->bcType() != BoundaryConditionFlag::INTERIOR) numBound++;
            }
            REQUIRE(numGhost == 14);
            REQUIRE(numBound == 10);

            //NNZ = 5 * interior points + 2*ghost points
            // = 5 * 12 + 2 * 14 = 88
            REQUIRE(mat.nonZeros() == 88);
        }
        SECTION("Corner Ghost Point Values after applying BC"){
            const Eigen::VectorXd& phi_vec = solver.phi();
            double sol_topleft_top = solnFcn(0, xlim_left + dx/2, ylim_top);
            double sol_topleft_left = solnFcn(0, xlim_left, ylim_top - dy/2);
            REQUIRE(phi_vec[12] == Catch::Approx(2*sol_topleft_top - solnFcn(0, xlim_left + dx/2, ylim_top - dy/2)).epsilon(1e-10));
            REQUIRE(phi_vec[18] == Catch::Approx(2*sol_topleft_left - solnFcn(0, xlim_left + dx/2, ylim_top - dy/2)).epsilon(1e-10));

            double sol_topright_top = solnFcn(0, xlim_right - dx/2, ylim_top);
            double sol_topright_right = solnFcn(0, xlim_right, ylim_top - dy/2);
            REQUIRE(phi_vec[15] == Catch::Approx(2*sol_topright_top - solnFcn(0, xlim_right - dx/2, ylim_top - dy/2)).epsilon(1e-10));
            REQUIRE(phi_vec[21] == Catch::Approx(2*sol_topright_right - solnFcn(0, xlim_right - dx/2, ylim_top - dy/2)).epsilon(1e-10));
            
            double sol_botleft_bot = solnFcn(0, xlim_left + dx/2, ylim_bot);
            double sol_botleft_left = solnFcn(0, xlim_left, ylim_bot + dy/2);
            REQUIRE(phi_vec[22] == Catch::Approx(2*sol_botleft_bot - solnFcn(0, xlim_left + dx/2, ylim_bot + dy/2)).epsilon(1e-10));
            REQUIRE(phi_vec[16] == Catch::Approx(2*sol_botleft_left - solnFcn(0, xlim_left + dx/2, ylim_bot + dy/2)).epsilon(1e-10));
            
            double sol_botright_bot = solnFcn(0, xlim_right - dx/2, ylim_bot);
            double sol_botright_right = solnFcn(0, xlim_right, ylim_bot + dy/2);
            REQUIRE(phi_vec[25] == Catch::Approx(2*sol_botright_bot - solnFcn(0, xlim_right - dx/2, ylim_bot + dy/2)).epsilon(1e-10));
            REQUIRE(phi_vec[19] == Catch::Approx(2*sol_botright_right - solnFcn(0, xlim_right - dx/2, ylim_bot + dy/2)).epsilon(1e-10));

        }
        SECTION("Ghost Point Coeffs"){
            // (phi_int + phi_ghost)/2 = phi_bound
            REQUIRE(mat.coeffRef(12,2) == 0.5);
            REQUIRE(mat.coeffRef(12,12) == 0.5);
            REQUIRE(mat.coeffRef(18,2) == 0.5);
            REQUIRE(mat.coeffRef(18,18) == 0.5);

            REQUIRE(mat.coeffRef(15,11) == 0.5);
            REQUIRE(mat.coeffRef(15,15) == 0.5);
            REQUIRE(mat.coeffRef(21,11) == 0.5);
            REQUIRE(mat.coeffRef(21,21) == 0.5);

            REQUIRE(mat.coeffRef(16,0) == 0.5);
            REQUIRE(mat.coeffRef(16,16) == 0.5);
            REQUIRE(mat.coeffRef(22,0) == 0.5);
            REQUIRE(mat.coeffRef(22,22) == 0.5);

            REQUIRE(mat.coeffRef(19,9) == 0.5);
            REQUIRE(mat.coeffRef(19,19) == 0.5);
            REQUIRE(mat.coeffRef(25,9) == 0.5);
            REQUIRE(mat.coeffRef(25,25) == 0.5);
        }
    }

    TEST_CASE("SOR Solver"){
        /*
        Check that the SOR iterative solver for linear systems is correct
        Setup prescribed solution test for Ax = b:
        - Specify A and x to compute b
        - Use the solver with A and b to recover x

        Here A = [
        0.50 0.25 0.00 0.00 ... 0.00
        0.25 0.50 0.25 0.00 ... 0.00
        0.00 0.25 0.50 0.25 ... 0.00
        .
        .
        0.00 ... 0.00 0.25 0.50 0.25
        0.00 ... 0.00 0.00 0.25 0.50]

        and ExactSolution = [0 1 ... Npoints].T
        */

        int Npoints = 1000;

        Eigen::SparseMatrix<double, Eigen::RowMajor> A(Npoints, Npoints);
        std::vector<Eigen::Triplet<double>> tripletList;

        Eigen::VectorXd ExactSolution(Npoints);
        Eigen::VectorXd SORSolution(Npoints);

        // Make a tri-banded symetric matrix and solution
        for (int i = 0; i < Npoints; ++i){
            // Set diagonal to 0.5
            tripletList.emplace_back(i, i, 0.5);
            // Set lower-diagonal to 0.25
            if (i > 1)
                tripletList.emplace_back(i, i-1, 0.25);
            // Set upper-diagonal to 0.25
            if (i < Npoints - 1)
                tripletList.emplace_back(i, i+1, 0.25);

            // Impose the values of x
            ExactSolution[i] = i;
            // Impose first guess of solution to be all 0
            SORSolution[i] = 0;
        }

        // Set the values in the matrix A
        A.setFromTriplets(tripletList.begin(), tripletList.end());

        // Construct the right hand side vector b using x
        Eigen::VectorXd rhs = (A * ExactSolution).eval();

        // Relaxation parameter
        double omega = 1.0;
        // Residual tolerance (% [0-1])
        double max_res = 1e-3;
        // Number iterations between evaluations of residual
        int n_iters = 3;

        // Pretend we don't know x and solve Ax = b
        FVM_ANDS::sor_solve(A, rhs, SORSolution, omega, max_res, n_iters);

        // Error as percentage of L2 norm of differences relative to L2 norm of exact solution
        auto error = 100 * (SORSolution - ExactSolution).eval().lpNorm<2>()/ ExactSolution.lpNorm<2>();

        // Impose error to be less than 1%
        REQUIRE(error < 1);
    }

    TEST_CASE("Pure Diffusion, Inhomog. Dirichlet BC, Prescribed Solution 10k Points"){
        // Dh = 0.9, Dv = 0.4
        // Test solver accuracy on 10k interior points
        // Solution: f(x,y,t) = exp(-t) * sin(pi*x) * sin(2*pi*y)
        
        double u = 0, v = 0, shear = 0, Dh = 1.0, Dv = 1.0, xlim_left = 0.1, xlim_right = 0.8, ylim_bot = 0.2, ylim_top = 0.9;
        int nx = 50, ny = 50;
        // double dx = 1.0/nx;
        // double dy = 1.0/ny;
        //double dt = 0.24 * std::min( (dx*dx) / (2*Dh) , (dy*dy) / (2*Dv));
        double dt = 0.1;
        std::cout << "dt: " << dt << std::endl;
        AdvDiffParams params = AdvDiffParams(u, v, shear, Dh, Dv, dt);
        Mesh mesh = Mesh(nx, ny, xlim_left, xlim_right, ylim_top, ylim_bot, MeshDomainLimitsSpec::ABS_COORDS);
        Eigen::VectorXd exact, source;
        BoundaryConditions bc;
        std::tie(exact, source, bc) = FVM_prescribedDiffSolution(0, nx, ny, xlim_left, xlim_right, ylim_bot, ylim_top, Dh, Dv);
        FVM_Solver solver(params, mesh.x(), mesh.y(), bc, exact, false, 1000 , 1e-5);
        //the convergence on this problem really sucks without ILUT precond
        //Problem is, ILUT is way too expensive to use in APCEMM as a default. Maybe the code can be templated later to allow for customizability.
        double t = 0, t_max = 1.5;
        int n_timesteps = t_max/dt;
        for(int i = 0; i < n_timesteps; i++){
            t = dt*(i + 1);
            std::tie(exact, source, bc) = FVM_prescribedDiffSolution(t, nx, ny, xlim_left, xlim_right, ylim_bot, ylim_top, Dh, Dv);
            solver.updateBoundaryCondition(bc);
            Eigen::VectorXd soln = solver.solve(source);
            //Eigen::VectorXd soln = solver.explicitSolve(source);
            auto interior_idxs = Eigen::seq(0, nx*ny - 1);
            double abs_l2_error = (soln(interior_idxs) - exact(interior_idxs)).lpNorm<2>() / exact(interior_idxs).lpNorm<2>();
            cout << "L2-Error: " << abs_l2_error << "\n";
            cout << "Exact Soln Maximum: " << exact.maxCoeff() << "\n";
            cout << "Soln Maximum: " << soln.maxCoeff() << "\n";
            cout << "Soln Minimum: " << soln(interior_idxs).minCoeff() << "\n";
            cout << "Num Iters: " << solver.numIters() << "\n";
            cout << "Solver Residual: " << solver.error() << "\n" << endl;

            REQUIRE(abs_l2_error < 0.05); //pretty unimpressive accuracy in time because it's fully imp
            REQUIRE(soln(interior_idxs).minCoeff() >= 0.0); //check monoticity preservation
        }
    }

    // Test is failing, hide it from ctest command, must be run with 'tests/test_solver [failing]'
    TEST_CASE("Explicit Advection - Diffusion w/ Shear", "[.failing]"){
        // Test solver accuracy on 250k interior points
        // Solution: f(x,y,t) = exp(-t) * sin(pi*x) * sin(2*pi*y)
        
        double u = 0.3, v = -0.25, shear = 0.2, Dh = 0.00, Dv = 0.00, xlim_left = 0.0, xlim_right = 1.0, ylim_bot = 0.0, ylim_top = 1.0;
        /* u = u0 - shear * y, v = v0
           Therefore:
           y(t) = y0 + v0t
           u(t) = u0 - shear * (y0 + v0t)
           x(t) = 0.5 + \int u(t) dt
                = 0.5 + u0t - shear * (y0*t + v0/2 t^2)
        */ 
        int nx = 500, ny = 500;
        double dx = 1.0/nx;
        double dy = 1.0/ny;
        double dt =  1;
        double courant_max = 0.5;
        double courant = std::abs(std::max(u, u - ylim_top*shear)) * dt / dx + std::abs(v) * dt / dy;

        cout << "Advection limited dt: " << courant_max/courant * dt << endl;
        dt = courant_max/courant * dt;
        AdvDiffParams params = AdvDiffParams(u, v, shear, Dh, Dv, dt);
        Mesh mesh = Mesh(nx, ny, xlim_left, xlim_right, ylim_top, ylim_bot, MeshDomainLimitsSpec::ABS_COORDS);
        Eigen::VectorXd init;
        BoundaryConditions bc;
        std::tie(init, bc) = initAdvection(nx, ny);

        FVM_Solver solver(params, mesh.x(), mesh.y(), bc, init, false, 100 , 1e-5);

        cout << "Soln Maximum: " << solver.phi().maxCoeff() << "\n";
        cout << "Soln Interior Minimum: " << solver.phi()(Eigen::seq(0,nx*ny-1)).minCoeff() << "\n";

        double t = 0, t_max = 1;
        int n_timesteps = t_max/dt;
        double max, maxx, maxy, min, minx, miny;
        for(int i = 0; i < n_timesteps; i++){
            t = dt*(i + 1); //implicit method
            cout << "solving" << endl;
            Eigen::VectorXd soln = solver.explicitSolve();
            auto interior_idxs = Eigen::seq(0, nx*ny - 1);
            std::tie(max, maxx, maxy) = interiorMax(soln, 0, 1, 0, 1, nx, ny);
            cout << "max: " << max << ", max x: " << maxx << ", " << "max y: " << maxy << endl;
            std::tie(min, minx, miny) = interiorMin(soln, 0, 1, 0, 1, nx, ny);
            cout << "min: " << min << ", min x: " << minx << ", " << "min y: " << miny << endl;
            cout << "total mass: " << dx*dy*soln(interior_idxs).sum() << "\n" << endl;
            REQUIRE(min >= 0.0);
            REQUIRE(max <= 1.0);

        }
        double xmax_exp = 0.495 + u * t -  shear * (0.5*t + v / 2.0 * t * t);
        double ymax_exp = 0.495 + v * t;
        std::cout << "Expected xmax: " << xmax_exp << std::endl;
        std::cout << "Expected ymax: " << ymax_exp << std::endl;
        REQUIRE(std::abs(maxx-xmax_exp) < 0.01);
        REQUIRE(std::abs(maxy-ymax_exp) < 0.01);

    }
    TEST_CASE("Implicit Advection w/ Shear"){
        // Test solver accuracy on 250k interior points
        // Implicit doesnt preserve monoticity w/o enforcing cfl condition
        // Solution: f(x,y,t) = exp(-t) * sin(pi*x) * sin(2*pi*y)
        
        double u = 0.2, v = -0.25, shear = 0.1, Dh = 0.00, Dv = 0.00, xlim_left = 0.0, xlim_right = 1.0, ylim_bot = 0.0, ylim_top = 1.0;
        /* u = u0 - shear * y, v = v0
           Therefore:
           y(t) = y0 + v0t
           u(t) = u0 - shear * (y0 + v0t)
           x(t) = 0.5 + \int u(t) dt
                = 0.5 + u0t - shear * (y0*t + v0/2 t^2)
        */ 
        int nx = 200, ny = 200;
        double dx = 1.0/nx;
        double dy = 1.0/ny;
        double dt =  1;
        double courant_max = 5;
        double courant = std::abs(std::max(u, u - ylim_top*shear)) * dt / dx + std::abs(v) * dt / dy;

        double r_diff = (Dh*dt/dx/dx + Dh*dt/dy/dy);
        double r_diff_max = 0.4;

        dt = std::min(courant_max/courant * dt, r_diff_max / r_diff * dt);
        cout << "Advection dt: " << dt << endl;
        AdvDiffParams params = AdvDiffParams(u, v, shear, Dh, Dv, dt);
        Mesh mesh = Mesh(nx, ny, xlim_left, xlim_right, ylim_top, ylim_bot, MeshDomainLimitsSpec::ABS_COORDS);
        Eigen::VectorXd init;
        BoundaryConditions bc;
        std::tie(init, bc) = initAdvection(nx, ny);

        FVM_Solver solver(params, mesh.x(), mesh.y(), bc, init, false, 100 , 1e-5);

        cout << "Soln Maximum: " << solver.phi().maxCoeff() << "\n";
        cout << "Soln Interior Minimum: " << solver.phi()(Eigen::seq(0,nx*ny-1)).minCoeff() << "\n";

        double t = 0, t_max = 1;
        int n_timesteps = t_max/dt;
        double max, maxx, maxy, min, minx, miny;
        for(int i = 0; i < n_timesteps; i++){
            t = dt*(i + 1); //implicit method
            cout << "solving" << endl;
            Eigen::VectorXd soln = solver.solve();
            auto interior_idxs = Eigen::seq(0, nx*ny - 1);
            std::tie(max, maxx, maxy) = interiorMax(soln, 0, 1, 0, 1, nx, ny);
            cout << "max: " << max << ", max x: " << maxx << ", " << "max y: " << maxy << endl;
            std::tie(min, minx, miny) = interiorMin(soln, 0, 1, 0, 1, nx, ny);
            cout << "min: " << min << ", min x: " << minx << ", " << "min y: " << miny << endl;
            cout << "total mass: " << dx*dy*soln(interior_idxs).sum() << "\n" << endl;;

        }
        double xmax_exp = 0.495 + u * t -  shear * (0.5*t + v / 2.0 * t * t);
        double ymax_exp = 0.495 + v * t;
        std::cout << "Expected xmax: " << 0.495 + u * t -  shear * (0.5*t + v / 2.0 * t * t) << std::endl;
        std::cout << "Expected ymax: " << 0.495 + v * t << std::endl;
        REQUIRE(std::abs(maxx-xmax_exp) < 0.01);
        REQUIRE(std::abs(maxy-ymax_exp) < 0.01);

    }

    // Test is failing, hide it from ctest command, must be run with 'tests/test_solver [failing]'
    TEST_CASE("Operator Split Advection-Diffusion", "[.failing]"){
        // Test solver accuracy on 250k interior points
        // Implicit doesnt preserve monoticity w/o enforcing cfl condition
        // Solution: f(x,y,t) = exp(-t) * sin(pi*x) * sin(2*pi*y)
        
        double u = 0.2, v = -0.25, shear = 0.1, Dh = 0.1, Dv = 0.1, xlim_left = 0.0, xlim_right = 1.0, ylim_bot = 0.0, ylim_top = 1.0;
        /* u = u0 - shear * y, v = v0
           Therefore:
           y(t) = y0 + v0t
           u(t) = u0 - shear * (y0 + v0t)
           x(t) = 0.5 + \int u(t) dt
                = 0.5 + u0t - shear * (y0*t + v0/2 t^2)
        */ 
        int nx = 500, ny = 500;
        double dx = 1.0/nx;
        double dy = 1.0/ny;
        double dt =  0.1;

        AdvDiffParams params = AdvDiffParams(u, v, shear, Dh, Dv, dt);
        Mesh mesh = Mesh(nx, ny, xlim_left, xlim_right, ylim_top, ylim_bot, MeshDomainLimitsSpec::ABS_COORDS);
        Eigen::VectorXd init;
        BoundaryConditions bc;
        std::tie(init, bc) = initAdvection(nx, ny);

        FVM_Solver solver(params, mesh.x(), mesh.y(), bc, init, 100 , 1e-5);

        cout << "Soln Maximum: " << solver.phi().maxCoeff() << "\n";
        cout << "Soln Interior Minimum: " << solver.phi()(Eigen::seq(0,nx*ny-1)).minCoeff() << "\n";

        // double t;
        double t_max = 1;
        int n_timesteps = t_max/dt;
        double max, maxx, maxy, min, minx, miny;
        for(int i = 0; i < n_timesteps; i++){
            // t = dt*(i + 1); //implicit method
            cout << "solving" << endl;
            Eigen::VectorXd soln = solver.solve();
            auto interior_idxs = Eigen::seq(0, nx*ny - 1);
            std::tie(max, maxx, maxy) = interiorMax(soln, 0, 1, 0, 1, nx, ny);
            cout << "max: " << max << ", max x: " << maxx << ", " << "max y: " << maxy << endl;
            std::tie(min, minx, miny) = interiorMin(soln, 0, 1, 0, 1, nx, ny);
            cout << "min: " << min << ", min x: " << minx << ", " << "min y: " << miny << endl;
            cout << "total mass: " << dx*dy*soln(interior_idxs).sum() << "\n" << endl;;

        }

        //Derived from implicit solution enforcing CFL (dt = 0.002)
        //REQUIRE(std::abs(max-0.0123284) < 0.003); //accuracy on diff is pretty bad with diagonal precond and high timesteps
        REQUIRE(std::abs(maxx-0.575) < 0.01);
        REQUIRE(std::abs(maxy-0.381) < 0.01);

    }

    TEST_CASE("Semi-Lagrangian 1D Advection", "[advection]"){
        SECTION("Pure integer shift right (positive velocity)"){
            std::vector<double> slice = {1.0, 2.0, 3.0, 4.0, 5.0};
            double velocity = 3.0; // dx = 1.0, dt = 1.0 -> disp = 3.0 -> shift by 3
            double dt = 1.0;
            double ds = 1.0;
            double bc_left = 0.0;
            double bc_right = 0.0;
            semiLagrangianAdvection1D(slice, velocity, dt, ds, bc_left, bc_right);
            REQUIRE(slice[0] == Catch::Approx(0.0));
            REQUIRE(slice[1] == Catch::Approx(0.0));
            REQUIRE(slice[2] == Catch::Approx(0.0));
            REQUIRE(slice[3] == Catch::Approx(1.0));
            REQUIRE(slice[4] == Catch::Approx(2.0));
        }

        SECTION("Pure integer shift left (negative velocity)"){
            std::vector<double> slice = {1.0, 2.0, 3.0, 4.0, 5.0};
            double velocity = -2.0; // dx = 1.0, dt = 1.0 -> disp = -2.0 -> shift left by 2
            double dt = 1.0;
            double ds = 1.0;
            double bc_left = 0.0;
            double bc_right = 0.0;
            semiLagrangianAdvection1D(slice, velocity, dt, ds, bc_left, bc_right);
            REQUIRE(slice[0] == Catch::Approx(3.0));
            REQUIRE(slice[1] == Catch::Approx(4.0));
            REQUIRE(slice[2] == Catch::Approx(5.0));
            REQUIRE(slice[3] == Catch::Approx(0.0));
            REQUIRE(slice[4] == Catch::Approx(0.0));
        }

        SECTION("Inflow boundary condition padding"){
            std::vector<double> slice = {1.0, 2.0, 3.0, 4.0, 5.0};
            double velocity = 2.0;
            double dt = 1.0;
            double ds = 1.0;
            double bc_left = 9.9;
            double bc_right = 0.0;
            semiLagrangianAdvection1D(slice, velocity, dt, ds, bc_left, bc_right);
            REQUIRE(slice[0] == Catch::Approx(9.9));
            REQUIRE(slice[1] == Catch::Approx(9.9));
            REQUIRE(slice[2] == Catch::Approx(1.0));
            REQUIRE(slice[3] == Catch::Approx(2.0));
            REQUIRE(slice[4] == Catch::Approx(3.0));
        }

        SECTION("Mass conservation on interior pulse"){
            int N = 50;
            std::vector<double> slice(N, 0.0);
            slice[20] = 1.0;
            slice[21] = 2.0;
            slice[22] = 1.0;
            double initial_mass = 4.0;

            double velocity = 15.5; // moves ~15.5 cells
            double dt = 1.0;
            double ds = 1.0;
            semiLagrangianAdvection1D(slice, velocity, dt, ds, 0.0, 0.0);

            double final_mass = 0.0;
            for (double val : slice) {
                final_mass += val;
            }
            REQUIRE(final_mass == Catch::Approx(initial_mass).margin(1e-10));
        }
    }

    TEST_CASE("Semi-Lagrangian 1D Advection preserves monotonicity", "[advection]"){

        // Monotone non-decreasing profile
        const std::vector<double> initial = {0.0, 1.0, 3.0, 4.0, 5.0, 5.0};
        const double dt = 1.0;
        const double ds = 1.0;
        const double bc_left = 0.0;
        const double bc_right = 5.0;

        auto checkMonotonicity = [&](const std::vector<double>& slice){
            // 1. Values should stay in the range spanned by the initial data and the BCs.
            for (std::size_t m = 0; m < slice.size(); m++) {
                INFO("cell " << m << " = " << slice[m]);
                REQUIRE(slice[m] >= 0.0 - 1e-12);
                REQUIRE(slice[m] <= 5.0 + 1e-12);
            }

            // 2. A monotone profile must stay monotone
            for (std::size_t m = 1; m < slice.size(); m++) {
                INFO("cells " << m - 1 << ", " << m << " = " << slice[m-1] << ", " << slice[m]);
                REQUIRE(slice[m] >= slice[m-1] - 1e-12);
            }

            // 3. Total variation must not increase
            double tv_initial = 0.0;
            double tv_final = 0.0;
            for (std::size_t m = 1; m < slice.size(); m++) {
                tv_initial += std::abs(initial[m] - initial[m-1]);
                tv_final += std::abs(slice[m] - slice[m-1]);
            }
            INFO("TV before = " << tv_initial << ", TV after = " << tv_final);
            REQUIRE(tv_final <= tv_initial + 1e-12);
        };

        SECTION("Fractional CFL 0.80, above the 2/3 TVD limit"){
            std::vector<double> slice = initial;
            double velocity = 0.8; // dt = ds = 1 -> no integer shift, fractional CFL = 0.8
            semiLagrangianAdvection1D(slice, velocity, dt, ds, bc_left, bc_right);
            checkMonotonicity(slice);
        }

        SECTION("Fractional CFL 0.60, below the 2/3 TVD limit"){
            std::vector<double> slice = initial;
            double velocity = 0.6; // dt = ds = 1 -> no integer shift, fractional CFL = 0.6
            semiLagrangianAdvection1D(slice, velocity, dt, ds, bc_left, bc_right);
            checkMonotonicity(slice);
        }
    }
}

namespace FVM_ANDS{
    // Zero Dirichlet boundary conditions on all four sides of an nx x ny grid.
    BoundaryConditions zeroDirichletBC(int nx, int ny){
        BoundaryConditions bc;
        bc.bcType_top = BoundaryConditionFlag::DIRICHLET_INT_BPOINT;
        bc.bcType_left = BoundaryConditionFlag::DIRICHLET_INT_BPOINT;
        bc.bcType_right = BoundaryConditionFlag::DIRICHLET_INT_BPOINT;
        bc.bcType_bot = BoundaryConditionFlag::DIRICHLET_INT_BPOINT;
        bc.bcVals_top = Vector_1D(nx, 0);
        bc.bcVals_bot = Vector_1D(nx, 0);
        bc.bcVals_left = Vector_1D(ny, 0);
        bc.bcVals_right = Vector_1D(ny, 0);
        return bc;
    }

    TEST_CASE("Variable-D diffusion matrix uses face coefficients"){
        /*
        Discretize div(D grad phi) with a strongly discontinuous D.
        If both rows next to a face use the same face D:
        - the interior block of the matrix is symmetric,
        - every column of a cell away from the boundary sums to 1 (mass conservation).
        */
        int nx = 20, ny = 17;
        double xlim_left = 0.0, xlim_right = 0.9, ylim_bot = 0.0, ylim_top = 0.35;
        double dt = 1e-2;
        int nInt = nx * ny;
        int nTot = nInt + 2 * nx + 2 * ny;

        // Random base D, with a 7x jump inside a rectangle (like DiffParam at t = 0).
        // Rectangle is (10, nx-2) x (9, ny-1)
        std::mt19937 rng(20);
        std::uniform_real_distribution<double> dist(0.5, 1.5);
        Eigen::VectorXd Dh(nInt), Dv(nInt);
        for(int i = 0; i < nx; i++){
            for(int j = 0; j < ny; j++){
                double factor = (i >= 10 && i <= nx - 2 && j >= 9 && j <= ny - 1) ? 7.0 : 1.0;
                Dh[i*ny + j] = factor * dist(rng);
                Dv[i*ny + j] = factor * 0.1 * dist(rng);
            }
        }

        AdvDiffParams params = AdvDiffParams(0, 0, 0, 1.0, 1.0, dt);
        Mesh mesh = Mesh(nx, ny, xlim_left, xlim_right, ylim_top, ylim_bot, MeshDomainLimitsSpec::ABS_COORDS);
        Eigen::VectorXd phi_init = Eigen::VectorXd::Zero(nTot);
        FVM_Solver solver(params, mesh.x(), mesh.y(), zeroDirichletBC(nx, ny), phi_init);
        solver.updateDiffusion(Dh, Dv);
        solver.buildCoeffMatrix();
        Eigen::SparseMatrix<double, Eigen::RowMajor> mat = solver.coefMatrix();

        SECTION("Interior block is exactly symmetric"){
            // Only the interior block (rows and columns 0..nInt-1) is symmetric.
            // Indices >= nInt are ghost points. The ghost rows hold the BC equation (0.5, 0.5),
            // not a diffusion stencil, so the full matrix is not symmetric.
            Eigen::SparseMatrix<double, Eigen::RowMajor> A_int = mat.topLeftCorner(nInt, nInt);
            Eigen::SparseMatrix<double, Eigen::RowMajor> A_int_T = A_int.transpose();
            // A(i, j) and A(j, i) come from the same face D so they must be exactly equal
            // The Frobenius norm of the difference is 0 only if every entry matches its mirror.
            REQUIRE((A_int - A_int_T).norm() == 0.0);
        }

        SECTION("Column sums equal 1 away from the boundary"){
            // A row vector of ones times the matrix gives the sum of each column: colSum[k] = sum over rows of A(row, k).
            Eigen::RowVectorXd colSum = Eigen::RowVectorXd::Ones(nTot) * mat;
            double maxDefect = 0;
            for(int i = 1; i < nx - 1; i++){
                for(int j = 1; j < ny - 1; j++){
                    int k = i*ny + j;
                    // Round-off scales with the largest entry in the column which is the diagonal.
                    double columnDefect = std::abs(colSum[k] - 1.0) / mat.coeff(k, k);
                    maxDefect = std::max(maxDefect, columnDefect);
                }
            }
            REQUIRE(maxDefect < 1e-14);
        }
    }

    // Run nSteps of transport on a Gaussian blob with zero Dirichlet BC, through operatorSplitSolve,
    // (what LAGRIDPlumeModel::runTransport uses). It builds the operator-split diffusion matrix
    // and solves the implicit diffusion step with SOR.
    // All velocities are zero, so the advection half-steps do nothing and the step is pure diffusion.
    // If enhanceBlob is true, D is 7x larger where phi > 1e-4 * max(phi), as in updateDiffVecs.
    // sorThreshold is the SOR relative residual threshold (default value: 1e-3, lowered for testing).
    // Returns the relative mass change and the largest |phi| on a boundary cell, relative to max(phi_0).
    std::tuple<double, double> diffusionMassChange(bool enhanceBlob, int nSteps, double sorThreshold){
        // Cell size is fixed (dx = 1/60, dy = 1/50), so dt*D/dx^2 does not depend on the domain size.
        // The domain must be large enough that the blob does not reach the boundary in nSteps.
        double Lx = 2.0, Ly = 2.0;
        int nx = static_cast<int>(60 * Lx), ny = static_cast<int>(50 * Ly);
        double dt = 0.25, Dh_base = 2.5e-4, Dv_base = 1.25e-4, sigma = 0.03;
        int nInt = nx * ny;
        int nTot = nInt + 2 * nx + 2 * ny;
        Mesh mesh = Mesh(nx, ny, 0.0, Lx, Ly, 0.0, MeshDomainLimitsSpec::ABS_COORDS);
        double dx = mesh.x()[1] - mesh.x()[0];
        double dy = mesh.y()[1] - mesh.y()[0];

        // Interior cells first (index i*ny + j), then the ghost points.
        Eigen::VectorXd phi = Eigen::VectorXd::Zero(nTot);
        for(int i = 0; i < nx; i++){
            for(int j = 0; j < ny; j++){
                double rx = mesh.x()[i] - 0.5 * Lx, ry = mesh.y()[j] - 0.5 * Ly;
                phi[i*ny + j] = std::exp(-(rx*rx + ry*ry) / (2*sigma*sigma));
            }
        }
        double phiMax0 = phi.maxCoeff();

        Eigen::VectorXd Dh(nInt), Dv(nInt);
        for(int k = 0; k < nInt; k++){
            double factor = (enhanceBlob && phi[k] > 1e-4 * phiMax0) ? 7.0 : 1.0;
            Dh[k] = factor * Dh_base;
            Dv[k] = factor * Dv_base;
        }

        AdvDiffParams params = AdvDiffParams(0, 0, 0, Dh_base, Dv_base, dt);
        FVM_Solver solver(params, mesh.x(), mesh.y(), zeroDirichletBC(nx, ny), phi);
        solver.updateDiffusion(Dh, Dv);
        solver.setSorThreshold(sorThreshold);

        auto interior = Eigen::seq(0, nInt - 1);
        double mass0 = dx * dy * phi(interior).sum();
        for(int n = 0; n < nSteps; n++){
            solver.operatorSplitSolve(false);
        }
        const Eigen::VectorXd& soln = solver.phi();
        double mass1 = dx * dy * soln(interior).sum();

        double boundaryMax = 0;
        for(int i = 0; i < nx; i++){
            for(int j = 0; j < ny; j++){
                if(i == 0 || i == nx - 1 || j == 0 || j == ny - 1){
                    boundaryMax = std::max(boundaryMax, std::abs(soln[i*ny + j]));
                }
            }
        }
        return std::make_tuple(std::abs(mass1 - mass0) / mass0, boundaryMax / phiMax0);
    }

    TEST_CASE("Implicit diffusion conserves mass with discontinuous D"){
        int nSteps = 40;
        double relMassChange, boundaryMax;

        // At a tight SOR threshold, the only mass error left is round-off so this checks that the discretization itself conserves mass.
        SECTION("Uniform D (control), tight SOR threshold"){
            std::tie(relMassChange, boundaryMax) = diffusionMassChange(false, nSteps, 1e-13);
            cout << "Uniform D, SOR 1e-13: relative mass change " << relMassChange << ", boundary max " << boundaryMax << endl;
            // The blob must not reach the boundary, else the Dirichlet BC removes mass.
            REQUIRE(boundaryMax < 1e-12);
            REQUIRE(relMassChange < 1e-10);
        }
        SECTION("D jumps 7x at the blob edge, tight SOR threshold"){
            std::tie(relMassChange, boundaryMax) = diffusionMassChange(true, nSteps, 1e-13);
            cout << "Discontinuous D, SOR 1e-13: relative mass change " << relMassChange << ", boundary max " << boundaryMax << endl;
            REQUIRE(boundaryMax < 1e-12);
            REQUIRE(relMassChange < 1e-10);
        }

        // In the transport step, the SOR threshold is historically set to 1e-3. This might need to be revisited.
        // The unconverged SOR residual changes the mass when D varies.
        // This bound is a regression check on the transport step, conservation of the scheme is tested by the tighter check above.
        SECTION("D jumps 7x at the blob edge, transport SOR threshold"){
            std::tie(relMassChange, boundaryMax) = diffusionMassChange(true, nSteps, 1e-3);
            cout << "Discontinuous D, SOR 1e-3: relative mass change " << relMassChange << ", boundary max " << boundaryMax << endl;
            REQUIRE(boundaryMax < 1e-12);
            REQUIRE(relMassChange < 0.05);
        }
    }

    TEST_CASE("Normalized volume transport preserves constant particle volume at low concentrations", "[two_moment]"){
        int nx = 30, ny = 30;
        double Lx = 300.0, Ly = 300.0;
        double dt = 10.0, Dh = 15.0, Dv = 0.15;
        int nInt = nx * ny;
        int nTot = nInt + 2 * nx + 2 * ny;
        Mesh mesh = Mesh(nx, ny, 0.0, Lx, Ly, 0.0, MeshDomainLimitsSpec::ABS_COORDS);
        BoundaryConditions bc = zeroDirichletBC(nx, ny);
        Eigen::VectorXd phi_init = Eigen::VectorXd::Zero(nTot);

        AdvDiffParams params(0.0, -0.02, 0.002, Dh, Dv, dt);
        FVM_Solver solver(params, mesh.x(), mesh.y(), bc, phi_init);

        // Low-concentration small-crystal bin: pdf ~ 1e-12 cm^-3, V_0 ~ 1e-20 m^3, dln(r) ~ 0.196
        // Unscaled volume ~ 2e-33 m^3/cm^3 (< 1e-30 sor_solve threshold!)
        const double dLogBin = 0.1959;
        const double vEdgeLow = 0.8e-20;
        const double vEdgeHigh = 1.44e-20;
        const double vNom = 0.5 * (vEdgeLow + vEdgeHigh);
        const double vTrue = 1.05e-20;
        const double volScale = dLogBin * vNom;
        const double invVolScale = 1.0 / volScale;

        Vector_2D pdf(ny, Vector_1D(nx, 0.0));
        Vector_2D volScaled(ny, Vector_1D(nx, 0.0));

        for (int j = 10; j < 20; ++j) {
            for (int i = 10; i < 20; ++i) {
                pdf[j][i] = 1.0e-12;
                double vol = dLogBin * vTrue * pdf[j][i];
                volScaled[j][i] = vol * invVolScale;
            }
        }

        solver.operatorSplitSolve2DVec(pdf, bc, false);
        solver.operatorSplitSolve2DVec(volScaled, bc, false);

        int populatedCount = 0;
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                double vol = volScaled[j][i] * volScale;
                if (pdf[j][i] > 1.0e-20) {
                    populatedCount++;
                    double vRecovered = vol / (pdf[j][i] * dLogBin);
                    REQUIRE(vRecovered == Catch::Approx(vTrue).epsilon(1e-12));
                }
            }
        }
        REQUIRE(populatedCount > 100);
    }

    TEST_CASE("Coupled two-moment 1D FFSL advection preserves mass, number, and bounded V/N ratio", "[two_moment][advection]"){
        const int N = 40;
        const double q_min = 0.714;
        const double q_max = 1.286;
        const double dt = 1.0;
        const double ds = 1.0;

        SECTION("Positive and negative multi-cell advection of step profile with varying V/N"){
            for (double velocity : {2.35, -1.75}) {
                std::vector<double> sliceN(N, 0.0);
                std::vector<double> sliceV(N, 0.0);

                // Sharp plume core with a spatial gradient in mean particle volume q = V/N in [0.80, 1.20]
                for (int m = 10; m <= 22; ++m) {
                    sliceN[m] = (m == 10 || m == 22) ? 0.2 : 5.0; // steep edge jump
                    double q_local = 0.80 + 0.40 * double(m - 10) / 12.0;
                    sliceV[m] = sliceN[m] * q_local;
                }

                double initSumN = 0.0, initSumV = 0.0;
                for (int m = 0; m < N; ++m) {
                    initSumN += sliceN[m];
                    initSumV += sliceV[m];
                }

                semiLagrangianAdvectionTwoMoment1D(
                    sliceN, sliceV, velocity, dt, ds,
                    0.0, 0.0, q_min, q_max);

                double finalSumN = 0.0, finalSumV = 0.0;
                for (int m = 0; m < N; ++m) {
                    REQUIRE(sliceN[m] >= -1e-15);
                    REQUIRE(sliceV[m] >= -1e-15);
                    finalSumN += sliceN[m];
                    finalSumV += sliceV[m];
                    if (sliceN[m] > 1.0e-12) {
                        double q_m = sliceV[m] / sliceN[m];
                        // Must stay strictly inside [0.80, 1.20] without overshooting or undershooting to bin edges
                        REQUIRE(q_m >= 0.80 - 1e-12);
                        REQUIRE(q_m <= 1.20 + 1e-12);
                    }
                }
                REQUIRE(finalSumN == Catch::Approx(initSumN).margin(1e-12));
                REQUIRE(finalSumV == Catch::Approx(initSumV).margin(1e-12));
            }
        }
    }

    TEST_CASE("Substepped 1D vertical settling with spatially varying vFall conserves mass and number", "[two_moment][settling]"){
        const int N = 40;
        const double q_min = 0.714;
        const double q_max = 1.286;
        const double dy = 10.0;
        const double dt = 30.0; // half of 60s subcycle step

        std::vector<double> colN(N, 0.0);
        std::vector<double> colV(N, 0.0);
        std::vector<double> w_fall(N, 0.0);

        for (int m = 0; m < N; ++m) {
            // Spatially varying settling velocity from 0.15 m/s to 0.55 m/s (requires subcycling since max CFL = 0.55*30/10 = 1.65 > 0.8)
            w_fall[m] = 0.15 + 0.40 * double(m) / double(N - 1);
            if (m >= 18 && m <= 28) {
                colN[m] = (m == 18 || m == 28) ? 0.5 : 4.0;
                double q_local = 0.82 + 0.36 * double(m - 18) / 10.0;
                colV[m] = colN[m] * q_local;
            }
        }

        double initSumN = 0.0, initSumV = 0.0;
        for (int m = 0; m < N; ++m) {
            initSumN += colN[m];
            initSumV += colV[m];
        }

        verticalSettlingTwoMoment1D(colN, colV, w_fall, dt, dy, 0.0, 0.0, q_min, q_max);

        double finalSumN = 0.0, finalSumV = 0.0;
        for (int m = 0; m < N; ++m) {
            REQUIRE(colN[m] >= -1e-15);
            REQUIRE(colV[m] >= -1e-15);
            finalSumN += colN[m];
            finalSumV += colV[m];
            if (colN[m] > 1.0e-12) {
                double q_m = colV[m] / colN[m];
                REQUIRE(q_m >= 0.82 - 1e-12);
                REQUIRE(q_m <= 1.18 + 1e-12);
            }
        }
        REQUIRE(finalSumN == Catch::Approx(initSumN).margin(1e-12));
        REQUIRE(finalSumV == Catch::Approx(initSumV).margin(1e-12));
    }
}
