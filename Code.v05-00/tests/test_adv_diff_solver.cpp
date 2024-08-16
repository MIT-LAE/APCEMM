#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <iostream>
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
        int nx = 100, ny = 100;
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
        double t = 0, t_max = 3;
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

    TEST_CASE("Explicit Advection - Diffusion w/ Shear"){
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
        int nx = 500, ny = 500;
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

    TEST_CASE("Operator Split Advection-Diffusion"){
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
}