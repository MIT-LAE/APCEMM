#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <SANDS/Solver.hpp>
#include "test_solver_helpfunctions.hpp"
//default all transport values to 0
//test on 10k points
int nx = 101;
int ny = 101;
double xlim_left = 0.0;
double xlim_right = 1.0;
double ylim_down = 0.0;
double ylim_up = 1.0;
double dt = 0.1;
double Dh = 1;
double Dv = 1;

double dx = (xlim_right + xlim_left) / (nx - 1);
double dy = (ylim_up + ylim_down) / (ny - 1);

SANDS::Solver solver = SANDS::Solver(nx, ny, xlim_left, xlim_right,
                                     ylim_up, ylim_down); 
TEST_CASE("Test SANDS, low scale/timestep"){
    solver.UpdateTimeStep(dt);
    solver.Initialize(0, true, "./");
    Vector_2D phi = Vector_2D(ny ,Vector_1D(nx, 0.00));
    Vector_2D exact = phi;
    Vector_2D cellAreas = Vector_2D(ny ,Vector_1D(nx, dx*dy));
    Vector_1D y = Vector_1D(ny);
    //init y coords
    for(int i = 0; i < ny; i++){
        y[i] = i * dy;
    }

    setVecToPrescribedSoln(phi, 0, dx, dy);
    setVecToPrescribedSoln(exact, 0, dx, dy);
    REQUIRE(l2norm(phi) == l2norm(exact)); //make sure init is fine
    SECTION("Do Nothing (No adv or diff)"){
        solver.UpdateTimeStep(dt);
        solver.UpdateAdv(0, 0);
        solver.UpdateDiff(0, 0);
        solver.UpdateShear(0, y);

        //Make sure it's stable when doing nothing.
        for(int i = 0; i < 10; i++){   
            solver.Run(phi, cellAreas, 1); 
        }
        REQUIRE(l2norm(phi)  == Catch::Approx(l2norm(exact)).epsilon(1e-10));
    }
    
    SECTION("Constant Velocity Advection"){
        dt = 0.01;
        solver.UpdateTimeStep(dt);
        solver.UpdateAdv(0.5, 0.5, false);
        solver.UpdateDiff(0, 0);
        solver.UpdateShear(0, y);

        initTestVecAdvection(phi, dx, dy);
        double max, x_max, y_max;
        std::tie(max, x_max, y_max) = vecMax2D(phi, dx, dy);
        double x_max_init = x_max;
        double y_max_init = y_max;

        REQUIRE(max == Catch::Approx(1.0).epsilon(0.01));
        for(int i = 0; i < 100; i++){   
            solver.Run(phi, cellAreas); 
        }
        //Ensure no artificial diffusion on this scale
        std::tie(max, x_max, y_max) = vecMax2D(phi, dx, dy);
        REQUIRE(max == Catch::Approx(1.0).epsilon(0.02));
        REQUIRE(y_max == Catch::Approx(y_max_init + 0.5).epsilon(0.02));
        REQUIRE(x_max == Catch::Approx(x_max_init + 0.5).epsilon(0.02));

    }
    SECTION("Shear"){
        dt = 0.01;
        solver.UpdateTimeStep(dt);
        solver.UpdateAdv(0, 0);
        solver.UpdateDiff(0, 0);
        solver.UpdateShear(0.5, y, false); 
        //u_top should be 0.25, u_bot should be -0.25

        initTestVecShear(phi, dx, dy);
        double max, x_max, y_max;
        std::tie(max, x_max, y_max) = vecMax2D(phi, dx, dy);
        double x_max_init = x_max;
        double y_max_init = y_max;

        REQUIRE(max == Catch::Approx(1.0).epsilon(0.01));
        REQUIRE(x_max == Catch::Approx(0.125).epsilon(0.05));
        for(int i = 0; i < 100; i++){   
            solver.Run(phi, cellAreas); 
        }
        //top should be at 0.125, bot at 0.625
        int ind_x_top = static_cast<int>(0.125/dx); 
        int ind_x_bot = static_cast<int>(0.625/dx);

        double maxbot = 0, maxtop = 0;
        double x_maxbot, x_maxtop;
        for(int i = 0; i < nx; i++){
            if(phi[ny-1][i] > maxbot){
                maxbot = phi[ny-1][i];
                x_maxbot = dx*i;
            }
            if(phi[0][i] > maxtop){
                maxtop = phi[0][i];
                x_maxtop = dx*i;
            }
        }
        REQUIRE(maxbot == Catch::Approx(1.0).epsilon(0.01));
        REQUIRE(maxtop == Catch::Approx(1.0).epsilon(0.01));
        REQUIRE(x_maxbot == Catch::Approx(0.625).epsilon(0.05));
        REQUIRE(x_maxtop == Catch::Approx(0.125).epsilon(0.05)); 
    }
}