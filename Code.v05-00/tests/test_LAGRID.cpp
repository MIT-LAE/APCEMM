#include "LAGRID/RemappingFunctions.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <iostream>
TEST_CASE("FreeCoordBoxGrid and Remapping") {
    Vector_1D dy = {1, 2, 3, 4};
    Vector_1D dx = {1, 2, 3, 4};
    Vector_1D x0 = {7.5, 5, 2.5, 0};
    Vector_2D phi = {
        {0, 0, 2, 0, 0},
        {0, 1, 2, 1, 0},
        {0, 1, 1, 1, 0},
        {0, 0, 1, 0, 0}
    };

    LAGRID::FreeCoordBoxGrid testgrid(dx, dy, phi, x0, 0, 0.5);

    SECTION("FreeCoordBoxGrid Constructor") {
        //basic sanity checks
        REQUIRE(testgrid.boxes.size() == 8);
        REQUIRE(testgrid.boundaryBoxIndices.size() == 6);

        //boundary node order
        REQUIRE(testgrid.boundaryBoxIndices[0] == 0);
        REQUIRE(testgrid.boundaryBoxIndices[1] == 1);
        REQUIRE(testgrid.boundaryBoxIndices[2] == 3);
        REQUIRE(testgrid.boundaryBoxIndices[4] == 6);

        //location and property of boxes
        REQUIRE(testgrid.boxes[0].topLeftY == 1);
        REQUIRE(testgrid.boxes[0].topLeftX == 9.5);
        REQUIRE(testgrid.boxes[0].botRightY == 0);
        REQUIRE(testgrid.boxes[0].botRightX == 10.5);

        REQUIRE(testgrid.boxes[3].topLeftY == 3);
        REQUIRE(testgrid.boxes[3].topLeftX == 11);
        REQUIRE(testgrid.boxes[3].botRightY == 1);
        REQUIRE(testgrid.boxes[3].botRightX == 13);

        REQUIRE(testgrid.boxes[4].topLeftY == 6);
        REQUIRE(testgrid.boxes[4].topLeftX == 5.5);
        REQUIRE(testgrid.boxes[4].botRightY == 3);
        REQUIRE(testgrid.boxes[4].botRightX == 8.5);

        //area function and mass
        REQUIRE(testgrid.boxes[4].area() == 9);
        REQUIRE(testgrid.boxes[4].mass == 9);
    }

    LAGRID::Remapping testremap = { -1, -1, 2, 2, 11, 6 };
    SECTION("Diffusion Loss Function") {
        REQUIRE(LAGRID::diffusionLossFunctionExact(testgrid, testremap) == 6);

        //TODO: Test after faster algorithm is actually developed or deemed necessary through test runs
        //REQUIRE(LAGRID::diffusionLossFunctionBoundaryEstimate(testgrid, testremap) == 6);
    }

    auto twoDGrid = LAGRID::mapToStructuredGrid(testgrid, testremap);
    twoDGrid = LAGRID::mapToStructuredGrid(testgrid, testremap);
    double mass_before = std::accumulate(testgrid.boxes.begin(), testgrid.boxes.end(), 0.0, 
                                                [](double sum, auto& b) -> double {
                                                    return sum + b.mass;
                                                });
    SECTION("Remapping") {
        REQUIRE(twoDGrid.yCoords.size() == 6);
        REQUIRE(twoDGrid.yCoords[0] == 0);
        REQUIRE(twoDGrid.yCoords[3] == 6);
        REQUIRE(twoDGrid.yCoords[5] == 10);
        REQUIRE(twoDGrid.xCoords.size() == 11);
        REQUIRE(twoDGrid.xCoords[0] == 0);
        REQUIRE(twoDGrid.xCoords[3] == 6);
        REQUIRE(twoDGrid.xCoords[10] == 20);

        double mass_after = 0;
        for (std::size_t j = 0; j < twoDGrid.phi.size(); j++){
            for (std::size_t i = 0; i < twoDGrid.phi[0].size(); i++) {
                mass_after += twoDGrid.phi[j][i] * twoDGrid.dx * twoDGrid.dy;
            }
        }
        REQUIRE(std::abs(mass_before - mass_after) < 1e-12);

        // Uncomment / play around to verify that SFINAE for twoDGridVariable constructor works, this won't compile.
        // auto a = LAGRID::twoDGridVariable(Vector_2D(), 2, Vector_1D());
    }

    SECTION("Adding Buffer") {
        twoDGrid.addBuffer(6, 10, 4, 8, 0.0);
        REQUIRE(twoDGrid.yCoords.size() == 12);
        REQUIRE(twoDGrid.xCoords.size() == 19);
        REQUIRE(twoDGrid.xCoords[0] == -6);
        REQUIRE(twoDGrid.yCoords[0] == -8);

        REQUIRE(twoDGrid.yCoords[3] == -2);
        REQUIRE(twoDGrid.yCoords[4] == 0);
        REQUIRE(twoDGrid.yCoords[11] == 14);

        REQUIRE(twoDGrid.xCoords[2] == -2);
        REQUIRE(twoDGrid.xCoords[3] == 0);
        REQUIRE(twoDGrid.xCoords[18] == 30);

        REQUIRE(twoDGrid.phi.size() == 12);
        REQUIRE(twoDGrid.phi[0].size() == 19);
        double mass = 0;
        for (std::size_t j = 0; j < twoDGrid.phi.size(); j++){
            for (std::size_t i = 0; i < twoDGrid.phi[0].size(); i++) {
                mass += twoDGrid.phi[j][i] * twoDGrid.dx * twoDGrid.dy;
            }
        }
        REQUIRE(std::abs(mass_before - mass) < 1e-12);
    }
}
