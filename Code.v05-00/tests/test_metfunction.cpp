#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <Util/PhysFunction.hpp>
#include <Util/PhysConstant.hpp>
#include <Util/MetFunction.hpp>
#include <iostream>

using namespace met;
using std::cout;
using std::endl;
TEST_CASE( "MetFunction", "[single-file]" ) {
	SECTION("Compute Lapse Rate"){
		//no idea where this equation comes from.

	}
	SECTION("Nearest Neighbor Float Array, Increasing"){
		float x1[] = {0.0, 1.0, 2.0, 3.0, 5.0};

		REQUIRE(nearestNeighbor(x1, -100000,5) == 0);
		REQUIRE(nearestNeighbor(x1, 0.5, 5) == 0);
		REQUIRE(nearestNeighbor(x1, 1.0, 5) == 1);
		REQUIRE(nearestNeighbor(x1, 2.0, 5) == 2);
		REQUIRE(nearestNeighbor(x1, 3.0, 5) == 3);
		REQUIRE(nearestNeighbor(x1, 5.0, 5) == 4);
		REQUIRE(nearestNeighbor(x1, 100000, 5) == 4);
	}

	SECTION("Nearest Neighbor Float Array, Decreasing"){
		float x1[] = {5.0, 3.0, 2.0, 1.0, 0.0};

		REQUIRE(nearestNeighbor(x1, -100000,5) == 4);
		REQUIRE(nearestNeighbor(x1, 0.5, 5) == 3);
		REQUIRE(nearestNeighbor(x1, 1.0, 5) == 3);
		REQUIRE(nearestNeighbor(x1, 2.0, 5) == 2);
		REQUIRE(nearestNeighbor(x1, 3.0, 5) == 1);
		REQUIRE(nearestNeighbor(x1, 5.0, 5) == 0);
		REQUIRE(nearestNeighbor(x1, 100000, 5) == 0);
	}
	SECTION("Nearest Neighbor Vector"){
		Vector_1D v = {0.0, 1.0, 3.0};
		float xq11 = 1.1;
		float xq12 = 0.1;
		float xq13 = 2.9;
		float xq14 = 0.0;
		float xq15 = 3.0;
		REQUIRE(nearestNeighbor(v,xq11) == 1);
		REQUIRE(nearestNeighbor(v,xq12) == 0);
		REQUIRE(nearestNeighbor(v,xq13) == 2);
		REQUIRE(nearestNeighbor(v,xq14) == 0);
		REQUIRE(nearestNeighbor(v,xq15) == 2);

	}
	SECTION("Linear Interpolation"){
		float x[] = {0.0, 1.0, 2.0, 4.0, 5.0};
		float x2[]= {5.0, 4.0, 2.0, 1.0, 0.0};
		float y[] = {0.0, 3.0, 4.0, -10.0, 1.0};
		float xq1 = 0.5;
		float xq2 = 1.2;
		float xq3 = 3.0;
		float xq4 = 0.0;
		float xq5 = 5.0;
		REQUIRE(linearInterp(x, y, xq1, 5) == Catch::Approx(1.5));
		REQUIRE(linearInterp(x, y, xq2, 5) == Catch::Approx(3.2));
		REQUIRE(linearInterp(x, y, xq3, 5) == Catch::Approx(-3.0));
		REQUIRE(linearInterp(x, y, xq4, 5) == Catch::Approx(0.0));
		REQUIRE(linearInterp(x, y, xq5, 5) == Catch::Approx(1.0));

		REQUIRE(linearInterp(x2, y, xq3, 5) == Catch::Approx(3.5));
		REQUIRE(linearInterp(x2, y, xq4, 5) == Catch::Approx(1.0));
		REQUIRE(linearInterp(x2, y, xq5, 5) == Catch::Approx(0.0));

		REQUIRE_THROWS_AS(linearInterp(x, y, -1, 5), std::range_error);
		REQUIRE_THROWS_AS(linearInterp(x, y, 6, 5), std::range_error);
		REQUIRE_THROWS_AS(linearInterp(x2, y, -1, 5), std::range_error);
		REQUIRE_THROWS_AS(linearInterp(x2, y, 6, 5), std::range_error);
	}
	SECTION("Saturation Depth Calculation"){
		float RHw[] = {50, 60, 70, 80, 90, 100, 110};
		float RHw2[] = {1, 2, 3, 4, 5, 6, 7};
		float RHw3[] = {150, 160, 170, 180, 190, 200, 210};
		float T[] = {240, 235, 230, 225, 220, 215, 210};
		float alt[] = {500, 600, 700, 800, 900, 1000, 1100};
		float RHi[7];
		for (int i = 0; i < 7; i++){
			RHi[i] = RHw[i]*physFunc::pSat_H2Ol(T[i])/physFunc::pSat_H2Os(T[i]);
		}
		//RHi = 69.2755 87.347 107.062 128.525 151.839 177.105 204.417
		REQUIRE(satdepth_calc(RHw, T, alt, 6, 7) == Catch::Approx(435.82044));

		//Confirm return value of 1 when current altitude is not Ice Supersaturated
		REQUIRE(satdepth_calc(RHw2, T, alt, 6, 7) == Catch::Approx(1.0));


		//Warnings for sat depth exceeding domain
		try{
			satdepth_calc(RHw3, T, alt, 6, 7);
		}
		catch(std::out_of_range e){
			std::string error = e.what();
			REQUIRE(error == "In met::satdepth_calc: No end of ice supersaturated layer found.");
		}
	}
}
