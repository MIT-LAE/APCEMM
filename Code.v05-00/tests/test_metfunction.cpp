#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <Util/PhysFunction.hpp>
#include <Util/PhysConstant.hpp>
#include <Util/MetFunction.hpp>

using namespace met;

TEST_CASE( "MetFunction", "[single-file]" ) {
	SECTION("Compute Lapse Rate"){
		//no idea where this equation comes from.

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
	SECTION("Linear Interpolation vectors"){
		Vector_1D x = {0.0, 1.0, 2.0, 4.0, 5.0};
		Vector_1D x2= {5.0, 4.0, 2.0, 1.0, 0.0};
		Vector_1D y = {0.0, 3.0, 4.0, -10.0, 1.0};
		double xq1 = 0.5;
		double xq2 = 1.2;
		double xq3 = 3.0;
		double xq4 = 0.0;
		double xq5 = 5.0;
		REQUIRE(linearInterp(x, y, xq1) == Catch::Approx(1.5));
		REQUIRE(linearInterp(x, y, xq2) == Catch::Approx(3.2));
		REQUIRE(linearInterp(x, y, xq3) == Catch::Approx(-3.0));
		REQUIRE(linearInterp(x, y, xq4) == Catch::Approx(0.0));
		REQUIRE(linearInterp(x, y, xq5) == Catch::Approx(1.0));

		REQUIRE(linearInterp(x2, y, xq3) == Catch::Approx(3.5));
		REQUIRE(linearInterp(x2, y, xq4) == Catch::Approx(1.0));
		REQUIRE(linearInterp(x2, y, xq5) == Catch::Approx(0.0));

		REQUIRE_THROWS_AS(linearInterp(x, y, -1), std::range_error);
		REQUIRE_THROWS_AS(linearInterp(x, y, 6), std::range_error);
		REQUIRE_THROWS_AS(linearInterp(x2, y, -1), std::range_error);
		REQUIRE_THROWS_AS(linearInterp(x2, y, 6), std::range_error);
	}
	SECTION("Linear Interpolation 2 points"){
		REQUIRE(linearInterp(0, 0, 1, 3, 0.333) == Catch::Approx(0.999));
		REQUIRE(linearInterp(1, 1, 2, 2, 1) == Catch::Approx(1));
		REQUIRE(linearInterp(1, 1, 2, 2, 2) == Catch::Approx(2));
		REQUIRE(linearInterp(-1, -1, -2, -2, -1.5) == Catch::Approx(-1.5));

		REQUIRE_THROWS_AS(linearInterp(1, 1, 3, 3, 0), std::range_error);
		REQUIRE_THROWS_AS(linearInterp(1, 2, 2, 3, 4), std::range_error);
		REQUIRE_THROWS_AS(linearInterp(-1, 2, -2, 3, 1), std::range_error);

		REQUIRE(linearInterp(0, 2, 1, 4, 2, true) == 6);
		REQUIRE(linearInterp(0, 2, -1, 0, -3, true) == -4);

	}
	SECTION("Interp Met Data"){
		Vector_1D alts({0,1,2,3,4});
		Vector_1D metData({0,2,4,6,8});

		Vector_1D alts2= {0, -1, -2, -3, -4};

		REQUIRE(linInterpMetData(alts, metData, 2.2) == Catch::Approx(4.4));
		REQUIRE(linInterpMetData(alts2, metData, -3.3) == Catch::Approx(6.6));
		REQUIRE(linInterpMetData(alts, metData, 2) == 4);
		REQUIRE(linInterpMetData(alts, metData, 3) == 6);
		REQUIRE(linInterpMetData(alts, metData, 0) == 0);
		REQUIRE(linInterpMetData(alts, metData, 4) == 8);
		REQUIRE(linInterpMetData(alts2, metData, -2) == 4);
		REQUIRE(linInterpMetData(alts2, metData, 0) == 0);
		REQUIRE(linInterpMetData(alts2, metData, -4) == 8);
		REQUIRE_THROWS_AS(linInterpMetData(alts, metData, 10), std::range_error);

	}
	SECTION("Saturation Depth Calculation"){
		Vector_1D RHw = {50, 60, 70, 80, 90, 100, 110};
		Vector_1D RHw2 = {1, 2, 3, 4, 5, 6, 7};
		Vector_1D RHw3 = {150, 160, 170, 180, 190, 200, 210};
		Vector_1D T = {240, 235, 230, 225, 220, 215, 210};
		Vector_1D alt = {500, 600, 700, 800, 900, 1000, 1100};
		Vector_1D RHi(7);
		Vector_1D RHi2(7);
		Vector_1D RHi3(7);
		for (int i = 0; i < 7; i++){
			RHi[i] = RHw[i]*physFunc::pSat_H2Ol(T[i])/physFunc::pSat_H2Os(T[i]);
			RHi2[i] = RHw2[i]*physFunc::pSat_H2Ol(T[i])/physFunc::pSat_H2Os(T[i]);
			RHi3[i] = RHw3[i]*physFunc::pSat_H2Ol(T[i])/physFunc::pSat_H2Os(T[i]);
		}
		//RHi = 69.2755 87.347 107.062 128.525 151.839 177.105 204.417
		REQUIRE(satdepth_calc(RHi, alt, 6, 1500) == Catch::Approx(435.82044));

		//Confirm return value of 1 when current altitude is not Ice Supersaturated
		REQUIRE(satdepth_calc(RHi2, alt, 6, 1500) == Catch::Approx(1.0));


		//Warnings for sat depth exceeding domain
		try{
			satdepth_calc(RHi3, alt, 6, 1500);
		}
		catch(std::out_of_range &e){
			std::string error = e.what();
			REQUIRE(error == "In met::satdepth_calc: No end of ice supersaturated layer found.");
		}
	}
	SECTION("Calculating dx and x coordinates on expanding/shrinking dy") {
		int nx = 5;
		Vector_1D dy_old = {3.0, 5.0};
		Vector_1D dy_new = {5.0, 3.0};
		Vector_1D dx_old = {3.0, 5.0};
		Vector_1D x0_old = {-7.5, -12.5};
		auto pair = met::calcNewXCoords(dy_old, dy_new, x0_old, dx_old, nx);
		REQUIRE(pair.dx_new == Vector_1D{5.0, 3.0});
		REQUIRE(pair.x0_new == Vector_1D{-12.5, -7.5});

		nx = 4;
		x0_old = {-6, -10};
		pair = met::calcNewXCoords(dy_old, dy_new, x0_old, dx_old, nx);
		REQUIRE(pair.dx_new == Vector_1D{5.0, 3.0});
		REQUIRE(pair.x0_new == Vector_1D{-10, -6});

	}
}
