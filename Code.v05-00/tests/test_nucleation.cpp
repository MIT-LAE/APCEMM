#include "AIM/Nucleation.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

using namespace AIM;

TEST_CASE ("Nucleation", "[single-file]" ) {

    SECTION( "Surface tension" ) {
        // All tests based on data from Figure 1, Vehkamaki et al. (2002)

        // x_m = 0.123
        REQUIRE( sigma(0.123, 192.767) ==  Catch::Approx(0.0866).margin(1e-4));
        REQUIRE( sigma(0.123, 301.697) ==  Catch::Approx(0.0722).margin(1e-4));

        // x_m = 0.503
        REQUIRE( sigma(0.503, 191.882) ==  Catch::Approx(0.0833).margin(1e-4));
        REQUIRE( sigma(0.503, 302.583) ==  Catch::Approx(0.0757).margin(1e-4));

        // x_m = 0.672
        REQUIRE( sigma(0.672, 190.774) ==  Catch::Approx(0.0789).margin(1e-4));
        REQUIRE( sigma(0.672, 303.469) ==  Catch::Approx(0.0745).margin(1e-4));
    }

    SECTION( "Density sulfuric acid solution" ) {
        // All tests based on data from Figure 2, Vehkamaki et al. (2002)
        // x_m = 0.123
        REQUIRE( rho(0.123, 192.289)/1000. ==  Catch::Approx(1.112).epsilon(1e-2) );
        REQUIRE( rho(0.123, 304.971)/1000. ==  Catch::Approx(1.0776).epsilon(1e-2));

        // x_m = 0.503
        REQUIRE( rho(0.503, 191.214)/1000. ==  Catch::Approx(1.482).epsilon(1e-2));
        REQUIRE( rho(0.503, 304.046)/1000. ==  Catch::Approx(1.389).epsilon(1e-2));

        // x_m = 0.765
        REQUIRE( rho(0.765, 190.520)/1000. ==  Catch::Approx(1.798).epsilon(1e-2));
        REQUIRE( rho(0.765, 303.584)/1000. ==  Catch::Approx(1.682).epsilon(1e-2));
    }

    SECTION( "Critical cluster mole fraction, sulfuric acid") {
        
        // When RH=1, sulfur_conc = 1, expression simplifies greatly
        REQUIRE(x_star(250.0, 1.0, 1.0) == Catch::Approx(0.740997 - 0.00266379*250.0));

        // Regression test
        REQUIRE(x_star(250.0, 0.5, 1.0e10) == Catch::Approx(0.3136351646));
    }

    SECTION( "H2SO4 Nucleation rate") {
        // Data obtained from Figure 2 in Maattanen et al. (2017)
        
        // [H2SO4] = 5e5 /cm3
        double x_m = x_star(200.0, 0.01e-2, 5.0e5);
        REQUIRE(nuclRate(200.0, x_m, 0.01e-2, 5.0e5) == Catch::Approx(0.000007315869129751891).epsilon(0.2));
        x_m = x_star(191.4772, 0.01e-2, 5.0e5);
        REQUIRE(nuclRate(191.4772, x_m, 0.01e-2, 5.0e5) == Catch::Approx(0.026366508987303607).epsilon(0.2));

        // [H2SO4] = 5e8 /cm3, 
        x_m = x_star(241.12, 1.0, 5.0e8);

        // This test acts more like order of magnitude test
        // (I suspect the graphical data extraction method used is very sensitive
        //  for a log plot like the one used)
        REQUIRE(nuclRate(241.12, x_m, 1.0, 5.0e8) == Catch::Approx(31084029.570249166).epsilon(0.5));
        
        x_m = x_star(288.94, 1.0, 5.0e8);
        REQUIRE(nuclRate(288.94, x_m, 1.0, 5.0e8) == Catch::Approx(0.04859226736697418).epsilon(0.2));


        // This test is based on data from Figure 8 in Vehkamaki et al. (2002)
        x_m = x_star(236.0, 0.55, 666084629.0809168);

        REQUIRE(nuclRate(236.0, x_m, 0.55, 666084629.0809168) == Catch::Approx( 142197696.75857794).epsilon(1e-1));
    }

    SECTION ( "Critical cluster molecule number" ) {
        // This test is based on data from Figure 8 in Vehkamaki et al. (2002)
        // Some deviation is expected as the plot shows theoretical cluster size
        double x_m = x_star(236.0, 0.55, 666084629.0809168);
        REQUIRE(nTot(236.,x_m, 0.55, 4795689.579) == Catch::Approx(10.25).epsilon(1e-1));
    }

    SECTION ( "Radius of cluster size" ) {
        // Regression test
        REQUIRE(radCluster(0.3, 1.0e9) == Catch::Approx(223.6e-9).epsilon(1e-4));
    }

    SECTION ( "Treshold sulfuric acid concentration" ) {
        // Based on Figure 7 in Vehkamaki et al. (2002)
        REQUIRE(nThresh(190.15, 13.977e-2) == Catch::Approx(67728.15277441413).epsilon(0.1));
        REQUIRE(nThresh(240.15, 97.14e-2) == Catch::Approx(2894266.1247167457).epsilon(0.1));
    }
}