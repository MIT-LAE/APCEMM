#include "AIM/buildKernel.hpp"
#include "Util/ForwardDecl.hpp"
#include "Util/PhysConstant.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

using namespace AIM;


TEST_CASE("Build kernel", "[single-file]") {
    
    SECTION("Brownian kernel 1D") {
        // Tests from Seinfeld and Pandis (2006) Table 13.3
        // Act more like order of magnitude tests

        double T = 288.15;
        double p = 101325.0;
        Vector_1D bin_centers(1);
        bin_centers[0] = 0.5*0.002e-6;
        double rho_1 =  1.0e3;
        double rho_2 = 1.0e3;
        double bin_R = 0.5*0.002e-6;

        Vector_1D kernel = buildBrownianKernel(T, p, bin_centers, rho_1, bin_R, rho_2);

        REQUIRE(1e6 * kernel[0] == Catch::Approx(8.9e-10).epsilon(0.5));

        // Test from Seinfeld and Pandis (2006) Table 13.3
        
        bin_R = 0.5*10.0e-6;

        kernel = buildBrownianKernel(T, p, bin_centers, rho_1, bin_R, rho_2);

        REQUIRE(1e6 * kernel[0] == Catch::Approx(8.5e-5).epsilon(1e-1));

    }

    SECTION("Brownian kernel 2D") {
        double T = 288.15;
        double p = 101325.0;
        Vector_1D bin_centers1(1);
        Vector_1D bin_centers2(1);
        bin_centers1[0] = 0.5*0.002e-6;
        bin_centers2[0] = 0.5*0.002e-6;
        double rho_1 =  1.0e3;
        double rho_2 = 1.0e3;
        // double bin_R = 0.5*0.002e-6;

        Vector_2D kernel = buildBrownianKernel(T, p, bin_centers1, rho_1, bin_centers2, rho_2);

        REQUIRE(1e6 * kernel[0][0] == Catch::Approx(8.9e-10).epsilon(0.5));

    }

    SECTION("Convective diffusion enhancement 1D") {
        // Based on figure 15.7 in Jacobson (2005)

        double T = 298.0;
        double p = 100000;
        Vector_1D bin_centers(1);
        bin_centers[0] = 0.01e-6;
        double rho_1 =  1.0e3;
        double rho_2 = 1.0e3;
        double bin_R = 0.08642108575959868e-6;

        Vector_1D brownian_kernel = buildBrownianKernel(T, p, bin_centers, rho_1, bin_R, rho_2);
        Vector_1D DEkernel = buildDEKernel(T, p, bin_centers, rho_1, bin_R, rho_2, brownian_kernel);
       
        REQUIRE(1e6 * DEkernel[0] == Catch::Approx(2.2060e-10).epsilon(0.2));

        bin_R = 1.102182525205009e-6;
        brownian_kernel = buildBrownianKernel(T, p, bin_centers, rho_1, bin_R, rho_2);
        DEkernel = buildDEKernel(T, p, bin_centers, rho_1, bin_R, rho_2, brownian_kernel);
       
        REQUIRE(1e6 * DEkernel[0] == Catch::Approx(3.0519786963413166e-8).epsilon(0.2));

    }

    SECTION("Convective diffusion enhancement 2D") {
        // Based on figure 15.7 in Jacobson (2005)

        double T = 298.0;
        double p = 100000;
        Vector_1D bin_centers1(1);
        Vector_1D bin_centers2(1);
        bin_centers1[0] = 0.01e-6;
        bin_centers2[0] = 0.08642108575959868e-6;
        double rho_1 =  1.0e3;
        double rho_2 = 1.0e3;

        Vector_2D brownian_kernel = buildBrownianKernel(T, p, bin_centers1, rho_1, bin_centers2, rho_2);
        Vector_2D DEkernel = buildDEKernel(T, p, bin_centers1, rho_1, bin_centers2, rho_2, brownian_kernel);

        REQUIRE(1e6 * DEkernel[0][0] == Catch::Approx(2.2060e-10).epsilon(0.2));
    }

    SECTION("Gravitational collection kernel 1D"){
        // Based on figure 15.7 in Jacobson (2005)

        double T = 298.0;
        double p = 100000;
        Vector_1D bin_centers(1);
        bin_centers[0] = 10.0e-6;
        double rho_1 =  1.0e3;
        double rho_2 = 1.0e3;
        double bin_R = 4.098982426133141e-6;

        Vector_1D GCKernel = buildGCKernel(T, p, bin_centers, rho_1, bin_R, rho_2);

        //REQUIRE(1.0e6 * GCKernel[0] == Catch::Approx(2.948256178847939e-7).epsilon(0.5));

    }

    SECTION("Gravitational collection kernel 2D"){
        // Based on figure 15.7 in Jacobson (2005)

        double T = 298.0;
        double p = 100000;
        Vector_1D bin_centers1(1);
        bin_centers1[0] = 10.0e-6;
        Vector_1D bin_centers2(1);
        bin_centers2[0] = 4.098982426133141e-6;
        double rho_1 =  1.0e3;
        double rho_2 = 1.0e3;

        Vector_2D GCKernel = buildGCKernel(T, p, bin_centers1, rho_1, bin_centers2, rho_2);

        //REQUIRE(1.0e6 * GCKernel[0][0] == Catch::Approx(2.948256178847939e-7).epsilon(0.5));

    }

    SECTION("Turbulent inertial motion kernel 1D") {
        // Based on Figure 15.7 from Jacobson (2005)
        
        double T = 298.0;
        double p = 100000;
        Vector_1D bin_centers(1);
        bin_centers[0] = 10.0e-6;
        double rho_1 =  1.0e3;
        double rho_2 = 1.0e3;
        double bin_R = 0.04303311824915588e-6;
        
        // Use this to account for the different epsilon value used
        double scaler = pow(5.0e-4, 0.75)/pow(physConst::EPSILON, 0.75);
        Vector_1D TIKernel = buildTIKernel(T, p, bin_centers, rho_1, bin_R, rho_2);

        REQUIRE(1.0e6 * scaler * TIKernel[0] == Catch::Approx(2.0976794828100117e-8).epsilon(0.1));
        
        bin_R = 4.591682295070266e-6;
        TIKernel = buildTIKernel(T, p, bin_centers, rho_1, bin_R, rho_2);
        REQUIRE(1.0e6 * scaler * TIKernel[0] == Catch::Approx(3.6018177927383583e-8).epsilon(0.1));

    }

    SECTION("Turbulent inertial kernel 2D"){
        // Based on figure 15.7 in Jacobson (2005)

        double T = 298.0;
        double p = 100000;
        Vector_1D bin_centers1(1);
        bin_centers1[0] = 10.0e-6;
        Vector_1D bin_centers2(1);
        bin_centers2[0] = 0.04303311824915588e-6;
        double rho_1 =  1.0e3;
        double rho_2 = 1.0e3;

       // Use this to account for the different epsilon value used
        double scaler = pow(5.0e-4, 0.75)/pow(physConst::EPSILON, 0.75);
        Vector_2D TIKernel = buildTIKernel(T, p, bin_centers1, rho_1, bin_centers2, rho_2);

        REQUIRE(1.0e6 * scaler * TIKernel[0][0] == Catch::Approx(2.0976794828100117e-8).epsilon(0.1));
        
        bin_centers2[0] = 4.591682295070266e-6;
        TIKernel = buildTIKernel(T, p, bin_centers1, rho_1, bin_centers2, rho_2);
        REQUIRE(1.0e6 * scaler * TIKernel[0][0] == Catch::Approx(3.6018177927383583e-8).epsilon(0.1));
    }

    SECTION("Turbulent shear kernel 1D") {
        // Based on figure 15.7 in Jacobson (2005)
        
        double T = 298.0;
        double p = 100000;
        Vector_1D bin_centers(1);
        bin_centers[0] = 10.0e-6;
        double rho_1 =  1.0e3;
        double rho_2 = 1.0e3;
        double bin_R = 0.020082459574328918e-6;
        
        // Use this to account for the different epsilon value used
        double scaler = pow(5.0e-4, 0.5)/pow(physConst::EPSILON, 0.5);
        Vector_1D TSKernel = buildTSKernel(T, p, bin_centers, rho_1, bin_R, rho_2);

        REQUIRE(1.0e6 * scaler * TSKernel[0] == Catch::Approx(6.700187503509577e-9).epsilon(0.1));
        
        bin_R = 5.227690396963787e-6;
        TSKernel = buildTSKernel(T, p, bin_centers, rho_1, bin_R, rho_2);
        REQUIRE(1.0e6 * scaler * TSKernel[0] == Catch::Approx(2.5118864315095718e-8).epsilon(0.1));


    }

    SECTION("Turbulent shear kernel 2D"){
        // Based on figure 15.7 in Jacobson (2005)

        double T = 298.0;
        double p = 100000;
        Vector_1D bin_centers1(1);
        bin_centers1[0] = 10.0e-6;
        Vector_1D bin_centers2(1);
        bin_centers2[0] = 0.020082459574328918e-6;
        double rho_1 =  1.0e3;
        double rho_2 = 1.0e3;

       // Use this to account for the different epsilon value used
        double scaler = pow(5.0e-4, 0.5)/pow(physConst::EPSILON, 0.5);
        Vector_2D TSKernel = buildTSKernel(T, p, bin_centers1, rho_1, bin_centers2, rho_2);

        REQUIRE(1.0e6 * scaler * TSKernel[0][0] == Catch::Approx(6.700187503509577e-9).epsilon(0.1));
        
        bin_centers2[0] = 5.227690396963787e-6;
        TSKernel = buildTSKernel(T, p, bin_centers1, rho_1, bin_centers2, rho_2);
        REQUIRE(1.0e6 * scaler * TSKernel[0][0] == Catch::Approx(2.5118864315095718e-8).epsilon(0.1));
    }

}