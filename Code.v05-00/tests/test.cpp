#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <Util/PhysFunction.hpp>
#include <Util/PhysConstant.hpp>

using namespace physFunc;

TEST_CASE( "Saturation pressures water", "[single-file]" ) {

    SECTION( "Triple point" ) {
        // Ensure triple point is on saturation curves
        REQUIRE( pSat_H2Ol(273.15) ==  Catch::Approx(611.73).epsilon(0.1));
        REQUIRE( pSat_H2Os(273.15) ==  Catch::Approx(611.73).epsilon(0.1));
    }
    // TODO: Add tests based on Sonntag 1994 paper

    // Check derivatives
    SECTION( "Derivatives" ) {
        REQUIRE((pSat_H2Os(273.15+0.01) - pSat_H2Os(273.15-0.01))/0.02 == Catch::Approx(dpSat_H2Os(273.15)));
    }
}

TEST_CASE( "Saturation pressures acids", "[single-file]") {

    SECTION( "Sulfuric acid" ) {
        // On substitution of T=360.15, the expression simplifies to p_atm*exp(-L)
        REQUIRE(pSat_H2SO4(360.15) == Catch::Approx(101325. * exp(-11.94) ).epsilon(0.01));
    }

    SECTION( "Nitric acid" ) {
        // T = 190, point extracted from Figure 1 in Hanson and Mauersberg(1988)
        REQUIRE(pSat_HNO3(190, 0.00001655285345714208 * 133.322) == \
                Catch::Approx(0.000007280218604164274 * 133.322).epsilon(0.1));
        
        // T = 205, point extracted from Figure 1 in Hanson and Mauersberg(1988)
        REQUIRE(pSat_HNO3(205, 0.001098067646521057 * 133.322) == \
                Catch::Approx(0.0000010516756971907947 * 133.322).epsilon(0.1));
    }
}

TEST_CASE( "Air properties", "[single-file]") {

    SECTION( "Density" ) {
        // T = 288.15 K, p = 101325 Pa, rho = 1.225 kg/m3
        REQUIRE(rhoAir(288.15, 101325.0) == Catch::Approx(1.225));
    }
    SECTION( "Dynamic viscosity" ) {
        // T = 273.15 K, p = 101325 Pa, mu = 1.7231e-5 Ns/m^2
        REQUIRE(dynVisc(273.15) ==  Catch::Approx(1.7231E-5).epsilon(0.01));
        // T = 223.15 K, p = 101325 Pa, mu = 1.474e-5 Ns/m^2
        REQUIRE(dynVisc(223.15) ==  Catch::Approx(1.474E-5).epsilon(0.05));
    }

    SECTION( "Kinematic viscosity" ) {
        // T = 273.15 K, p = 101325 Pa, nu = 1.7231e-5 m^2/s
        REQUIRE(kinVisc(273.15, 101325) ==  Catch::Approx(1.338e-5).epsilon(0.01));
        // T = 223.15 K, p = 101325 Pa, nu = 9.319e-6 m^2/s
        REQUIRE(kinVisc(223.15, 101325) ==  Catch::Approx(9.319e-6).epsilon(0.05));
        // T = 220.0 K, p = 0.2454 bar, nu = 3.706e-5 m^2/s
        REQUIRE(kinVisc(223.25, 26500.) ==  Catch::Approx(3.5251e-5).epsilon(0.05));
    }

    SECTION( "Thermal speed" ) {
        // At sea-level in the International Standard Atmosphere
        REQUIRE(thermalSpeed(288.15) == Catch::Approx(458.9446));
        // At 10 km altitude in the International Standard Atmosphere
        REQUIRE(thermalSpeed(223.25) == Catch::Approx(403.9696));
    }

    SECTION( "Mean free path") {
        // At sea-level in the International Standard Atmosphere
        REQUIRE(lambda(288.15, 101325.) == Catch::Approx(6.6332e-8).epsilon(0.1));
        // At 10 km altitude in the International Standard Atmosphere
        REQUIRE(lambda(223.2521, 26500.) == Catch::Approx(1.9651e-7).epsilon(0.15));
    }
}

TEST_CASE("Particle dynamics", "[single-file]") {

    SECTION ("Mass of sphere") {
        REQUIRE(mass_sphere(pow(3/(4*physConst::PI), 1./3.), 1.0) == Catch::Approx(1.0)); 
    }

    SECTION ("Knudsen number") {
        // Example 20.1 from Jacobson (2005)
        REQUIRE(Kn(0.5e-6, 220., 25.0e2) == Catch::Approx(3.54).epsilon(0.1));
    }

    SECTION ("Fall speed") {
        // Example 20.1 from Jacobson (2005)
        REQUIRE(vFall(0.5e-6, 1.5e3, 220.0, 25.0e2) == Catch::Approx(0.037e-2).epsilon(0.1));
        REQUIRE(vFall(0.5e-6, 1.5e3, 220.0, 1013.0e2) == Catch::Approx(0.0063e-2).epsilon(0.1));
    }

    SECTION ("Slip flow correction") {
        // Example 20.1 from Jacobson (2005)
        REQUIRE(slip_flowCorrection(3.54) == Catch::Approx(6.58).epsilon(0.01));
    }

    SECTION ("Diffusion coefficient particle") {
        // Based off figure 9.8 in Seinfeld and Pandis (2006)
        REQUIRE(partDiffCoef(0.0005e-6, 293.15, 101325.) == Catch::Approx(5.0e-6).epsilon(0.1));
        REQUIRE(partDiffCoef(5e-6, 293.15, 101325.) == Catch::Approx(2.5e-12).margin(1e-13));
    }

    SECTION ("Particle mean free path") {

        // Based off table 9.5 in Seinfeld and Pandis (2006)

        // Determine particle mass
        double m_p = (8*physConst::kB * 293.15) / (pow(4965.0e-2,2) * physConst::PI);
        REQUIRE(lambda_p(0.001e-6, m_p, 293.15, 101325.) == Catch::Approx(6.59e-8).epsilon(0.1));

        m_p = (8*physConst::kB * 293.15) / (pow(4.96e-5,2) * physConst::PI);
        REQUIRE(lambda_p(10.0e-6, m_p, 293.15, 101325.) == Catch::Approx(6.08e-8).epsilon(0.1));
    }

    SECTION ("Delta_p") {
        // Regression test
        double m_p =  (8*physConst::kB * 293.15) / (pow(4.96e-5,2) * physConst::PI);
        REQUIRE(delta_p(1.0e-6, m_p, 293.15, 101325.) == Catch::Approx(3.909e-7).epsilon(0.01));
    }

    SECTION ("Reynolds number") {
        // Based off Figure 15.6 in Jacobson (2005)
        REQUIRE(1e10 * Reynolds_p(0.5 * 0.0197e-6, 1000.0, 292., 999.e2) == Catch::Approx(1.79).epsilon(0.01));
        REQUIRE(Reynolds_p(0.5 * 95.05469403923057e-6, 1000.0, 292., 999.e2) == Catch::Approx(1.707).epsilon(0.01));
    }

    SECTION( "Schmidt number" ) {
        // Regression test
        REQUIRE(Schmidt_p(1.0e-6, 288.15, 101325.0) == Catch::Approx(1151911.59));
    }

    SECTION( "Stokes number" ) {
        // Needs to be zero when particles have equal fall speeds
        REQUIRE(Stokes_p(1.0e-6, 1.0, 1.0e-6, 1.0, 298., 101325.) == Catch::Approx(0.0));
        
        // Needs to be equal to V_fall_max**2/(r_max*g) if one particle falls much
        // faster than the other
        double V = vFall(1e-8, 1000.0, 298.15, 101325.0);
        double St = pow(V,2)/((100.0e-6) * physConst::g);
        REQUIRE(Stokes_p(1.0e-8, 1000.0, 100.0e-6, 1e-8, 298.15, 101325.0) == Catch::Approx(St).epsilon(0.01));

        // Regression test
        REQUIRE(Stokes_p(1.0e-6, 1000.0, 1.0e-6, 10.0, 298.15, 101325.0) == Catch::Approx(1.65778e-5));
    }

    SECTION( "Aggregation efficiency" ){

        // Values are taken from Figure 2.5 from Ludlum (1980)

        // Reynolds number << 1
        double r_j = 20e-6;
        double r_i = 10e-6;
        double rho_j = 1000.0;
        double rho_i = 1000.0;
        double T = 293.15;
        double P = 101325.0;
        REQUIRE(E_agg(r_i, rho_i, r_j, rho_j, T, P) == Catch::Approx(0.25).epsilon(0.1));

        // Reynolds number > 1
        T = 100.0;
        P = 201325.0;
        REQUIRE(E_agg(r_i, rho_i, r_j, rho_j, T, P) == Catch::Approx(0.75).epsilon(0.1));

    }

    SECTION ( "Diffusion coefficients gases" ) {
        // HNO3 and H2SO4 tests based off Table 9.1 in Seinfeld and Pandis (2006) 
        double T = 298.;
        double P = 101325.;

        REQUIRE(DiffCoef_HNO3(T, P) == Catch::Approx(0.11e-4).epsilon(0.1));
        REQUIRE(DiffCoef_H2SO4(T, P) == Catch::Approx(0.10e-4).epsilon(0.1));

        REQUIRE(DiffCoef_H2O(273.15, P) == Catch::Approx(2.11e-5));

        // Based on Table 13.1 in Pruppacher and Klett (2010)
        REQUIRE(CorrDiffCoef_H2O(0.01e-6, 283.15, 800.0e2) == Catch::Approx(5.2e-8).epsilon(0.01));
        REQUIRE(CorrDiffCoef_H2O(10.0e-6, 283.15, 800.0e2) == Catch::Approx(0.19e-4).epsilon(0.05));

        // Regression tests for HNO3 and H2SO4 corrected diffusion coefficients
        REQUIRE(CorrDiffCoef_H2SO4(0.01e-6, 283.15, 800.0e2) == Catch::Approx(6.207e-7).epsilon(0.01));
        REQUIRE(CorrDiffCoef_HNO3(0.01e-6, 283.15, 800.0e2) == Catch::Approx(7.662e-7).epsilon(0.01));

    }

    SECTION ("Thermal Conductivity") {
        // Based on Table 13.1 in Pruppacher and Klett (2010)
        REQUIRE(ThermalCond(0.01e-6, 283.15, 800.0e2) == Catch::Approx(4.184 * 100. * 1.9e-6).epsilon(0.05));
        REQUIRE(ThermalCond(10.0e-6, 283.15, 800.0e2) == Catch::Approx(4.184 * 100. * 5.9e-5).epsilon(0.05));
    }
}

TEST_CASE("Microphysics", "[single-file]") {

    SECTION ("Latent heat of sublimation") {
        // Values from Table 2.1 from Rogers and Yau (1989)
        REQUIRE(LHeatSubl_H2O(273.15-40) == Catch::Approx(2839.0e3).epsilon(0.01));
        REQUIRE(LHeatSubl_H2O(273.15) == Catch::Approx(2834.0e3).epsilon(0.01));
        REQUIRE(LHeatSubl_H2O(273.15-20) == Catch::Approx(2838.0e3).epsilon(0.01));
    }

    SECTION ("Kelvin effect") {
        // Regression test
        REQUIRE(Kelvin(1.0e-9) == Catch::Approx(1.6487212707));
    }

    SECTION ("Growth rate") {
        double T = 273.15;
        double p_h2o = 1.01 * pSat_H2Ol(T);
        double p = 101325.0;
        double r = 0.75e-6;

        double H2O = 1e-6 * physConst::Na * p_h2o /(physConst::R * T);
        
        double volume_rate = growthRate(r, T, p, H2O);
        
        // Regression test
        REQUIRE( 1.e12 * volume_rate == Catch::Approx(1.3065441974));

        // TODO: add data-based test

    }
}