#include <stdexcept>
#include <string>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "Core/Emission.hpp"
#include "Core/Engine.hpp"
#include "Core/Fuel.hpp"

namespace {

const std::string engineFileName = std::string(APCEMM_TESTS_DIR) + "/../defaults/ENG_EI.txt";

Engine testEngine() {
    /* Cruise-like conditions: T [K], p [Pa], RHw [%], Mach [-] */
    return Engine("GEnx-2B67B", engineFileName, 217.0, 22000.0, 40.0, 0.8);
}

/* Check every observable member, so that a copy that drops any member fails */
void requireSameEmission(const Emission &a, const Emission &b) {
    REQUIRE(a.getCO2() == b.getCO2());
    REQUIRE(a.getH2O() == b.getH2O());
    REQUIRE(a.getNOx() == b.getNOx());
    REQUIRE(a.getNO() == b.getNO());
    REQUIRE(a.getNO2() == b.getNO2());
    REQUIRE(a.getHNO2() == b.getHNO2());
    REQUIRE(a.getSO2() == b.getSO2());
    REQUIRE(a.getCO() == b.getCO());
    REQUIRE(a.getHC() == b.getHC());
    REQUIRE(a.getCH4() == b.getCH4());
    REQUIRE(a.getC2H6() == b.getC2H6());
    REQUIRE(a.getPRPE() == b.getPRPE());
    REQUIRE(a.getALK4() == b.getALK4());
    REQUIRE(a.getCH2O() == b.getCH2O());
    REQUIRE(a.getALD2() == b.getALD2());
    REQUIRE(a.getGLYX() == b.getGLYX());
    REQUIRE(a.getMGLY() == b.getMGLY());
    REQUIRE(a.getSO2toSO4() == b.getSO2toSO4());
    REQUIRE(a.getSoot() == b.getSoot());
    REQUIRE(a.getSootRad() == b.getSootRad());
    REQUIRE(a.getEngineName() == b.getEngineName());
    REQUIRE(a.getFuelChem() == b.getFuelChem());
}

} // namespace

TEST_CASE("Emission copy keeps every member", "[Emission]") {
    const Engine engine = testEngine();
    const Fuel fuel("C12H24", 600.0);
    const Emission source(engine, fuel, 0.02);

    REQUIRE(source.getSO2toSO4() == 0.02);
    REQUIRE(source.getSO2() > 0.0);

    SECTION("Copy construction") {
        const Emission copy(source);
        requireSameEmission(copy, source);
    }

    SECTION("Copy assignment") {
        /* Target built with a different ratio: an assignment that drops SO2toSO4 leaves 0.0 */
        Emission target(engine, fuel, 0.0);
        target = source;
        requireSameEmission(target, source);
    }
}

TEST_CASE("Fuel sulfur content", "[Fuel]") {
    SECTION("Zero FSC is kept and gives zero SO2") {
        const Fuel fuel("C12H24", 0.0);
        REQUIRE(fuel.getFSC() == 0.0);
        const Emission emission(testEngine(), fuel, 0.02);
        REQUIRE(emission.getSO2() == 0.0);
    }

    SECTION("Fuel::fromEISO2 gives back EI(SO2)") {
        const Fuel fuel = Fuel::fromEISO2("C12H24", 1.2);
        REQUIRE(fuel.EI_SO2() == Catch::Approx(1.2).epsilon(1e-15));
        const Emission emission(testEngine(), fuel, 0.02);
        REQUIRE(emission.getSO2() == fuel.EI_SO2());
    }

    SECTION("EI(SO2) per ppm of sulfur uses the SO2 to S molar mass ratio") {
        /* 1000 ppm = 1 g S/kg fuel, which gives about 2 g SO2/kg fuel */
        REQUIRE(Fuel("C12H24", 1000.0).EI_SO2() == Catch::Approx(MW_SO2 / MW_S).epsilon(1e-15));
    }

    SECTION("Negative FSC throws") {
        REQUIRE_THROWS_AS(Fuel("C12H24", -1.0), std::invalid_argument);
    }
}
