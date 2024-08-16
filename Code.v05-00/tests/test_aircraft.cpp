#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <string>
#include <Core/Aircraft.hpp>
#include <Util/ForwardDecl.hpp>
#include "Util/PhysConstant.hpp"


std::string engineFileName = std::string(APCEMM_TESTS_DIR)+"/../../input_data/ENG_EI.txt";
TEST_CASE("Vortex Losses Survival Fraction"){
    /* Comparing against results from Unterstrasser (2016)
    
     Issue is that params like z_emit, N_BV, and beta (the parameter chosen in calculating z_desc) are arbitrarily chosen and hard coded. 
     Also, it's assumed that all soot particles become ice.
     How to quantify this uncertainty? Therefore tests are just order of magnitude sanity checks.
     Function works perfectly when the aforementioned parameters are manually prescribed. -Michael
     */
   SECTION("Unterstrasser Table A2 #25"){
          //B747, EI_iceno = 2.8e14, 217 K, RH = 120%, N_BV = 0.015 s-1, z_atm = 148m, z_desc = 361m, z_emit = 98m
          Aircraft aircraft("B747", engineFileName, 200000.0, 217.0, 22000.0, 148.0, 0.015);
          double EI_ice_target = 2.8E14;
          double EISootRad = 20.0E-9;
          double volume = 4.0/3.0 * physConst::PI * pow(EISootRad, 3.0);
          double mass = physConst::RHO_SOOT*volume*1.0E3; // convert to grams
          double EI_soot = EI_ice_target * mass; 
          double z_atm = 148;

          REQUIRE(aircraft.VortexLosses(EI_soot, EISootRad, z_atm) == Catch::Approx(0.506).margin(0.2));
   }
      SECTION("Unterstrasser Table A2 #80"){
          //B747, EI_iceno = 1.22e15, 217 K, RH = 120%, N_BV = 0.015 s-1, z_atm = 148m, z_desc = 361m, z_emit = 98m
          Aircraft aircraft("B747", engineFileName, 200000.0, 217.0, 22000.0, 148.0, 0.013);
          double EI_ice_target = 1.22E15;
          double EISootRad = 20.0E-9;
          double volume = 4.0/3.0 * physConst::PI * pow(EISootRad, 3.0);
          double mass = physConst::RHO_SOOT*volume*1.0E3; //convert to grams
          double EI_soot = EI_ice_target * mass; 
          double z_atm = 148;

          REQUIRE(aircraft.VortexLosses(EI_soot, EISootRad, z_atm) == Catch::Approx(0.218).margin(0.2));
   }

}