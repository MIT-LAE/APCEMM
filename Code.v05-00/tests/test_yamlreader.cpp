#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <fstream>
#include <Util/YamlUtils.hpp>
#include <YamlInputReader/YamlInputReader.hpp>
#include <YamlInputReader/YamlPathUtils.hpp>
#include <Core/Input.hpp>
#include "Core/OutputFilenames.hpp"
#include "APCEMM.h"

using namespace YamlInputReader;

string YAML_DIR = "/input-yamls";

//APCEMM_TEST_DIR is a preprocessor macro
TEST_CASE("YamlInputReader Helper Functions"){
    SECTION("Yaml Compiles"){
        string filename = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test.yaml";
        YAML::Node n = YAML::LoadFile(filename);
    }
    SECTION("String Trim"){
        REQUIRE(trim("") == "");
        REQUIRE(trim(" ") == "");
        REQUIRE(trim("    a b c") == "a b c");
        REQUIRE(trim(" a b c ") == "a b c");
        REQUIRE(trim("a b c    ") == "a b c");
    }
    SECTION("String Split"){
        vector<string> tokens;

        string teststr = "10:20:30";
        tokens = split(teststr, ":");
        REQUIRE(tokens.size() == 3);
        REQUIRE(tokens[0] == "10");
        REQUIRE(tokens[1] == "20");
        REQUIRE(tokens[2] == "30");

        REQUIRE(split(""," ").size() == 0);
        REQUIRE(split("abcde", " ")[0] == "abcde");

        //bad spacing test
        teststr = "10  20   30";
        tokens = split(teststr, " ");
        REQUIRE(tokens.size() == 3);
        REQUIRE(tokens[0] == "10");
        REQUIRE(tokens[1] == "20");
        REQUIRE(tokens[2] == "30");
        //zero delimiter width exception
        string err;
        try{
            tokens = split(" ", "");
        }
        catch(std::invalid_argument &e){
            err = string(e.what());
        }
        REQUIRE(err.find("In YamlInputReader::split: Delimiter length cannot be zero") == 0);
    }
    SECTION("String to double"){
            REQUIRE(strToDouble(" +20.0 ") == 20.0);
            REQUIRE(strToDouble(" 20. ") == 20.);
            REQUIRE(strToDouble("-30.0") == -30.0);
            REQUIRE(strToDouble(" 2e2 ") == 2e2);
            REQUIRE(strToDouble(" 2e+2 ") == 2e+2);
            REQUIRE(strToDouble(" 2e-4 ") == 2e-4);
            //exceptions
            string error = "";
            try{
                strToDouble("asdf");
            }
            catch (std::invalid_argument &e){
                error = e.what();
            }
            REQUIRE (error.find("Something went wrong with processing a number") == 0);
            error = "";
            try{
                strToDouble("10 2e3");
            }
            catch (std::invalid_argument &e){
                error = e.what();
            }
            REQUIRE (error.find("Something went wrong with processing a number") == 0);
            error = "";
            try{
                strToDouble("1e3e3");
            }
            catch (std::invalid_argument &e){
                error = e.what();
            }
            REQUIRE (error.find("Something went wrong with processing a number") == 0);
            error = "";
            try{
                strToDouble("2e-4.");
            }
            catch (std::invalid_argument &e){
                error = e.what();
            }
            REQUIRE (error.find("Something went wrong with processing a number") == 0);
            error = "";
            try{
                strToDouble("10e3.3");
            }
            catch (std::invalid_argument &e){
                error = e.what();
            }
            REQUIRE (error.find("Something went wrong with processing a number") == 0);
            error = "";
        
    }
    SECTION("Parse boolean string"){
        REQUIRE(parseBoolString(" t ") == true);
        REQUIRE(parseBoolString(" tRuE") == true);
        REQUIRE(parseBoolString(" 1 ") == true);
        REQUIRE(parseBoolString("yES ") == true);
        REQUIRE(parseBoolString(" f ") == false);
        REQUIRE(parseBoolString(" FAlse") == false);
        REQUIRE(parseBoolString("0") == false);
        REQUIRE(parseBoolString("no") == false);
        string err;
        try{
            parseBoolString("");
        }
        catch(std::invalid_argument &e){
            err = e.what();
        }
        REQUIRE(err.find("Unable to read boolean value") == 0);
    }
    SECTION("Parse a scalar parameter"){
        REQUIRE(parseScalarParam(" 20 ") == 20);
        REQUIRE(parseScalarParam("-30.5") == -30.5);
        REQUIRE(parseScalarParam(" 20.0E-09 ") == 20.0e-9);

        // Both formats the removed parameter sweep accepted are rejected, and print
        // the migration message.
        REQUIRE_THROWS_WITH(
            parseScalarParam("200 220 240", "Pressure [hPa] (double)"),
            Catch::Matchers::ContainsSubstring("Several values given at Pressure [hPa] (double)") &&
            Catch::Matchers::ContainsSubstring("exactly one simulation per process") &&
            Catch::Matchers::ContainsSubstring("one input file per value")
        );
        REQUIRE_THROWS_WITH(
            parseScalarParam("200:20:240", "Pressure [hPa] (double)"),
            Catch::Matchers::ContainsSubstring("Several values given at Pressure [hPa] (double)") &&
            Catch::Matchers::ContainsSubstring("one input file per value")
        );
        // The Monte Carlo min:max format goes the same way.
        REQUIRE_THROWS_WITH(
            parseScalarParam("200:240"),
            Catch::Matchers::ContainsSubstring("one input file per value")
        );

        REQUIRE_THROWS_WITH(
            parseScalarParam("asdf"),
            Catch::Matchers::ContainsSubstring("Something went wrong with processing a number")
        );
    }
    SECTION("Parse a list of values"){
        Vector_1D vec = parseVectorDoubleString(" 20 30.4 40.0 50 ");
        REQUIRE(vec.size() == 4);
        REQUIRE(vec[0] == 20);
        REQUIRE(vec[1] == 30.4);
        REQUIRE(vec[3] == 50);

        // A single value is still a list of one
        REQUIRE(parseVectorDoubleString("1").size() == 1);

        REQUIRE_THROWS_WITH(
            parseVectorDoubleString("1 asdf ed 3"),
            Catch::Matchers::ContainsSubstring("Something went wrong with processing a number")
        );
    }
    SECTION("Parse a list of integer indices"){
        // Species and aerosol timeseries indices keep the space-separated form
        vector<int> indices = parseVectorIntString(" 1 3 5 ");
        REQUIRE(indices.size() == 3);
        REQUIRE(indices[0] == 1);
        REQUIRE(indices[1] == 3);
        REQUIRE(indices[2] == 5);

        REQUIRE(parseVectorIntString("1") == vector<int>{1});

        REQUIRE_THROWS_WITH(
            parseVectorIntString("1 2.5", "Species indices to include (list of ints)"),
            Catch::Matchers::ContainsSubstring("Decimals not allowed in int inputs")
        );
    }
}
TEST_CASE("Read Yaml File"){
    string filename = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test.yaml";

    YAML::Node data = YAML::LoadFile(filename);
    // The menu readers are called here without going through mergeYamlInputFiles, which
    // is what normally resolves paths. readPath rejects a path left relative, so do that
    // one step by hand.
    resolvePathsInPlace(data, std::filesystem::path(filename).parent_path());
    SECTION("Read Simulation Menu"){
        OptInput input;
        readSimMenu(input, data["SIMULATION MENU"]);
        // If in DEBUG, SIMULATION_OMP_NUM_THREADS is set to 1 which fails the test
        #ifndef DEBUG
            REQUIRE(input.SIMULATION_OMP_NUM_THREADS == 8);
        #endif
        REQUIRE(input.SIMULATION_OVERWRITE == true);
    }
    SECTION("Read Param Menu"){
        // readParamMenu fills an Input object directly, with hPa converted to
        // Pa and the SO2 to SO4 conversion turned from a percentage to a ratio.
        Input scenario;
        readParamMenu(scenario, data["PARAMETER MENU"]);

        REQUIRE(scenario.simulationTime() == 24);
        REQUIRE(scenario.pressure_Pa() == 22000);
        REQUIRE(scenario.horizDiff() == 15.0);
        REQUIRE(scenario.vertiDiff() == 0.15);
        REQUIRE(scenario.nBV() == 0.013);

        REQUIRE(scenario.longitude_deg() == -15);
        REQUIRE(scenario.latitude_deg() == 60);
        REQUIRE(scenario.emissionDOY() == 81);
        REQUIRE(scenario.emissionTime() == 8);

        REQUIRE(scenario.backgNOx() == 5100);
        REQUIRE(scenario.backgHNO3() == 81.5);
        REQUIRE(scenario.backgO3() == 100);
        REQUIRE(scenario.backgCO() == 40);
        REQUIRE(scenario.backgCH4() == 1.76);
        REQUIRE(scenario.backgSO2() == 7.25);

        REQUIRE(scenario.EI_SO2() == 0.1);
        REQUIRE(scenario.EI_SO2TOSO4() == 0.05);
        REQUIRE(scenario.EI_Soot() == 0.06);

        REQUIRE(scenario.sootRad() == 20.0e-9);
        REQUIRE(scenario.fuelFlow() == 2.8);
        REQUIRE(scenario.aircraftMass() == 2.00e5);
        REQUIRE(scenario.flightSpeed() == 250.0);
        REQUIRE(scenario.numEngines() == 4);
        REQUIRE(scenario.wingspan() == 69.8);
        REQUIRE(scenario.coreExitTemp() == 547.3);
        REQUIRE(scenario.bypassArea() == 1.804);
    }
    SECTION("Read Transport Menu"){
        OptInput input;
        readTransportMenu(input, data["TRANSPORT MENU"]);
        REQUIRE(input.TRANSPORT_TRANSPORT == true);
        REQUIRE(input.TRANSPORT_TIMESTEP == 10);
    }
    SECTION("Read Aerosol Menu"){
        OptInput input;
        readAeroMenu(input, data["AEROSOL MENU"]);
        REQUIRE(input.AEROSOL_GRAVSETTLING == true);
        REQUIRE(input.AEROSOL_COAGULATION_SOLID == true);
        REQUIRE(input.AEROSOL_ICE_GROWTH == true);
        REQUIRE(input.AEROSOL_ICE_GROWTH_TIMESTEP == 10);
    }
    SECTION("Read Met Menu"){
        OptInput input;
        string error;
        try{
            readMetMenu(input, data["METEOROLOGY MENU"]);
        }
        catch(std::invalid_argument &e){
            error = e.what();
        }
    
        REQUIRE(input.MET_FILENAME == "/path/to/met/input");
        REQUIRE(input.MET_DT == 1.0);
        REQUIRE(input.MET_TEMPTIMESERIES == true);
        REQUIRE(input.MET_INTERPTEMPDATA == true);
        REQUIRE(input.MET_RHTIMESERIES == true);
        REQUIRE(input.MET_INTERPRHDATA == true);
        REQUIRE(input.MET_SHEARTIMESERIES == true);
        REQUIRE(input.MET_INTERPSHEARDATA == true);
        REQUIRE(input.MET_VERTVELOCTIMESERIES == true);
        REQUIRE(input.MET_INTERPVERTVELOC == true);
        REQUIRE(input.MET_ENABLE_TEMP_PERTURB == true);
        REQUIRE(input.MET_TEMP_PERTURB_AMPLITUDE == 2.0);

    }
    SECTION("Read Diagnostic Menu"){
        OptInput input;
        readDiagMenu(input, data["DIAGNOSTIC MENU"]);
        REQUIRE(input.TS_AERO == true);
        REQUIRE(input.TS_AERO_FILENAME == "ts_aerosol_hhmm.nc");
        REQUIRE(input.TS_AERO_FREQ == 10);
    }
    SECTION("Read Advanced Options Menu") {
        OptInput input;
        readAdvancedMenu(input, data["ADVANCED OPTIONS MENU"]);
        REQUIRE(input.ADV_GRID_NX == 200);
        REQUIRE(input.ADV_GRID_NY == 180);
        REQUIRE(input.ADV_GRID_XLIM_LEFT == 1.0e+3);
        REQUIRE(input.ADV_GRID_XLIM_RIGHT == 1.0e+3);
        REQUIRE(input.ADV_GRID_YLIM_DOWN == 1.5e+3);
        REQUIRE(input.ADV_GRID_YLIM_UP == 300);
        REQUIRE(input.ADV_AMBIENT_LAPSERATE == -3.0);
        REQUIRE(input.ADV_TROPOPAUSE_PRESSURE == 2.0e+4);

    }

}
TEST_CASE("Read one input file into a scenario"){
    string filename = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test1.yaml";
    OptInput input;
    Input scenario;
    YAML::Node merged = YamlInputReader::mergeYamlInputFiles({filename});
    YamlInputReader::populateInput(input, scenario, merged);
    REQUIRE(scenario.simulationTime() == 24);
    REQUIRE(scenario.pressure_Pa() == 22000);
    // Not set by test1.yaml, so it keeps the compiled default
    REQUIRE(scenario.horizDiff() == 15.0);
    REQUIRE(scenario.vertiDiff() == 0.15);
    REQUIRE(scenario.nBV() == 0.013);
    REQUIRE(scenario.longitude_deg() == -15);
    REQUIRE(scenario.latitude_deg() == 60);
    REQUIRE(scenario.emissionDOY() == 81);
    REQUIRE(scenario.emissionTime() == 8);
    REQUIRE(scenario.backgNOx() == 5100);
    REQUIRE(scenario.backgHNO3() == 81.5);
    REQUIRE(scenario.backgO3() == 100);
    REQUIRE(scenario.backgCO() == 40);
    REQUIRE(scenario.backgCH4() == 1.76);
    REQUIRE(scenario.backgSO2() == 7.25);
    REQUIRE(scenario.EI_SO2() == 0.1);
    REQUIRE(scenario.EI_SO2TOSO4() == 0.05);
    REQUIRE(scenario.EI_Soot() == 0.06);
    REQUIRE(scenario.sootRad() == 20.0e-9);
    REQUIRE(scenario.fuelFlow() == 2.8);
    REQUIRE(scenario.aircraftMass() == 2.00e5);
    REQUIRE(scenario.flightSpeed() == 250.0);
    REQUIRE(scenario.numEngines() == 4);
    REQUIRE(scenario.wingspan() == 69.8);
    REQUIRE(scenario.coreExitTemp() == 547.3);
    REQUIRE(scenario.bypassArea() == 1.804);
}

TEST_CASE("Reject input files written for the removed multi-case runs"){
    OptInput input;
    Input scenario;

    SECTION("A PARAM SWEEP SUBMENU names the removed option"){
        // A removed key is a key validation error, so the merge phase raises it.
        string filename = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test15.yaml";
        REQUIRE_THROWS_WITH(
            YamlInputReader::mergeYamlInputFiles({filename}),
            Catch::Matchers::ContainsSubstring("Removed option found") &&
            Catch::Matchers::ContainsSubstring("PARAM SWEEP SUBMENU") &&
            Catch::Matchers::ContainsSubstring("exactly one simulation per process") &&
            Catch::Matchers::ContainsSubstring("one input file per value") &&
            Catch::Matchers::ContainsSubstring("test15.yaml") &&
            !Catch::Matchers::ContainsSubstring("Unknown key found")
        );
    }

    SECTION("A swept PARAMETER MENU entry gets the same message"){
        // The key is valid, but its value is a swept list, so the merge phase
        // accepts the file and the populate phase raises the error.
        string filename = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test16.yaml";
        YAML::Node merged;
        REQUIRE_NOTHROW(merged = YamlInputReader::mergeYamlInputFiles({filename}));
        REQUIRE_THROWS_WITH(
            YamlInputReader::populateInput(input, scenario, merged),
            Catch::Matchers::ContainsSubstring("Several values given at Pressure [hPa] (double)") &&
            Catch::Matchers::ContainsSubstring("exactly one simulation per process") &&
            Catch::Matchers::ContainsSubstring("one input file per value")
        );
    }

}

TEST_CASE("Reject input files naming an option that was removed for being unused"){
    // Test that obsolete options are correctly reporting an error with a descriptive error message
    SECTION("A removed leaf key is named in the error"){
        string filename = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test17.yaml";
        REQUIRE_THROWS_WITH(
            YamlInputReader::mergeYamlInputFiles({filename}),
            Catch::Matchers::ContainsSubstring("Removed option found") &&
            Catch::Matchers::ContainsSubstring("SIMULATION MENU -> Use threaded FFT (T/F)") &&
            Catch::Matchers::ContainsSubstring("never changed a result") &&
            Catch::Matchers::ContainsSubstring("test17.yaml") &&
            // Check that the message overrides the default "Unknow key found"
            !Catch::Matchers::ContainsSubstring("Unknown key found")
        );
    }
}

TEST_CASE("mergeYamlNodes keeps each key exactly once"){
    // Minimal test before using real files
    // Override contains the same keys as defaults, ensure that only one remains after merge
    // and that it is the value of the "later" file (the override)
    YAML::Node defaults = YAML::Load("MENU:\n  alpha: 1\n  beta: 2\n  gamma: 3\n");
    YAML::Node overrides = YAML::Load("MENU:\n  beta: 20\n  delta: 4\n");

    YAML::Node merged = mergeYamlNodes(defaults, overrides);

    REQUIRE(merged.size() == 1);          // one MENU, not two
    REQUIRE(merged["MENU"].size() == 4);  // alpha beta gamma delta, not 6

    // Check that later file wins, earlier-only keys survive.
    REQUIRE(merged["MENU"]["alpha"].as<int>() == 1);   // defaults only
    REQUIRE(merged["MENU"]["beta"].as<int>() == 20);   // overridden
    REQUIRE(merged["MENU"]["gamma"].as<int>() == 3);   // defaults only
    REQUIRE(merged["MENU"]["delta"].as<int>() == 4);   // overrides only
}

TEST_CASE("mergeYamlNodes rejects non-scalar keys"){
    // YAML's allows for non scalar keys.
    // APCEMM inputs never use it, and the merge matches keys by string representation
    // of the key. If we had a non scalar key, the merging would break
    YAML::Node plain = YAML::Load("MENU:\n  alpha: 1\n");
    YAML::Node complexKey = YAML::Load("? [a, b]\n: 1\n");

    // Bad key in overrides
    REQUIRE_THROWS_WITH(
        mergeYamlNodes(plain, complexKey),
        Catch::Matchers::ContainsSubstring("map keys must be scalars")
    );

    // Bad key in defaults
    REQUIRE_THROWS_WITH(
        mergeYamlNodes(complexKey, plain),
        Catch::Matchers::ContainsSubstring("map keys must be scalars")
    );

    // A nested non-scalar key is caught too.
    YAML::Node nested = YAML::Load("MENU:\n  ? {x: 1}\n  : 2\n");
    REQUIRE_THROWS_WITH(
        mergeYamlNodes(plain, nested),
        Catch::Matchers::ContainsSubstring("map keys must be scalars")
    );
}

TEST_CASE("mergeYamlInputFiles resolves defaults and overrides into one node"){
    string filename1 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test1.yaml";
    string filename2 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test2.yaml";
    YAML::Node merged = YamlInputReader::mergeYamlInputFiles({filename1, filename2});

    YAML::Node metParams = merged["PARAMETER MENU"]["METEOROLOGICAL PARAMETERS SUBMENU"];
    // Set by both files: the later one wins.
    REQUIRE(metParams["Pressure [hPa] (double)"].as<double>() == 320);
    // Set by the first file only
    REQUIRE(metParams["Verti. diff. [m^2/s] (double)"].as<double>() == 0.15);
    // Set by the second file only
    REQUIRE(metParams["Horiz. diff. coeff. [m^2/s] (double)"].as<double>() == 17.0);
    // Set by neither file, so the compiled-in default stays
    REQUIRE(merged["SIMULATION MENU"]["External EPM NetCDF file"].as<string>() == "=MISSING=");
    REQUIRE(merged["SIMULATION MENU"]["RANDOM NUMBER GENERATION SUBMENU"]["Force seed value (T/F)"].as<string>() == "F");
}

TEST_CASE("recordEffectiveSeed updates seed submenu"){
    string filename = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test1.yaml";
    YAML::Node merged = YamlInputReader::mergeYamlInputFiles({filename});
    // Verify that the default is an unforced seed otherwise the rest of the test is useless
    REQUIRE(merged["SIMULATION MENU"]["RANDOM NUMBER GENERATION SUBMENU"]["Force seed value (T/F)"].as<string>() == "F");

    SECTION("The seed submenu is rewritten in place"){
        YamlInputReader::recordEffectiveSeed(merged, 1755792000);

        // Verify that the node is correctly updated
        YAML::Node seedSubmenu = merged["SIMULATION MENU"]["RANDOM NUMBER GENERATION SUBMENU"];
        REQUIRE(seedSubmenu["Force seed value (T/F)"].as<string>() == "T");
        REQUIRE(seedSubmenu["Seed value (positive int)"].as<string>() == "1755792000");

        // Sanity check that if we recreate an input from the new node
        // the seed menu is correct
        OptInput input;
        Input scenario;
        YamlInputReader::populateInput(input, scenario, merged);
        REQUIRE(input.SIMULATION_FORCE_SEED == true);
        REQUIRE(input.SIMULATION_SEED_VALUE == 1755792000);
    }
}

TEST_CASE("writeYaml round trip"){
    string filename = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test1.yaml";
    YAML::Node merged = YamlInputReader::mergeYamlInputFiles({filename});

    std::filesystem::path outputDir = std::filesystem::temp_directory_path() / "apcemm-writeyaml-test";
    std::filesystem::remove_all(outputDir);

    SECTION("Round trip through the file"){
        std::filesystem::create_directories(outputDir);
        YamlInputReader::writeYaml(merged, outputDir, OutputFiles::MERGED_YAML);

        std::filesystem::path written = outputDir / OutputFiles::MERGED_YAML;
        REQUIRE(std::filesystem::exists(written));

        YAML::Node reloaded = YAML::LoadFile(written.string());
        REQUIRE(YAML::Dump(reloaded) == YAML::Dump(merged));
    }

    SECTION("Missing output directory raises"){
        REQUIRE_THROWS_WITH(
            YamlInputReader::writeYaml(merged, outputDir, OutputFiles::MERGED_YAML),
            Catch::Matchers::ContainsSubstring("Could not open") &&
            Catch::Matchers::ContainsSubstring(OutputFiles::MERGED_YAML)
        );
    }

    std::filesystem::remove_all(outputDir);
}

TEST_CASE("Merge Input Files"){
    string filename1 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test1.yaml";
    string filename2 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test2.yaml";
    OptInput input;
    Input scenario;
    YAML::Node merged = YamlInputReader::mergeYamlInputFiles({filename1, filename2});
    YamlInputReader::populateInput(input, scenario, merged);
    REQUIRE(scenario.simulationTime() == 24);
    REQUIRE(scenario.pressure_Pa() == 32000);
    REQUIRE(scenario.horizDiff() == 17.0);
    REQUIRE(scenario.vertiDiff() == 0.15);
    REQUIRE(scenario.nBV() == 0.017);
    REQUIRE(scenario.longitude_deg() == -45);
    REQUIRE(scenario.latitude_deg() == 65);
    REQUIRE(scenario.emissionDOY() == 143);
    REQUIRE(scenario.emissionTime() == 12);
    REQUIRE(scenario.backgNOx() == 5100);
    REQUIRE(scenario.backgHNO3() == 81.5);
    REQUIRE(scenario.backgO3() == 100);
    REQUIRE(scenario.backgCO() == 40);
    REQUIRE(scenario.backgCH4() == 1.76);
    REQUIRE(scenario.backgSO2() == 7.25);
    REQUIRE(scenario.EI_SO2() == 0.1);
    REQUIRE(scenario.EI_SO2TOSO4() == 0.05);
    REQUIRE(scenario.EI_Soot() == 0.06);
    REQUIRE(scenario.sootRad() == 20.0e-9);
    REQUIRE(scenario.fuelFlow() == 2.9);
    REQUIRE(scenario.aircraftMass() == 2.30e5);
    REQUIRE(scenario.flightSpeed() == 250.0);
    REQUIRE(scenario.numEngines() == 4);
    REQUIRE(scenario.wingspan() == 69.8);
    REQUIRE(scenario.coreExitTemp() == 547.3);
    REQUIRE(scenario.bypassArea() == 1.804);
}

TEST_CASE("Validate Input Files"){
    OptInput input;
    Input scenario;
    string validFile = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test1.yaml";
    string filename1 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test3.yaml";
    string filename2 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test4.yaml";
    string filename3 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test5.yaml";

    SECTION("Invalid key (scalar) at the root level") {
        // Check that it detects the invalid key, that it points to the correct file and that
        // it prints out the name of the invalid key (here a scalar)
        REQUIRE_THROWS_WITH(
            YamlInputReader::mergeYamlInputFiles({validFile, filename1}),
            Catch::Matchers::ContainsSubstring("Unknown key found") &&
            Catch::Matchers::ContainsSubstring("test3.yaml") &&
            Catch::Matchers::ContainsSubstring("INVALID YAML INPUT")
    );
    }

    SECTION("Invalid key (map) in a submenu"){
        // Check that it detects the invalid key and that it prints out the name
        // of the invalid key (here a map)
        REQUIRE_THROWS_WITH(
            YamlInputReader::mergeYamlInputFiles({filename2}),
            Catch::Matchers::ContainsSubstring("Unknown key found") &&
            Catch::Matchers::ContainsSubstring("INVALID YAML KEY")
        );
    }

    SECTION("Valid key but wrong type (map instead of scalar)"){
        // Here we have a key that is supposed to be a scalar but instead is a map
        // Check that if detects this and prints the name correctly
        REQUIRE_THROWS_WITH(
            YamlInputReader::mergeYamlInputFiles({filename3}),
            Catch::Matchers::ContainsSubstring("is a map in provided YAML but not in the default input.yaml") &&
            Catch::Matchers::ContainsSubstring("Met input file path (string)")
        );
    }

    string filename4 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test6.yaml";
    string filename5 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test7.yaml";
    string filename6 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test8.yaml";
    string filename7 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test9.yaml";

    SECTION("Valid key but wrong type (scalar instead of map)"){
        // A scalar where the default holds a submenu would replace the whole
        // submenu. Check that if detects this and prints the name correctly
        REQUIRE_THROWS_WITH(
            YamlInputReader::mergeYamlInputFiles({filename4}),
            Catch::Matchers::ContainsSubstring("is a value in provided YAML but a map in the default input.yaml") &&
            Catch::Matchers::ContainsSubstring("SIMULATION MENU -> OUTPUT SUBMENU") &&
            Catch::Matchers::ContainsSubstring("test6.yaml")
        );
    }

    SECTION("Valid key but wrong type (sequence instead of map)"){
        REQUIRE_THROWS_WITH(
            YamlInputReader::mergeYamlInputFiles({filename5}),
            Catch::Matchers::ContainsSubstring("is a value in provided YAML but a map in the default input.yaml") &&
            Catch::Matchers::ContainsSubstring("SIMULATION MENU")
        );
    }

    SECTION("Null values mean no override and stay valid"){
        YAML::Node merged;
        REQUIRE_NOTHROW(merged = YamlInputReader::mergeYamlInputFiles({filename6}));
        REQUIRE_NOTHROW(YamlInputReader::populateInput(input, scenario, merged));
        REQUIRE(input.SIMULATION_FORCE_SEED == false);
        REQUIRE(input.SIMULATION_SEED_VALUE == 0);
    }

    SECTION("Empty input file stays valid"){
        // The file stays empty on purpose, so it cannot set the required output folder
        // and cannot reach populateInput. What it tests is that the merge accepts it and
        // leaves every compiled default in place.
        YAML::Node merged;
        REQUIRE_NOTHROW(merged = YamlInputReader::mergeYamlInputFiles({filename7}));
        REQUIRE(merged["SIMULATION MENU"]["RANDOM NUMBER GENERATION SUBMENU"]["Seed value (positive int)"].as<string>() == "0");
        REQUIRE(merged["SIMULATION MENU"]["OUTPUT SUBMENU"]["Output folder (string)"].as<string>() == "=MISSING=");
    }

    SECTION("Document root is a bare value"){
        string filename8 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test10.yaml";
        REQUIRE_THROWS_WITH(
            YamlInputReader::mergeYamlInputFiles({filename8}),
            Catch::Matchers::ContainsSubstring("The document root is a value in provided YAML") &&
            Catch::Matchers::ContainsSubstring("test10.yaml")
        );
    }

    SECTION("Repeated menu heading at the root"){
        // yaml-cpp keeps both entries but reads the first one, so the second
        // heading would be dropped silently.
        string filename9 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test11.yaml";
        REQUIRE_THROWS_WITH(
            YamlInputReader::mergeYamlInputFiles({filename9}),
            Catch::Matchers::ContainsSubstring("Duplicate key found") &&
            Catch::Matchers::ContainsSubstring("SIMULATION MENU") &&
            Catch::Matchers::ContainsSubstring("test11.yaml")
        );
    }

    SECTION("Non-scalar key at the root"){
        // The validator reads every key as a string. It must reject a non-scalar
        // key with its own message, not with a yaml-cpp "bad conversion".
        string filename11 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test13.yaml";
        REQUIRE_THROWS_WITH(
            YamlInputReader::mergeYamlInputFiles({filename11}),
            Catch::Matchers::ContainsSubstring("map keys must be scalars") &&
            Catch::Matchers::ContainsSubstring("test13.yaml")
        );
    }

    SECTION("Non-scalar key inside a submenu"){
        // Same check, one level down, where the recursion reaches the key.
        string filename12 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test14.yaml";
        REQUIRE_THROWS_WITH(
            YamlInputReader::mergeYamlInputFiles({filename12}),
            Catch::Matchers::ContainsSubstring("map keys must be scalars") &&
            Catch::Matchers::ContainsSubstring("test14.yaml")
        );
    }

    SECTION("Repeated key inside a submenu"){
        // The error names the full path of the key, not just the key itself.
        string filename10 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test12.yaml";
        REQUIRE_THROWS_WITH(
            YamlInputReader::mergeYamlInputFiles({filename10}),
            Catch::Matchers::ContainsSubstring("Duplicate key found") &&
            Catch::Matchers::ContainsSubstring("SIMULATION MENU -> OUTPUT SUBMENU -> Output folder (string)") &&
            Catch::Matchers::ContainsSubstring("test12.yaml")
        );
    }
}

namespace {
    namespace fs = std::filesystem;

    // Test that "a relative path resolves against the directory of the file
    // that declares it", so the test needs input files in more than one directory. 
    // Write these test file in /tmp/
    struct TwoDirectoryInputs {
        fs::path root;

        explicit TwoDirectoryInputs(const string& name)
            : root(fs::temp_directory_path() / name) {
            fs::remove_all(root);
            fs::create_directories(root / "caseA" / "data");
            fs::create_directories(root / "caseB");
        }
        // Ensure clean up test files in /tmp/
        ~TwoDirectoryInputs() { fs::remove_all(root); }

        void write(const string& relativePath, const string& contents) const {
            std::ofstream out(root / relativePath);
            out << contents;
        }
        string path(const string& relativePath) const {
            return (root / relativePath).string();
        }
        // What the resolution rule must produce for a path written in a given file.
        string resolved(const string& relativePath) const {
            return fs::weakly_canonical(root / relativePath).generic_string();
        }
    };

    // Get a node value from menu/submenu/key keys
    string pathIn(const YAML::Node& merged, const string& submenu, const string& key) {
        return merged["SIMULATION MENU"][submenu][key].as<string>();
    }
}

TEST_CASE("Relative path resolution"){
    // Setup test dir
    TwoDirectoryInputs inputs("apcemm-path-resolution-test");
    // Temp data
    inputs.write("caseA/data/init.txt", "background\n");
    // First input file in directory caseA/
    inputs.write("caseA/base.yaml",
        "SIMULATION MENU:\n"
        "  OUTPUT SUBMENU:\n"
        "    Output folder (string): out/\n"
        "  Input background condition (string): data/init.txt\n");
    // Second input file is in directory caseB/
    inputs.write("caseB/override.yaml",
        "PARAMETER MENU:\n"
        "  Plume Process [hr] (double): 3\n");

    // The input files are in different directories to test relative path resolution
    const string base = inputs.path("caseA/base.yaml");
    const string overrideFile = inputs.path("caseB/override.yaml");

    SECTION("Path resolution works"){
        // Before the resolution used to resolve with the directory of the last file
        YAML::Node merged = YamlInputReader::mergeYamlInputFiles({base, overrideFile});
        REQUIRE(merged["SIMULATION MENU"]["Input background condition (string)"].as<string>()
                == inputs.resolved("caseA/data/init.txt"));
        // 'out/' is written with a trailing '/', resolution drops it.
        REQUIRE(pathIn(merged, "OUTPUT SUBMENU", "Output folder (string)")
                == inputs.resolved("caseA/out"));
    }

    SECTION("Each file's own paths resolve against its own directory"){
        inputs.write("caseB/override.yaml",
            "SIMULATION MENU:\n"
            "  OUTPUT SUBMENU:\n"
            "    Output folder (string): out/\n");
        YAML::Node merged = YamlInputReader::mergeYamlInputFiles({base, overrideFile});
        // Declared by the first file, resolved under the first file's directory.
        REQUIRE(merged["SIMULATION MENU"]["Input background condition (string)"].as<string>()
                == inputs.resolved("caseA/data/init.txt"));
        // Declared by both files, so the last one overrides the first (normal override behavior)
        // and resolves against the second input file dir
        REQUIRE(pathIn(merged, "OUTPUT SUBMENU", "Output folder (string)")
                == inputs.resolved("caseB/out"));
    }

    SECTION("Absolute paths and the sentinels are not changed"){
        const string absolute = inputs.resolved("caseA/data/init.txt");
        inputs.write("caseA/absolute.yaml",
            "SIMULATION MENU:\n"
            "  OUTPUT SUBMENU:\n"
            "    Output folder (string): " + absolute + "\n"
            "  Input background condition (string): " + absolute + "\n"
            "  Input engine emissions (string): =DEFAULT=\n"
            "  External EPM NetCDF file: =MISSING=\n");
        YAML::Node merged = YamlInputReader::mergeYamlInputFiles({inputs.path("caseA/absolute.yaml")});
        REQUIRE(pathIn(merged, "OUTPUT SUBMENU", "Output folder (string)") == absolute);
        REQUIRE(merged["SIMULATION MENU"]["Input background condition (string)"].as<string>() == absolute);
        REQUIRE(merged["SIMULATION MENU"]["Input engine emissions (string)"].as<string>() == "=DEFAULT=");
        REQUIRE(merged["SIMULATION MENU"]["External EPM NetCDF file"].as<string>() == "=MISSING=");
    }
}

TEST_CASE("readPath rejects a path that was never resolved"){
    // Check that readPath correctly guards against left over relative paths
    YAML::Node node = YAML::Load(
        "absolute: /tmp/somewhere\n"
        "relative: data/init.txt\n"
        "missing: =MISSING=\n"
        "default: =DEFAULT=\n");
    REQUIRE(readPath(node, "absolute") == "/tmp/somewhere");
    REQUIRE(readPath(node, "default") == "=DEFAULT=");
    REQUIRE_THROWS_WITH(
        readPath(node, "relative"),
        Catch::Matchers::ContainsSubstring("still relative") &&
        Catch::Matchers::ContainsSubstring("PATH_KEYS")
    );
}

namespace {
    // Load in the defaults via loading an empty file. This is so that we can mutate
    // the default yaml for tests here
    YAML::Node compiledDefaults() {
        const string emptyFile = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test9.yaml";
        return YamlInputReader::mergeYamlInputFiles({emptyFile});
    }
}

TEST_CASE("checkDefaultPaths guards the compiled-in defaults"){
    const string pathKey = "Input background condition (string)";

    SECTION("Real defaults pass"){
        REQUIRE_NOTHROW(checkDefaultPaths(compiledDefaults()));
    }

    SECTION("A PATH_KEYS entry that is no longer a leaf is rejected"){
        // Renaming a key in defaults/input.yaml means one of the elements in
        // PATH_KEYS is not being used, so the path resolution would not occur
        // for the renamed key.
        // This test ensures that checkDefaultPaths throws if one of the PATH_KEYS
        // is not present in the default file. Most likely, this means that a key
        // was renamed in default/input.yaml but not updated in PATH_KEYS
        YAML::Node defaults = compiledDefaults();
        YAML::Node simMenu = defaults["SIMULATION MENU"];
        REQUIRE(simMenu.remove(pathKey));
        REQUIRE_THROWS_WITH(
            checkDefaultPaths(defaults),
            Catch::Matchers::ContainsSubstring("not a scalar leaf") &&
            Catch::Matchers::ContainsSubstring("SIMULATION MENU -> " + pathKey)
        );
    }

    SECTION("A relative path in the defaults is rejected"){
        // Paths in the defaults cannot be relative, instead they should be
        // =MISSING= or =DEFAULT=. Ensure checkDefaultPaths throws if it
        // finds a relative path
        YAML::Node defaults = compiledDefaults();
        defaults["SIMULATION MENU"][pathKey] = "relative/init.txt";
        REQUIRE_THROWS_WITH(
            checkDefaultPaths(defaults),
            Catch::Matchers::ContainsSubstring("relative path") &&
            Catch::Matchers::ContainsSubstring("relative/init.txt") &&
            Catch::Matchers::ContainsSubstring("SIMULATION MENU -> " + pathKey)
        );
    }
}

TEST_CASE("Not setting output folder is an error"){
    // 'Output folder (string)' is required and the defaults is set to '=MISSING='.
    // Check that populateInput throws if the output dir is unspecified.
    // Without this check APCEMM creates a directory named '=MISSING=' in the
    // working directory.
    OptInput input;
    Input scenario;
    YAML::Node merged = compiledDefaults();
    REQUIRE_THROWS_WITH(
        YamlInputReader::populateInput(input, scenario, merged),
        Catch::Matchers::ContainsSubstring("No output folder set") &&
        Catch::Matchers::ContainsSubstring("Output folder (string)") &&
        Catch::Matchers::ContainsSubstring("=MISSING=")
    );
}
