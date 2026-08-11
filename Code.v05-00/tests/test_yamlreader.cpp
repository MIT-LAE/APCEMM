#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Util/YamlUtils.hpp>
#include <YamlInputReader/YamlInputReader.hpp>
#include <Core/Input.hpp>
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
    SECTION("Parse Parameter Sweep/Monte Carlo Sim input"){
        string teststr;
        Vector_1D vec;
        teststr = " 20:40 ";
        //Test monte carlo case
        vec = parseParamSweepInput(teststr, "", true, 10);
        REQUIRE(vec.size() == 10);
        REQUIRE(vec[0] != vec[1]);
        for (double d: vec){
            REQUIRE (d >= 20);
            REQUIRE (d <= 40);
        }
        //Parameter sweep min:step:max case
        teststr = " 20 : 5.0 : 40 "; // results in 20 25 30 35 40 (length 5)
        vec = parseParamSweepInput(teststr);
        REQUIRE (vec.size() == 5);
        REQUIRE (vec[0] == 20);
        REQUIRE (vec[1] == 25);
        REQUIRE (vec[4] == 40);
        //Standard a b c d .... case
        teststr = " 20 30.4 40.0 50 ";
        vec = parseParamSweepInput(teststr);
        REQUIRE (vec.size() == 4);
        REQUIRE (vec[0] == 20);
        REQUIRE (abs(vec[1] - 30) < 1e-40);
        REQUIRE (vec[3] == 50);

        // some exceptions
        string err;
        try{
            vec = parseParamSweepInput("1:2:3:4", "", true, 10);
        }
        catch (std::invalid_argument &e){
            err = e.what();
        }
        REQUIRE(err.find("Monte Carlo Simulation requires") == 0);
        try{
            vec = parseParamSweepInput("1:2");
        }
        catch (std::invalid_argument &e){
            err = e.what();
        }
        REQUIRE(err.find("Parameter sweep step input requires") == 0);
        try{
            vec = parseParamSweepInput("1 asdf ed 3");
        }
        catch (std::invalid_argument &e){
            err = e.what();
        }
        REQUIRE(err.find("Something went wrong with processing a number") == 0);
        try{
            vec = parseParamSweepInput("10:10:20 0");
        }
        catch(std::invalid_argument &e){
            if(string(e.what()).find("Something went wrong with processing a number") == 0) err = "stod fail 1";
        }
        REQUIRE(err == "stod fail 1");
        try{
            vec = parseParamSweepInput("10:10:20..0");
        }
        catch(std::invalid_argument &e){
            if(string(e.what()).find("Something went wrong with processing a number") == 0) err = "stod fail 2";
        }
        REQUIRE(err == "stod fail 2");
    }
}
TEST_CASE("Read Yaml File"){
    string filename = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test.yaml";

    YAML::Node data = YAML::LoadFile(filename);
    SECTION("Read Simulation Menu"){
        OptInput input;
        string err;
        try{
            readSimMenu(input, data["SIMULATION MENU"]);
        }
        catch (std::invalid_argument &e){
            err = e.what();
        }
        // If in DEBUG, SIMULATION_OMP_NUM_THREADS is set to 1 which fails the test
        #ifndef DEBUG
            REQUIRE(input.SIMULATION_OMP_NUM_THREADS == 8);
        #endif
        REQUIRE(input.SIMULATION_PARAMETER_SWEEP == true);
        REQUIRE(input.SIMULATION_MONTECARLO == true);
        REQUIRE(input.SIMULATION_MCRUNS == 2);
        REQUIRE(input.SIMULATION_OVERWRITE == true);
        REQUIRE(input.SIMULATION_THREADED_FFT == true);
        REQUIRE(input.SIMULATION_USE_FFTW_WISDOM == true);
        REQUIRE(input.SIMULATION_SAVE_FORWARD == true);
        REQUIRE(input.SIMULATION_FORWARD_FILENAME == "APCEMM_Case_*");
        REQUIRE(input.SIMULATION_ADJOINT == true);
        REQUIRE(input.SIMULATION_ADJOINT_FILENAME == "APCEMM_ADJ_Case_*");
        REQUIRE(input.SIMULATION_BOXMODEL == true);
        REQUIRE(input.SIMULATION_BOX_FILENAME == "APCEMM_BOX_CASE_*");
        REQUIRE(err == "In Simulation Menu: Parameter sweep and Monte Carlo cannot have the same value!");

    }
    SECTION("Read Param Menu"){
        OptInput input;
        string err;
        readParamMenu(input, data["PARAMETER MENU"]);

        REQUIRE(input.PARAMETER_PARAM_MAP["PLUMEPROCESS"][0] == 24);
        REQUIRE(input.PARAMETER_PARAM_MAP["PRESSURE"].size() == 3);
        REQUIRE(input.PARAMETER_PARAM_MAP["PRESSURE"][0] == 22000);
        REQUIRE(input.PARAMETER_PARAM_MAP["PRESSURE"][1] == 23000);
        REQUIRE(input.PARAMETER_PARAM_MAP["DH"][0] == 15.0);
        REQUIRE(input.PARAMETER_PARAM_MAP["DV"][0] == 0.15);
        REQUIRE(input.PARAMETER_PARAM_MAP["NBV"][0] == 0.013);

        REQUIRE(input.PARAMETER_PARAM_MAP["LONGITUDE"][0] == -15);
        REQUIRE(input.PARAMETER_PARAM_MAP["LATITUDE"][0] == 60);
        REQUIRE(input.PARAMETER_PARAM_MAP["EDAY"][0] == 81);
        REQUIRE(input.PARAMETER_PARAM_MAP["ETIME"][0] == 8);

        REQUIRE(input.PARAMETER_PARAM_MAP["BACKG_NOX"][0] == 5100);
        REQUIRE(input.PARAMETER_PARAM_MAP["BACKG_HNO3"][0] == 81.5);
        REQUIRE(input.PARAMETER_PARAM_MAP["BACKG_O3"][0] == 100);
        REQUIRE(input.PARAMETER_PARAM_MAP["BACKG_CO"][0] == 40);
        REQUIRE(input.PARAMETER_PARAM_MAP["BACKG_CH4"][0] == 1.76);
        REQUIRE(input.PARAMETER_PARAM_MAP["BACKG_SO2"][0] == 7.25);

        REQUIRE(input.PARAMETER_PARAM_MAP["EI_NOX"][0] == 10);
        REQUIRE(input.PARAMETER_PARAM_MAP["EI_CO"][0] == 1);
        REQUIRE(input.PARAMETER_PARAM_MAP["EI_UHC"][0] == 0.6);
        REQUIRE(input.PARAMETER_PARAM_MAP["EI_SO2"][0] == 0.1);
        REQUIRE(input.PARAMETER_PARAM_MAP["EI_SO2TOSO4"][0] == 0.05);
        REQUIRE(input.PARAMETER_PARAM_MAP["EI_SOOT"][0] == 0.06);

        REQUIRE(input.PARAMETER_PARAM_MAP["EI_SOOTRAD"][0] == 20.0e-9);
        REQUIRE(input.PARAMETER_PARAM_MAP["FF"][0] == 2.8);
        REQUIRE(input.PARAMETER_PARAM_MAP["AMASS"][0] == 2.00e5);
        REQUIRE(input.PARAMETER_PARAM_MAP["FSPEED"][0] == 250.0);
        REQUIRE(input.PARAMETER_PARAM_MAP["NUMENG"][0] == 4);
        REQUIRE(input.PARAMETER_PARAM_MAP["WINGSPAN"][0] == 69.8);
        REQUIRE(input.PARAMETER_PARAM_MAP["COREEXITTEMP"][0] == 547.3);
        REQUIRE(input.PARAMETER_PARAM_MAP["BYPASSAREA"][0] == 1.804);
    }
    SECTION("Read Transport Menu"){
        OptInput input;
        readTransportMenu(input, data["TRANSPORT MENU"]);
        REQUIRE(input.TRANSPORT_TRANSPORT == true);
        REQUIRE(input.TRANSPORT_FILL == true);
        REQUIRE(input.TRANSPORT_TIMESTEP == 10);
        REQUIRE(input.TRANSPORT_UPDRAFT == true);
        REQUIRE(input.TRANSPORT_UPDRAFT_TIMESCALE == 3600);
        REQUIRE(input.TRANSPORT_UPDRAFT_VELOCITY == 5);
    }
    SECTION("Read Chemistry Menu"){
        OptInput input;
        readChemMenu(input, data["CHEMISTRY MENU"]);
        REQUIRE(input.CHEMISTRY_CHEMISTRY == true);
        REQUIRE(input.CHEMISTRY_HETCHEM == true);
        REQUIRE(input.CHEMISTRY_TIMESTEP == 10);
        REQUIRE(input.CHEMISTRY_JRATE_FOLDER == "/net/d04/data/fritzt/APCEMM_Data/J-Rates");
    }
    SECTION("Read Aerosol Menu"){
        OptInput input;
        readAeroMenu(input, data["AEROSOL MENU"]);
        REQUIRE(input.AEROSOL_GRAVSETTLING == true);
        REQUIRE(input.AEROSOL_COAGULATION_SOLID == true);
        REQUIRE(input.AEROSOL_COAGULATION_LIQUID == true);
        REQUIRE(input.AEROSOL_COAGULATION_TIMESTEP == 60);
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
    
        REQUIRE(input.MET_LOADMET == true);
        REQUIRE(input.MET_FILENAME == "/path/to/met/input");
        REQUIRE(input.MET_DT == 1.0);
        REQUIRE(input.MET_LOADTEMP == true);
        REQUIRE(input.MET_TEMPTIMESERIES == true);
        REQUIRE(input.MET_INTERPTEMPDATA == true);
        REQUIRE(input.MET_LOADRH == true);
        REQUIRE(input.MET_RHTIMESERIES == true);
        REQUIRE(input.MET_INTERPRHDATA == true);
        REQUIRE(input.MET_LOADSHEAR == true);
        REQUIRE(input.MET_SHEARTIMESERIES == true);
        REQUIRE(input.MET_INTERPSHEARDATA == true);
        REQUIRE(input.MET_LOADVERTVELOC == true);
        REQUIRE(input.MET_VERTVELOCTIMESERIES == true);
        REQUIRE(input.MET_INTERPVERTVELOC == true);
        REQUIRE(input.MET_ENABLE_TEMP_PERTURB == true);
        REQUIRE(input.MET_TEMP_PERTURB_AMPLITUDE == 2.0);
        REQUIRE(input.MET_TEMP_PERTURB_TIMESCALE == 10);

    }
    SECTION("Read Diagnostic Menu"){
        OptInput input;
        readDiagMenu(input, data["DIAGNOSTIC MENU"]);
        REQUIRE(input.DIAG_FILENAME == "trac_avg.apcemm.hhmm");
        REQUIRE(input.TS_SPEC == true);
        REQUIRE(input.TS_FILENAME == "ts_hhmm.nc");
        REQUIRE(input.TS_SPECIES.size() == 3);
        REQUIRE(input.TS_SPECIES[2] == 3);
        REQUIRE(input.TS_FREQ == 10);
        REQUIRE(input.TS_AERO == true);
        REQUIRE(input.TS_AERO_FILENAME == "ts_aerosol_hhmm.nc");
        REQUIRE(input.TS_AEROSOL.size() == 3);
        REQUIRE(input.TS_AEROSOL[2] == 5);
        REQUIRE(input.TS_AERO_FREQ == 10);
        REQUIRE(input.PL_PL == true);
        REQUIRE(input.PL_O3 == true);
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
        REQUIRE(input.ADV_CSIZE_DEPTH_BASE == 180.0);
        REQUIRE(input.ADV_CSIZE_DEPTH_SCALING_FACTOR == 0.5);
        REQUIRE(input.ADV_CSIZE_WIDTH_BASE == 100.0);
        REQUIRE(input.ADV_CSIZE_WIDTH_SCALING_FACTOR == 0.5);
        REQUIRE(input.ADV_AMBIENT_LAPSERATE == -3.0);
        REQUIRE(input.ADV_TROPOPAUSE_PRESSURE == 2.0e+4);

    }

}
TEST_CASE("Generate All Cases"){
    OptInput input;
    input.PARAMETER_PARAM_MAP = {{"test1", {1}}};
    vector<std::unordered_map<string,double>> combinations = generateCases(input);
    REQUIRE(combinations.size() == 1);

    input.PARAMETER_PARAM_MAP = {{"test1", {1}}, {"test2", {2}}};
    combinations = generateCases(input);
    REQUIRE(combinations.size() == 1);
    REQUIRE(combinations[0]["test1"] == 1);
    REQUIRE(combinations[0]["test2"] == 2);

    input.PARAMETER_PARAM_MAP = {{"test1", {1, 2, 3}}, {"test2", {4}}, {"test3", {5, 6}}};
    combinations = generateCases(input);
    //should result in {(1,4,5), (2,4,5), (3,4,5), (1,4,6), (2,4,6), (3,4,6)}
    REQUIRE(combinations.size() == 6);

    REQUIRE(combinations[0]["test1"] == 1);
    REQUIRE(combinations[0]["test2"] == 4);
    REQUIRE(combinations[0]["test3"] == 5);
    REQUIRE(combinations[1]["test1"] == 2);
    REQUIRE(combinations[1]["test3"] == 5);
    REQUIRE(combinations[5]["test1"] == 3);
    REQUIRE(combinations[5]["test3"] == 6);

    input.PARAMETER_PARAM_MAP = {{"test1", {1, 2, 3}}, {"test2", {4, 0}}, {"test3", {5, 6, 7, 8}}, {"test4", {9, 10, 11}}};
    combinations = generateCases(input);
    REQUIRE(combinations.size() == 72);
}
TEST_CASE("Generate Input Objects"){
    string filename = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test1.yaml";
    OptInput input;
    YamlInputReader::readYamlInputFiles(input, {filename});
    vector<std::unordered_map<string,double>> cases = generateCases(input);
    REQUIRE(cases.size() == 6);
    Input caseInput = Input(0, cases, "", "", "", "", "");
    REQUIRE(caseInput.simulationTime() == 24);
    REQUIRE(caseInput.pressure_Pa() == 22000);
    REQUIRE(caseInput.horizDiff() == 15.0);
    REQUIRE(caseInput.vertiDiff() == 0.15);
    REQUIRE(caseInput.nBV() == 0.013);
    REQUIRE(caseInput.longitude_deg() == -15);
    REQUIRE(caseInput.latitude_deg() == 60);
    REQUIRE(caseInput.emissionDOY() == 81);
    REQUIRE(caseInput.emissionTime() == 8);
    REQUIRE(caseInput.backgNOx() == 5100);
    REQUIRE(caseInput.backgHNO3() == 81.5);
    REQUIRE(caseInput.backgO3() == 100);
    REQUIRE(caseInput.backgCO() == 40);
    REQUIRE(caseInput.backgCH4() == 1.76);
    REQUIRE(caseInput.backgSO2() == 7.25);
    REQUIRE(caseInput.EI_NOx() == 10);
    REQUIRE(caseInput.EI_CO() == 1);
    REQUIRE(caseInput.EI_HC() == 0.6);
    REQUIRE(caseInput.EI_SO2() == 0.1);
    REQUIRE(caseInput.EI_SO2TOSO4() == 0.05);
    REQUIRE(caseInput.EI_Soot() == 0.06);
    REQUIRE(caseInput.sootRad() == 20.0e-9);
    REQUIRE(caseInput.fuelFlow() == 2.8);
    REQUIRE(caseInput.aircraftMass() == 2.00e5);
    REQUIRE(caseInput.flightSpeed() == 250.0);
    REQUIRE(caseInput.numEngines() == 4);
    REQUIRE(caseInput.wingspan() == 69.8);
    REQUIRE(caseInput.coreExitTemp() == 547.3);
    REQUIRE(caseInput.bypassArea() == 1.804);
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

TEST_CASE("Merge Input Files"){
    string filename1 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test1.yaml";
    string filename2 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test2.yaml";
    OptInput input;
    YamlInputReader::readYamlInputFiles(input, {filename1, filename2});
    vector<std::unordered_map<string,double>> cases = generateCases(input);
    REQUIRE(cases.size() == 6);
    Input caseInput = Input(0, cases, "", "", "", "", "");
    REQUIRE(caseInput.simulationTime() == 24);
    REQUIRE(caseInput.pressure_Pa() == 32000);
    REQUIRE(caseInput.horizDiff() == 17.0);
    REQUIRE(caseInput.vertiDiff() == 0.15);
    REQUIRE(caseInput.nBV() == 0.017);
    REQUIRE(caseInput.longitude_deg() == -45);
    REQUIRE(caseInput.latitude_deg() == 65);
    REQUIRE(caseInput.emissionDOY() == 143);
    REQUIRE(caseInput.emissionTime() == 12);
    REQUIRE(caseInput.backgNOx() == 5100);
    REQUIRE(caseInput.backgHNO3() == 81.5);
    REQUIRE(caseInput.backgO3() == 100);
    REQUIRE(caseInput.backgCO() == 40);
    REQUIRE(caseInput.backgCH4() == 1.76);
    REQUIRE(caseInput.backgSO2() == 7.25);
    REQUIRE(caseInput.EI_NOx() == 10);
    REQUIRE(caseInput.EI_CO() == 1);
    REQUIRE(caseInput.EI_HC() == 0.6);
    REQUIRE(caseInput.EI_SO2() == 0.1);
    REQUIRE(caseInput.EI_SO2TOSO4() == 0.05);
    REQUIRE(caseInput.EI_Soot() == 0.06);
    REQUIRE(caseInput.sootRad() == 20.0e-9);
    REQUIRE(caseInput.fuelFlow() == 2.9);
    REQUIRE(caseInput.aircraftMass() == 2.30e5);
    REQUIRE(caseInput.flightSpeed() == 250.0);
    REQUIRE(caseInput.numEngines() == 4);
    REQUIRE(caseInput.wingspan() == 69.8);
    REQUIRE(caseInput.coreExitTemp() == 547.3);
    REQUIRE(caseInput.bypassArea() == 1.804);
}

TEST_CASE("Validate Input Files"){
    OptInput input;
    string validFile = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test1.yaml";
    
    string filename1 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test3.yaml";
    string filename2 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test4.yaml";
    string filename3 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test5.yaml";

    SECTION("Invalid key (scalar) at the root level") {
        // Check that it detects the invalid key, that it points to the correct file and that
        // it prints out the name of the invalid key (here a scalar)
        REQUIRE_THROWS_WITH(
            YamlInputReader::readYamlInputFiles(input, {validFile, filename1}), 
            Catch::Matchers::ContainsSubstring("Unknown key found") &&
            Catch::Matchers::ContainsSubstring("test3.yaml") &&
            Catch::Matchers::ContainsSubstring("INVALID YAML INPUT")
    );
    }

    SECTION("Invalid key (map) in a submenu"){
        // Check that it detects the invalid key and that it prints out the name
        // of the invalid key (here a map)
        REQUIRE_THROWS_WITH(
            YamlInputReader::readYamlInputFiles(input, {filename2}),
            Catch::Matchers::ContainsSubstring("Unknown key found") &&
            Catch::Matchers::ContainsSubstring("INVALID YAML KEY")
        );
    }

    SECTION("Valid key but wrong type (map instead of scalar)"){
        // Here we have a key that is supposed to be a scalar but instead is a map
        // Check that if detects this and prints the name correctly
        REQUIRE_THROWS_WITH(
            YamlInputReader::readYamlInputFiles(input, {filename3}),
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
            YamlInputReader::readYamlInputFiles(input, {filename4}),
            Catch::Matchers::ContainsSubstring("is a value in provided YAML but a map in the default input.yaml") &&
            Catch::Matchers::ContainsSubstring("SIMULATION MENU -> OUTPUT SUBMENU") &&
            Catch::Matchers::ContainsSubstring("test6.yaml")
        );
    }

    SECTION("Valid key but wrong type (sequence instead of map)"){
        REQUIRE_THROWS_WITH(
            YamlInputReader::readYamlInputFiles(input, {filename5}),
            Catch::Matchers::ContainsSubstring("is a value in provided YAML but a map in the default input.yaml") &&
            Catch::Matchers::ContainsSubstring("SIMULATION MENU")
        );
    }

    SECTION("Null values mean no override and stay valid"){
        // Reading must succeed and leave the compiled defaults in place. If a null
        // cleared out the submenu, readSimMenu would fail on the missing keys.
        REQUIRE_NOTHROW(YamlInputReader::readYamlInputFiles(input, {filename6}));
        REQUIRE(input.SIMULATION_FORWARD_FILENAME == "APCEMM_Case_*");
    }

    SECTION("Empty input file stays valid"){
        REQUIRE_NOTHROW(YamlInputReader::readYamlInputFiles(input, {filename7}));
        REQUIRE(input.SIMULATION_FORWARD_FILENAME == "APCEMM_Case_*");
    }

    SECTION("Document root is a bare value"){
        string filename8 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test10.yaml";
        REQUIRE_THROWS_WITH(
            YamlInputReader::readYamlInputFiles(input, {filename8}),
            Catch::Matchers::ContainsSubstring("The document root is a value in provided YAML") &&
            Catch::Matchers::ContainsSubstring("test10.yaml")
        );
    }

    SECTION("Repeated menu heading at the root"){
        // yaml-cpp keeps both entries but reads the first one, so the second
        // heading would be dropped silently.
        string filename9 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test11.yaml";
        REQUIRE_THROWS_WITH(
            YamlInputReader::readYamlInputFiles(input, {filename9}),
            Catch::Matchers::ContainsSubstring("Duplicate key found") &&
            Catch::Matchers::ContainsSubstring("SIMULATION MENU") &&
            Catch::Matchers::ContainsSubstring("test11.yaml")
        );
    }

    SECTION("Repeated key inside a submenu"){
        // The error names the full path of the key, not just the key itself.
        string filename10 = string(APCEMM_TESTS_DIR) + YAML_DIR + "/test12.yaml";
        REQUIRE_THROWS_WITH(
            YamlInputReader::readYamlInputFiles(input, {filename10}),
            Catch::Matchers::ContainsSubstring("Duplicate key found") &&
            Catch::Matchers::ContainsSubstring("SIMULATION MENU -> OUTPUT SUBMENU -> Output folder (string)") &&
            Catch::Matchers::ContainsSubstring("test12.yaml")
        );
    }
}
