#include "YamlInputReader/YamlInputReader.hpp"
#include "Core/Input_Mod.hpp"
#include "YamlInputReader/YamlPathUtils.hpp"
#include "APCEMM.h"
#include "Util/ForwardDecl.hpp"
#include "Util/YamlUtils.hpp"
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <set>
#include <string>
#include <vector>


// Read default configuration from CMake-generated include file.
const std::string default_input =
#include "Defaults/Input.hpp"
;


bool ichar_equals(char a, char b) {
  return std::tolower(static_cast<unsigned char>(a)) ==
         std::tolower(static_cast<unsigned char>(b));
}

bool iequals(std::string_view lhs, std::string_view rhs) {
  return std::ranges::equal(lhs, rhs, ichar_equals);
}


namespace YamlInputReader{
    // Helper function to get all keys from a YAML Map node
    std::set<std::string> getYamlKeys(const YAML::Node& node) {
        std::set<std::string> keys;
        if (!node.IsMap()) {
            return keys;
        }
        for (const auto& it : node) {
            keys.insert(getScalarKey(it.first));
        }
        return keys;
    }

    // Reject a map node that uses the same key twice.
    // YAML standard does not allow it, but yaml-cpp keeps both entries
    // while node[key] returns only the first one.
    // This will be fixed in a newer version of yaml-cpp (post 2026/03/12), but
    // we use v0.8.0 which is much older.
    // This means the user must group all the entries of a menu under a single heading.
    // "source" names the file the node comes from, so the message points at the
    // right file.
    void checkNoDuplicateKeys(const YAML::Node& node, const std::string& source, const std::string& currentPath = "") {
        std::set<std::string> seenKeys;
        for (const auto& it : node) {
            const std::string key = getScalarKey(it.first);
            const bool isNewKey = seenKeys.insert(key).second;
            if (!isNewKey) {
                std::string errorPath = buildKeyPath(currentPath, key);
                throw std::runtime_error("Duplicate key found in " + source + ": '" + errorPath + "'. Each key must appear only once. Group all the entries of a menu under a single heading.");
            }
        }
    }

    // Same check, applied to every map of a document.
    // Iterate the entries directly instead of node[key]: with a repeated key,
    // node[key] only ever returns the first entry, so the second subtree would
    // never be visited.
    void checkNoDuplicateKeysRecursive(const YAML::Node& node, const std::string& source, const std::string& currentPath = "") {
        if (!node.IsMap()) {
            return;
        }
        checkNoDuplicateKeys(node, source, currentPath);
        for (const auto& it : node) {
            const std::string key = getScalarKey(it.first);
            const std::string nextPath = buildKeyPath(currentPath, key);
            checkNoDuplicateKeysRecursive(it.second, source, nextPath);
        }
    }

    // Keys that previous versions accepted and that we now ignore with a warning. Checked
    // before the generic unknown-key error so an outdated input file gets a warning
    // naming the option that is deprecated instead of failing with "Unknown key found".
    bool checkDeprecatedKey(const std::string& key, const std::string& errorPath) {
        static const std::set<std::string> deprecatedKeys = {
            "Chemistry Timestep [min] (double)",
            "Coag. timestep [min] (double)",
            "Temp. Perturb. Timescale (min)",
        };
        if (deprecatedKeys.contains(key)) {
            std::cout << "WARNING: Deprecated option found: '" << errorPath << "'. This option is no longer used and has no effect." << std::endl;
            return true;
        }
        return false;
    }

    // Keys that previous versions accepted and that we now reject. Checked
    // before the generic unknown-key error so an outdated input file gets a message
    // naming the option that went away instead of "Unknown key found".
    void checkRemovedKey(const std::string& key, const std::string& errorPath) {
        static const std::set<std::string> removedKeys = {
            "PARAM SWEEP SUBMENU",
            "Parameter sweep (T/F)",
            "Run Monte Carlo (T/F)",
            "Num Monte Carlo runs (int)",
        };
        if (removedKeys.contains(key)) {
            throw std::runtime_error("Removed option found: '" + errorPath + "'. Delete it from your input file. " + ONE_RUN_PER_PROCESS_MESSAGE);
        }
    }

    void validateYamlKeys(const YAML::Node& defaultNode, const YAML::Node& userNode, const std::string& currentPath = "") {
        // Values in user yaml will replace the default node so we ensure that the types are compatible (value vs map)
        // If the userNode is a value, defaultNode should also be a value
        if (!userNode.IsMap()) {
            // A missing or null user value means "no override": mergeYamlNodes keeps the default.
            if (defaultNode.IsMap() && userNode.IsDefined() && !userNode.IsNull()) {
                // Edge case: top level of YAML is just a value (e.g empty file with only a value)
                if (currentPath.empty()) {
                    throw std::runtime_error("The document root is a value in provided YAML but a map in the default input.yaml (should be a set of menus).");
                }
                throw std::runtime_error("Invalid key: '" + currentPath + "' is a value in provided YAML but a map in the default input.yaml (should be a submenu).");
            }
            // User node is not a map, and default is not either so they are compatible
            // Because the user node is not a map it has no keys to check so we are done
            // validating this branch.
            return;
        }

        // Second compatibility check that is the reciprocal of ^
        // If the userNode is a map, then defaultNode should also be a map,
        // otherwise the user is trying to add a structure that doesn't exist
        if (!defaultNode.IsMap()) {
            throw std::runtime_error("Invalid key: '" + currentPath + "' is a map in provided YAML but not in the default input.yaml (should be a value).");
        }

        // Check this level before looking at the keys: getYamlKeys() returns a set,
        // which hides a repeated key.
        checkNoDuplicateKeys(userNode, "the input file", currentPath);

        auto defaultKeys = getYamlKeys(defaultNode);
        auto userKeys = getYamlKeys(userNode);

        for (const auto& key : userKeys) {
            if (!defaultKeys.contains(key)) {
                // The key from the user's YAML does not exist in the default YAML.
                std::string errorPath = buildKeyPath(currentPath, key);
                if (checkDeprecatedKey(key, errorPath)) {
                    continue;
                }
                checkRemovedKey(key, errorPath);
                throw std::runtime_error("Unknown key found: '" + errorPath + "'");
            }

            // Recurse into every map/value to check their validity
            const YAML::Node nextUserNode = userNode[key];
            const YAML::Node nextDefaultNode = defaultNode[key];
            std::string nextPath = buildKeyPath(currentPath, key);

            validateYamlKeys(nextDefaultNode, nextUserNode, nextPath);
        }
    }

    YAML::Node mergeYamlInputFiles(const vector<string>& filenames){
        YAML::Node defaultData = YAML::Load(default_input);
        YAML::Node mergedData = YAML::Load(default_input);

        // The defaults are compiled in and do not depend on the input files, so
        // check them once, before the loop. Errors here should not happen as the
        // default should not be modified and hopefully not be distributed with
        // errors in it, but this is cheap to verify.
        try {
            checkNoDuplicateKeysRecursive(defaultData, "the default input.yaml");
            checkDefaultPaths(defaultData);
        } catch (const std::runtime_error& e) {
            throw std::runtime_error(std::string(e.what()) + " This is not a problem with your input file: check the state of your APCEMM repository.");
        }

        for (auto filename: filenames) {
            YAML::Node userData = YAML::LoadFile(filename);
        
            // Validate the user's YAML file against the default structure
            try {
                validateYamlKeys(defaultData, userData);
            } catch (const std::runtime_error& e) {
                throw std::runtime_error("Invalid field in YAML input file '" + filename + "': " + e.what());
            }
            // Resolve relative paths against the location of the current input file before merging
            resolvePathsInPlace(userData, std::filesystem::path(filename).parent_path());
            mergedData = mergeYamlNodes(mergedData, userData);
        }

        return mergedData;
    }

    void populateInput(OptInput& input, Input& scenario, const YAML::Node& mergedData){
        try {
            readSimMenu(input, mergedData["SIMULATION MENU"]);
        }
        catch (const std::exception& e) {
            throw std::runtime_error("Something went wrong in reading the SIMULATION MENU! Please double-check your input file with the reference in Code.v05-00/defaults/input.yaml\n  Exception: " + std::string(e.what()));
        }

        try {
            readParamMenu(scenario, mergedData["PARAMETER MENU"]);
        }
        catch (const std::exception& e) {
            throw std::runtime_error("Something went wrong in reading the PARAMETER MENU! Please double-check your input file with the reference in Code.v05-00/defaults/input.yaml\n  Exception: " + std::string(e.what()));
        }

        // Inputs have successfully been parsed, now check that they seem physical
        scenario.checkInputValidity();

        try {
            readTransportMenu(input, mergedData["TRANSPORT MENU"]);
        }
        catch (const std::exception& e) {
            throw std::runtime_error("Something went wrong in reading the TRANSPORT MENU! Please double-check your input file with the reference in Code.v05-00/defaults/input.yaml\n  Exception: " + std::string(e.what()));
        }
        
        try {
            readChemMenu(input, mergedData["CHEMISTRY MENU"]);
        }
        catch (const std::exception& e) {
            throw std::runtime_error("Something went wrong in reading the CHEMISTRY MENU! Please double-check your input file with the reference in Code.v05-00/defaults/input.yaml\n  Exception: " + std::string(e.what()));
        }

        try {
            readAeroMenu(input, mergedData["AEROSOL MENU"]);  
        }
        catch (const std::exception& e) {
            throw std::runtime_error("Something went wrong in reading the AEROSOL MENU! Please double-check your input file with the reference in Code.v05-00/defaults/input.yaml\n  Exception: " + std::string(e.what()));
        }

        try {
            readMetMenu(input, mergedData["METEOROLOGY MENU"]);
        }
        catch (const std::exception& e) {
            throw std::runtime_error("Something went wrong in reading the METEOROLOGY MENU! Please double-check your input file with the reference in Code.v05-00/defaults/input.yaml\n  Exception: " + std::string(e.what()));
        }

        try {
            readDiagMenu(input, mergedData["DIAGNOSTIC MENU"]);
        }
        catch (const std::exception& e) {
            throw std::runtime_error("Something went wrong in reading the DIAGNOSTIC MENU! Please double-check your input file with the reference in Code.v05-00/defaults/input.yaml\n  Exception: " + std::string(e.what()));
        }

        try {
            readAdvancedMenu(input, mergedData["ADVANCED OPTIONS MENU"]);
        }
        catch (const std::exception& e) {
            throw std::runtime_error("Something went wrong in reading the ADVANCED OPTIONS MENU! Please double-check your input file with the reference in Code.v05-00/defaults/input.yaml\n  Exception: " + std::string(e.what()));
        }
    }

    // Overwrite both seed fields with the actual seed used so the saved merged input file
    // can exactly reproduce the run.
    void recordEffectiveSeed(YAML::Node& node, unsigned int seed){
        const string submenuKey = "RANDOM NUMBER GENERATION SUBMENU";
        const string forceKey = "Force seed value (T/F)";
        const string seedKey = "Seed value (positive int)";

        // Double check that these keys exist otherwise we'd be creating
        // new nodes instead of updating in place. const[] look-up
        // does not create a new node if it does not exist
        // Only happens if the RNG submenu is changed and keys are not updated here
        const YAML::Node& readOnly = node;
        const YAML::Node existing = readOnly["SIMULATION MENU"][submenuKey];
        if (!existing.IsDefined() || !existing.IsMap() || !existing[forceKey].IsDefined() || !existing[seedKey].IsDefined()){
            throw std::runtime_error("Cannot record the effective seed: 'SIMULATION MENU -> "
                                     + submenuKey + "' with keys '" + forceKey + "' and '"
                                     + seedKey + "' is missing from the merged input.");
        }

        // Update seed values here
        YAML::Node seedSubmenu = node["SIMULATION MENU"][submenuKey];
        seedSubmenu[forceKey] = "T";
        seedSubmenu[seedKey] = seed;
    }

    // Output dir must exist before calling this
    void writeYaml(const YAML::Node& node, const std::filesystem::path& outputDir, const string& filename){
        const std::filesystem::path fullPath = outputDir / filename;
        std::ofstream out(fullPath);
        if (!out) {
            throw std::runtime_error("Could not open '" + fullPath.string() + "' for writing merged YAML.");
        }
        out << node << std::endl;
        if (!out) {
            throw std::runtime_error("Could not write merged YAML at '" + fullPath.string() + "'.");
        }
    }

    void readSimMenu(OptInput& input, const YAML::Node& simNode){
        input.SIMULATION_OMP_NUM_THREADS = parseIntString(simNode["OpenMP Num Threads (positive int)"].as<string>(), "OpenMP Num Threads (positive int)");
        if(input.SIMULATION_OMP_NUM_THREADS < 1){
            throw std::invalid_argument("OpenMP Num Threads (under SIMULATION MENU) cannot be less than 1!");
        }

        #ifdef DEBUG
            // In DEBUG, force single threaded code to avoid non-reproducibility due to multithreading
            std::cout << "Compiled in DEBUG mode: setting OMP_NUM_THREADS = 1 to force single threading" << std::endl;
            input.SIMULATION_OMP_NUM_THREADS = 1;
        #endif

        std::cout << ">>> APCEMM OpenMP Num Threads is set to " << input.SIMULATION_OMP_NUM_THREADS << std::endl;
        if (input.SIMULATION_OMP_NUM_THREADS > 4){
            std::cout << ">>> Hint: APCEMM performance does not scale past 4 threads..." << std::endl;
        }

        YAML::Node outputSubmenu = simNode["OUTPUT SUBMENU"];
        std::string outputFolder = readPath(outputSubmenu, "Output folder (string)");
        // Ensure that outputFolder has a real value set by the user
        if (isPathSentinel(outputFolder)) {
            throw std::invalid_argument("No output folder set: 'Output folder (string)' under SIMULATION MENU -> "
                                        "OUTPUT SUBMENU is '" + outputFolder + "'. APCEMM has no default output "
                                        "folder, so set one in your input file. A relative path resolves against "
                                        "the directory of the file that sets it.");
        }
        // Ensure path to save directory is terminated by "/"
        if ( outputFolder.back() != '/' ) {outputFolder = outputFolder + "/";}
        input.SIMULATION_OUTPUT_FOLDER = outputFolder;
        input.SIMULATION_OVERWRITE = parseBoolString(outputSubmenu["Overwrite if folder exists (T/F)"].as<string>(), "Overwrite if folder exists (T/F)");
        input.SIMULATION_THREADED_FFT = parseBoolString(simNode["Use threaded FFT (T/F)"].as<string>(), "Use threaded FFT (T/F)");

        YAML::Node fftwWisdomSubmenu = simNode["FFTW WISDOM SUBMENU"];
        input.SIMULATION_USE_FFTW_WISDOM = parseBoolString(fftwWisdomSubmenu["Use FFTW WISDOM (T/F)"].as<string>(), "Use FFTW WISDOM (T/F)");
        input.SIMULATION_DIRECTORY_W_WRITE_PERMISSION = readPath(fftwWisdomSubmenu, "Dir w/ write permission (string)");
        input.SIMULATION_INPUT_BACKG_COND = readPath(simNode, "Input background condition (string)");
        input.SIMULATION_INPUT_ENG_EI = readPath(simNode, "Input engine emissions (string)");

        YAML::Node saveForwardSubmenu = simNode["SAVE FORWARD RESULTS SUBMENU"];
        input.SIMULATION_SAVE_FORWARD = parseBoolString(saveForwardSubmenu["Save forward results (T/F)"].as<string>(), "Save forward results (T/F)");
        input.SIMULATION_FORWARD_FILENAME = saveForwardSubmenu["netCDF filename format (string)"].as<string>();

        YAML::Node adjointSubmenu = simNode["ADJOINT OPTIMIZATION SUBMENU"];
        input.SIMULATION_ADJOINT = parseBoolString(adjointSubmenu["Turn on adjoint optim. (T/F)"].as<string>(), "Turn on adjoint optim. (T/F)");
        input.SIMULATION_ADJOINT_FILENAME = adjointSubmenu["netCDF filename format (string)"].as<string>();

        YAML::Node boxModelSubmenu = simNode["BOX MODEL SUBMENU"];
        input.SIMULATION_BOXMODEL = parseBoolString(boxModelSubmenu["Run box model (T/F)"].as<string>(), "Run box model (T/F)");
        input.SIMULATION_BOX_FILENAME = boxModelSubmenu["netCDF filename format (string)"].as<string>();

        YAML::Node seedSubmenu = simNode["RANDOM NUMBER GENERATION SUBMENU"];
        input.SIMULATION_FORCE_SEED = parseBoolString(seedSubmenu["Force seed value (T/F)"].as<string>(), "Force seed value (T/F)");
        input.SIMULATION_SEED_VALUE = parseScalarUIntParam(seedSubmenu["Seed value (positive int)"].as<string>(), "Seed value (positive int)");

        string epm =
            simNode["EPM type (original/external/new)"].as<string>();
        if (iequals(epm, "original")) {
            input.SIMULATION_EPM_TYPE = epm_type::EPM_ORIGINAL;
          } else if (iequals(epm, "external")) {
            input.SIMULATION_EPM_TYPE = epm_type::EPM_EXTERNAL;
          } else if (iequals(epm, "new")) {
            input.SIMULATION_EPM_TYPE = epm_type::EPM_NEW_PHYSICS;
        } else {
            throw std::invalid_argument("Invalid EPM type specified in SIMULATION MENU: " + epm);
        }

        input.SIMULATION_EXTERNAL_EPM_NETCDF_FILENAME = readPath(simNode, "External EPM NetCDF file");
    }
    void readParamMenu(Input& scenario, const YAML::Node& paramNode){

        scenario.set_simulationTime(parseScalarParam(paramNode["Plume Process [hr] (double)"].as<string>(), "Plume Process [hr] (double)"));

        YAML::Node metParamSubmenu = paramNode["METEOROLOGICAL PARAMETERS SUBMENU"];
        //convert hPa to Pa because the solver uses Pa as the default unit
        scenario.set_pressure_Pa(100 * parseScalarParam(metParamSubmenu["Pressure [hPa] (double)"].as<string>(), "Pressure [hPa] (double)"));
        scenario.set_horizDiff(parseScalarParam(metParamSubmenu["Horiz. diff. coeff. [m^2/s] (double)"].as<string>(), "Horiz. diff. coeff. [m^2/s] (double)"));
        scenario.set_vertiDiff(parseScalarParam(metParamSubmenu["Verti. diff. [m^2/s] (double)"].as<string>(), "Verti. diff. [m^2/s] (double)"));
        scenario.set_nBV(parseScalarParam(metParamSubmenu["Brunt-Vaisala Frequency [s^-1] (double)"].as<string>(), "Brunt-Vaisala Frequency [s^-1] (double)"));

        YAML::Node locTimeSubmenu = paramNode["LOCATION AND TIME SUBMENU"];
        scenario.set_longitude_deg(parseScalarParam(locTimeSubmenu["LON [deg] (double)"].as<string>(), "LON [deg] (double)"));
        scenario.set_latitude_deg(parseScalarParam(locTimeSubmenu["LAT [deg] (double)"].as<string>(), "LAT [deg] (double)"));
        scenario.set_emissionDOY(parseScalarUIntParam(locTimeSubmenu["Emission day [1-365] (int)"].as<string>(), "Emission day [1-365] (int)"));
        scenario.set_emissionTime(parseScalarParam(locTimeSubmenu["Emission time [hr] (double)"].as<string>(), "Emission time [hr] (double)"));

        YAML::Node backMixRatioSubmenu = paramNode["BACKGROUND MIXING RATIOS SUBMENU"];
        scenario.set_backgNOx(parseScalarParam(backMixRatioSubmenu["NOx [ppt] (double)"].as<string>(), "NOx [ppt] (double)"));
        scenario.set_backgHNO3(parseScalarParam(backMixRatioSubmenu["HNO3 [ppt] (double)"].as<string>(), "HNO3 [ppt] (double)"));
        scenario.set_backgO3(parseScalarParam(backMixRatioSubmenu["O3 [ppb] (double)"].as<string>(), "O3 [ppb] (double)"));
        scenario.set_backgCO(parseScalarParam(backMixRatioSubmenu["CO [ppb] (double)"].as<string>(), "CO [ppb] (double)"));
        scenario.set_backgCH4(parseScalarParam(backMixRatioSubmenu["CH4 [ppm] (double)"].as<string>(), "CH4 [ppm] (double)"));
        scenario.set_backgSO2(parseScalarParam(backMixRatioSubmenu["SO2 [ppt] (double)"].as<string>(), "SO2 [ppt] (double)"));

        YAML::Node eiSubmenu = paramNode["EMISSION INDICES SUBMENU"];
        scenario.set_EI_NOx(parseScalarParam(eiSubmenu["NOx [g(NO2)/kg_fuel] (double)"].as<string>(), "NOx [g(NO2)/kg_fuel] (double)"));
        scenario.set_EI_CO(parseScalarParam(eiSubmenu["CO [g/kg_fuel] (double)"].as<string>(), "CO [g/kg_fuel] (double)"));
        scenario.set_EI_HC(parseScalarParam(eiSubmenu["UHC [g/kg_fuel] (double)"].as<string>(), "UHC [g/kg_fuel] (double)"));
        scenario.set_EI_SO2(parseScalarParam(eiSubmenu["SO2 [g/kg_fuel] (double)"].as<string>(), "SO2 [g/kg_fuel] (double)"));
        //Convert % to ratio
        scenario.set_EI_SO2TOSO4(parseScalarParam(eiSubmenu["SO2 to SO4 conv [%] (double)"].as<string>(), "SO2 to SO4 conv [%] (double)") / 100.0);
        scenario.set_EI_Soot(parseScalarParam(eiSubmenu["Soot [g/kg_fuel] (double)"].as<string>(), "Soot [g/kg_fuel] (double)"));

        scenario.set_sootRad(parseScalarParam(paramNode["Soot Radius [m] (double)"].as<string>(), "Soot Radius [m] (double)"));
        scenario.set_fuelFlow(parseScalarParam(paramNode["Total fuel flow [kg/s] (double)"].as<string>(), "Total fuel flow [kg/s] (double)"));
        scenario.set_aircraftMass(parseScalarParam(paramNode["Aircraft mass [kg] (double)"].as<string>(), "Aircraft mass [kg] (double)"));
        scenario.set_flightSpeed(parseScalarParam(paramNode["Flight speed [m/s] (double)"].as<string>(), "Flight speed [m/s] (double)"));
        scenario.set_numEngines(parseScalarParam(paramNode["Num. of engines [2/4] (int)"].as<string>(), "Num. of engines [2/4] (int)"));
        scenario.set_wingspan(parseScalarParam(paramNode["Wingspan [m] (double)"].as<string>(), "Wingspan [m] (double)"));
        scenario.set_coreExitTemp(parseScalarParam(paramNode["Core exit temp. [K] (double)"].as<string>(), "Core exit temp. [K] (double)"));
        scenario.set_bypassArea(parseScalarParam(paramNode["Exit bypass area [m^2] (double)"].as<string>(), "Exit bypass area [m^2] (double)"));
    }
    void readTransportMenu(OptInput& input, const YAML::Node& transportNode){
        input.TRANSPORT_TRANSPORT = parseBoolString(transportNode["Turn on Transport (T/F)"].as<string>(), "Turn on Transport (T/F)");
        input.TRANSPORT_FILL = parseBoolString(transportNode["Fill Negative Values (T/F)"].as<string>(), "Fill Negative Values (T/F)");
        input.TRANSPORT_TIMESTEP = parseDoubleString(transportNode["Transport Timestep [min] (double)"].as<string>(), "Transport Timestep [min] (double)");

        YAML::Node updraftSubmenu = transportNode["PLUME UPDRAFT SUBMENU"];
        input.TRANSPORT_UPDRAFT = parseBoolString(updraftSubmenu["Turn on plume updraft (T/F)"].as<string>(), "Turn on plume updraft (T/F)");
        input.TRANSPORT_UPDRAFT_TIMESCALE = parseDoubleString(updraftSubmenu["Updraft timescale [s] (double)"].as<string>(), "Updraft timescale [s] (double)");
        input.TRANSPORT_UPDRAFT_VELOCITY = parseDoubleString(updraftSubmenu["Updraft veloc. [cm/s] (double)"].as<string>(), "Updraft veloc. [cm/s] (double)");
    }
    void readChemMenu(OptInput& input, const YAML::Node& chemNode){
        input.CHEMISTRY_CHEMISTRY = parseBoolString(chemNode["Turn on Chemistry (T/F)"].as<string>(), "Turn on Chemistry (T/F)");
        input.CHEMISTRY_HETCHEM = parseBoolString(chemNode["Perform hetero. chem. (T/F)"].as<string>(), "Perform hetero. chem. (T/F)");
        input.CHEMISTRY_JRATE_FOLDER = readPath(chemNode, "Photolysis rates folder (string)");
    }
    void readAeroMenu(OptInput& input, const YAML::Node& aeroNode){
        input.AEROSOL_GRAVSETTLING = parseBoolString(aeroNode["Turn on grav. settling (T/F)"].as<string>(), "Turn on grav. settling (T/F)");
        input.AEROSOL_COAGULATION_SOLID = parseBoolString(aeroNode["Turn on solid coagulation (T/F)"].as<string>(), "Turn on solid coagulation (T/F)");
        input.AEROSOL_COAGULATION_LIQUID = parseBoolString(aeroNode["Turn on liquid coagulation (T/F)"].as<string>(), "Turn on liquid coagulation (T/F)");
        input.AEROSOL_ICE_GROWTH = parseBoolString(aeroNode["Turn on ice growth (T/F)"].as<string>(), "Turn on ice growth (T/F)");
        input.AEROSOL_ICE_GROWTH_TIMESTEP = parseDoubleString(aeroNode["Ice growth timestep [min] (double)"].as<string>(), "Ice growth timestep [min] (double)");
    }
    void readMetMenu(OptInput& input, const YAML::Node& metNode){
        YAML::Node metInputSubmenu = metNode["METEOROLOGICAL INPUT SUBMENU"];
        input.MET_LOADMET = parseBoolString(metInputSubmenu["Use met. input (T/F)"].as<string>(), "Use met. input (T/F)");
        input.MET_FILENAME = readPath(metInputSubmenu, "Met input file path (string)");
        input.MET_DT = parseDoubleString(metInputSubmenu["Time series data timestep [hr] (double)"].as<string>(), "Time series data timestep [hr] (double)");
        input.MET_LOADTEMP = parseBoolString(metInputSubmenu["Init temp. from met. (T/F)"].as<string>(), "Init temp. from met. (T/F)");
        input.MET_TEMPTIMESERIES = parseBoolString(metInputSubmenu["Temp. time series input (T/F)"].as<string>(), "Temp. time series input (T/F)");
        input.MET_INTERPTEMPDATA = parseBoolString(metInputSubmenu["Interpolate temp. met. data (T/F)"].as<string>(), "Interpolate temp. met. data (T/F)");
        input.MET_LOADRH = parseBoolString(metInputSubmenu["Init RH from met. (T/F)"].as<string>(), "Init RH from met. (T/F)");
        input.MET_RHTIMESERIES = parseBoolString(metInputSubmenu["RH time series input (T/F)"].as<string>(), "RH time series input (T/F)");
        input.MET_INTERPRHDATA = parseBoolString(metInputSubmenu["Interpolate RH met. data (T/F)"].as<string>(), "Interpolate RH met. data (T/F)");
        input.MET_LOADSHEAR = parseBoolString(metInputSubmenu["Init wind shear from met. (T/F)"].as<string>(), "Init wind shear from met. (T/F)");
        input.MET_SHEARTIMESERIES = parseBoolString(metInputSubmenu["Wind shear time series input (T/F)"].as<string>(), "Wind shear time series input (T/F)");
        input.MET_INTERPSHEARDATA = parseBoolString(metInputSubmenu["Interpolate shear met. data (T/F)"].as<string>(), "Interpolate shear met. data (T/F)");
        input.MET_LOADVERTVELOC = parseBoolString(metInputSubmenu["Init vert. veloc. from met. data (T/F)"].as<string>(), "Init vert. veloc. from met. data (T/F)");
        input.MET_VERTVELOCTIMESERIES = parseBoolString(metInputSubmenu["Vert. veloc. time series input (T/F)"].as<string>(), "Vert. veloc. time series input (T/F)");
        input.MET_INTERPVERTVELOC = parseBoolString(metInputSubmenu["Interpolate vert. veloc. met. data (T/F)"].as<string>(), "Interpolate vert. veloc. met. data (T/F)");
        
        YAML::Node tempPerturbMenu = metNode["TEMPERATURE PERTURBATION SUBMENU"];
        input.MET_ENABLE_TEMP_PERTURB = parseBoolString( tempPerturbMenu["Enable Temp. Pert. (T/F)"].as<string>(), "Enable Temp. Pert. (T/F)" );
        input.MET_TEMP_PERTURB_AMPLITUDE = parseDoubleString( tempPerturbMenu["Temp. Perturb. Amplitude (double)"].as<string>(), "Temp. Perturb. Amplitude (double)" );
    }
    void readDiagMenu(OptInput& input, const YAML::Node& diagNode){
        input.DIAG_FILENAME = diagNode["netCDF filename format (string)"].as<string>();

        YAML::Node specTsSubmenu = diagNode["SPECIES TIMESERIES SUBMENU"];
        input.TS_SPEC = parseBoolString(specTsSubmenu["Save species timeseries (T/F)"].as<string>(), "Save species timeseries (T/F)");
        input.TS_FILENAME = specTsSubmenu["Inst timeseries file (string)"].as<string>();
        input.TS_SPECIES = parseVectorIntString(specTsSubmenu["Species indices to include (list of ints)"].as<string>(), "Species indices to include (list of ints)");
        input.TS_FREQ = parseDoubleString(specTsSubmenu["Save frequency [min] (double)"].as<string>(), "Save frequency [min] (double)");

        YAML::Node aeroTsSubmenu = diagNode["AEROSOL TIMESERIES SUBMENU"];
        input.TS_AERO = parseBoolString(aeroTsSubmenu["Save aerosol timeseries (T/F)"].as<string>(), "Save aerosol timeseries (T/F)");
        input.TS_AERO_FILENAME = aeroTsSubmenu["Inst timeseries file (string)"].as<string>();
        input.TS_AEROSOL = parseVectorIntString(aeroTsSubmenu["Aerosol indices to include (list of ints)"].as<string>(), "Aerosol indices to include (list of ints)");
        input.TS_AERO_FREQ = parseDoubleString(aeroTsSubmenu["Save frequency [min] (double)"].as<string>(), "Save frequency [min] (double)");

        YAML::Node plSubmenu = diagNode["PRODUCTION & LOSS SUBMENU"];
        input.PL_PL = parseBoolString(plSubmenu["Turn on P/L diag (T/F)"].as<string>(), "Turn on P/L diag (T/F)");
        input.PL_O3 = parseBoolString(plSubmenu["Save O3 P/L (T/F)"].as<string>(), "Save O3 P/L (T/F)");
    }

    void readAdvancedMenu(OptInput& input, const YAML::Node& advancedNode) {
        YAML::Node gridSubmenu = advancedNode["GRID SUBMENU"];
        input.ADV_GRID_NX = parseUIntString(gridSubmenu["NX (positive int)"].as<string>(), "NX (positive int)");
        input.ADV_GRID_NY = parseUIntString(gridSubmenu["NY (positive int)"].as<string>(), "NY (positive int)");
        input.ADV_GRID_XLIM_RIGHT = parseDoubleString(gridSubmenu["XLIM_RIGHT (positive double)"].as<string>(), "XLIM_RIGHT (positive double)");
        input.ADV_GRID_XLIM_LEFT = parseDoubleString(gridSubmenu["XLIM_LEFT (positive double)"].as<string>(), "XLIM_LEFT (positive double)");
        input.ADV_GRID_YLIM_UP = parseDoubleString(gridSubmenu["YLIM_UP (positive double)"].as<string>(), "YLIM_UP (positive double)");
        input.ADV_GRID_YLIM_DOWN = parseDoubleString(gridSubmenu["YLIM_DOWN (positive double)"].as<string>(), "YLIM_DOWN (positive double)");
        
        YAML::Node csizeSubmenu = advancedNode["INITIAL CONTRAIL SIZE SUBMENU"];
        input.ADV_CSIZE_DEPTH_BASE = parseDoubleString(csizeSubmenu["Base Contrail Depth [m] (double)"].as<string>(), "Base Contrail Depth [m] (double)");
        input.ADV_CSIZE_DEPTH_SCALING_FACTOR = parseDoubleString(csizeSubmenu["Contrail Depth Scaling Factor [-] (double)"].as<string>(), "Contrail Depth Scaling Factor [-] (double)");
        input.ADV_CSIZE_WIDTH_BASE = parseDoubleString(csizeSubmenu["Base Contrail Width [m] (double)"].as<string>(), "Base Contrail Width [m] (double)");
        input.ADV_CSIZE_WIDTH_SCALING_FACTOR = parseDoubleString(csizeSubmenu["Contrail Width Scaling Factor [-] (double)"].as<string>(), "Contrail Width Scaling Factor [-] (double)");

        input.ADV_AMBIENT_LAPSERATE = parseDoubleString(advancedNode["Ambient Lapse Rate [K/km] (double)"].as<string>(), "Ambient Lapse Rate [K/km] (double)");
        input.ADV_TROPOPAUSE_PRESSURE = parseDoubleString(advancedNode["Tropopause Pressure [Pa] (double)"].as<string>(), "Tropopause Pressure [Pa] (double)");

        if (input.ADV_GRID_XLIM_LEFT < 0 ||
            input.ADV_GRID_XLIM_RIGHT < 0 ||
            input.ADV_GRID_YLIM_UP < 0 ||
            input.ADV_GRID_YLIM_DOWN < 0) {
            throw std::invalid_argument("No values in GRID SUBMENU can be less than zero!");
        }

        YAML::Node earlyPlumeSubmenu = advancedNode["EARLY PLUME SUBMENU"];
        input.ADV_EP_N_REF = parseDoubleString(earlyPlumeSubmenu["Reference ice crystal count [#/m] (double)"].as<string>(), "Reference ice crystal count [#/m] (double)");
        input.ADV_EP_WINGSPAN_REF = parseDoubleString(earlyPlumeSubmenu["Reference wingspan [m] (double)"].as<string>(), "Reference wingspan [m] (double)");
        input.ADV_EP_N_POSTJET_OVERRIDE = parseBoolString(earlyPlumeSubmenu["Override post-jet ice crystal count (T/F)"].as<string>(), "Override post-jet ice crystal count (T/F)");
        input.ADV_EP_N_POSTJET = parseDoubleString(earlyPlumeSubmenu["Post-jet ice crystal count [#/m] (double)"].as<string>(), "Post-jet ice crystal count [#/m] (double)");
        input.ADV_SAVE_PSD_GRID = parseBoolString(advancedNode["Save gridded particle size distribution (T/F)"].as<string>(), "Save gridded particle size distribution (T/F)");
    }

    // Wrapper around parseDoubleString to have a nice rejection message for
    // deprecated sweep style inputs (e.g. "200 220 240" and "200:20:240")
    double parseScalarParam(const string paramString, const string paramLocation){
        const string s = trim(paramString);
        if(s.find(':') != string::npos || split(s, " ").size() > 1){
            throw std::invalid_argument("Several values given at " + paramLocation + ". " + ONE_RUN_PER_PROCESS_MESSAGE);
        }
        return parseDoubleString(s, paramLocation);
    }

    // Same as parseScalarParam but returns a UInt
    UInt parseScalarUIntParam(const string paramString, const string paramLocation){
        double valueDouble = parseScalarParam(paramString, paramLocation);
        if (valueDouble < 0 || valueDouble > std::numeric_limits<UInt>::max()) {
            throw std::invalid_argument("Value out of range [0, maxUInt] at " + paramLocation);
        }
        if(std::fmod(valueDouble, 1.0) > 1e-40) {
            throw (std::invalid_argument("Decimals not allowed in int inputs at " + paramLocation + "!"));
        }
        return static_cast<UInt>(valueDouble);
    }

    // Space-separated list, used by the species and aerosol index lists of the
    // DIAGNOSTIC MENU. Those are lists of outputs so we still support it for now.
    Vector_1D parseVectorDoubleString(const string paramString, const string paramLocation){
        const vector<string> tokens = split(trim(paramString), " ");
        Vector_1D values;
        for (const string& token: tokens){
            values.push_back(parseDoubleString(token, paramLocation));
        }
        return values;
    }

    vector<string> split(const string str, const string delimiter){
        if (delimiter.length() == 0){
            throw std::invalid_argument("In YamlInputReader::split: Delimiter length cannot be zero!");
        }
        string s = str;
        vector<string> tokens;
        size_t pos = 0;
        std::string token;
        while ((pos = s.find(delimiter)) != std::string::npos) {
            token = s.substr(0, pos);
            //Ignores tokens with length 0
            if(token.length() > 0) tokens.push_back(token);
            s.erase(0, pos + delimiter.length());
        }
        if(s.length() > 0) tokens.push_back(s);
        return tokens;
    }

    vector<int> parseVectorIntString(const string paramString, const string paramLocation){
        Vector_1D vec1d = parseVectorDoubleString(paramString, paramLocation);
        vector<int> vecint;
        for (double d: vec1d){
            if(std::fmod(d,1.0) > 1e-40) throw (std::invalid_argument("Decimals not allowed in int inputs at " + paramLocation + "!"));
            vecint.push_back((int)(d));
        }
        return vecint;
    }
}
