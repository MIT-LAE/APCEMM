#include "YamlInputReader/YamlInputReader.hpp"

using std::cout;
using std::endl;
namespace YamlInputReader{
    void readYamlInputFile(OptInput& input, string filename){
        INPUT_FILE_PATH = std::filesystem::path(filename);
        YAML::Node data = YAML::LoadFile(filename);

        try {
            readSimMenu(input, data["SIMULATION MENU"]);
        }
        catch (...) {
            std::cout << "Something went wrong in reading the SIMULATION MENU! Please double-check your input file with the reference in SampleRunDir!";
            exit(1);
        }

        try {
            readParamMenu(input, data["PARAMETER MENU"]);
        }
        catch (...) {
            std::cout << "Something went wrong in reading the PARAMETER MENU! Please double-check your input file with the reference in SampleRunDir!";
            exit(1);
        }

        try {
            readTransportMenu(input, data["TRANSPORT MENU"]);
        }
        catch (...) {
            std::cout << "Something went wrong in reading the TRANSPORT MENU! Please double-check your input file with the reference in SampleRunDir!";
            exit(1);
        }
        
        try {
            readChemMenu(input, data["CHEMISTRY MENU"]);
        }
        catch (...) {
            std::cout << "Something went wrong in reading the CHEMISTRY MENU! Please double-check your input file with the reference in SampleRunDir!";
            exit(1);
        }

        try {
            readAeroMenu(input, data["AEROSOL MENU"]);  
        }
        catch (...) {
            std::cout << "Something went wrong in reading the AEROSOL MENU! Please double-check your input file with the reference in SampleRunDir!";
            exit(1);
        }

        try {
            readMetMenu(input, data["METEOROLOGY MENU"]);
        }
        catch (const std::invalid_argument& e) {
            std::cerr << "ERROR: " << e.what() << std::endl;
            exit(1);
        }
        catch (...) {
            std::cout << "Something went wrong in reading the METEOROLOGY MENU! Please double-check your input file with the reference in SampleRunDir!";
            exit(1);
        }

        try {
            readDiagMenu(input, data["DIAGNOSTIC MENU"]);
        }
        catch (...) {
            std::cout << "Something went wrong in reading the DIAGNOSTIC MENU! Please double-check your input file with the reference in SampleRunDir!";
            exit(1);
        }

        try {
            readAdvancedMenu(input, data["ADVANCED OPTIONS MENU"]);
        }
        catch (...) {
            std::cout << "Something went wrong in reading the ADVANCED OPTIONS MENU! Please double-check your input file with the reference in SampleRunDir!";
            exit(1);
        }
    }
    void readSimMenu(OptInput& input, const YAML::Node& simNode){
        input.SIMULATION_OMP_NUM_THREADS = parseIntString(simNode["OpenMP Num Threads (positive int)"].as<string>(), "OpenMP Num Threads (positive int)");
        if(input.SIMULATION_OMP_NUM_THREADS < 1){
            throw std::invalid_argument("OpenMP Num Threads (under SIMULATION MENU) cannot be less than 1!");
        }

        YAML::Node paramSweepSubmenu = simNode["PARAM SWEEP SUBMENU"];
        input.SIMULATION_PARAMETER_SWEEP = parseBoolString(paramSweepSubmenu["Parameter sweep (T/F)"].as<string>(), "Parameter sweep (T/F");
        input.SIMULATION_MONTECARLO = parseBoolString(paramSweepSubmenu["Run Monte Carlo (T/F)"].as<string>(), "Run Monte Carlo (T/F)");
        input.SIMULATION_MCRUNS =  parseIntString(paramSweepSubmenu["Num Monte Carlo runs (int)"].as<string>(), "Num Monte Carlo runs (int)");

        YAML::Node outputSubmenu = simNode["OUTPUT SUBMENU"];
        input.SIMULATION_OUTPUT_FOLDER = parseFileSystemPath(outputSubmenu["Output folder (string)"].as<string>());
        input.SIMULATION_OVERWRITE = parseBoolString(outputSubmenu["Overwrite if folder exists (T/F)"].as<string>(), "Overwrite if folder exists (T/F)");
        input.SIMULATION_THREADED_FFT = parseBoolString(simNode["Use threaded FFT (T/F)"].as<string>(), "Use threaded FFT (T/F)");

        YAML::Node fftwWisdomSubmenu = simNode["FFTW WISDOM SUBMENU"];
        input.SIMULATION_USE_FFTW_WISDOM = parseBoolString(fftwWisdomSubmenu["Use FFTW WISDOM (T/F)"].as<string>(), "Use FFTW WISDOM (T/F)");
        input.SIMULATION_DIRECTORY_W_WRITE_PERMISSION = parseFileSystemPath(fftwWisdomSubmenu["Dir w/ write permission (string)"].as<string>());
        input.SIMULATION_INPUT_BACKG_COND = parseFileSystemPath(simNode["Input background condition (string)"].as<string>());
        input.SIMULATION_INPUT_ENG_EI = parseFileSystemPath(simNode["Input engine emissions (string)"].as<string>());

        YAML::Node saveForwardSubmenu = simNode["SAVE FORWARD RESULTS SUBMENU"];
        input.SIMULATION_SAVE_FORWARD = parseBoolString(saveForwardSubmenu["Save forward results (T/F)"].as<string>(), "Save forward results (T/F)");
        input.SIMULATION_FORWARD_FILENAME = saveForwardSubmenu["netCDF filename format (string)"].as<string>();

        YAML::Node adjointSubmenu = simNode["ADJOINT OPTIMIZATION SUBMENU"];
        input.SIMULATION_ADJOINT = parseBoolString(adjointSubmenu["Turn on adjoint optim. (T/F)"].as<string>(), "Turn on adjoint optim. (T/F)");
        input.SIMULATION_ADJOINT_FILENAME = adjointSubmenu["netCDF filename format (string)"].as<string>();

        YAML::Node boxModelSubmenu = simNode["BOX MODEL SUBMENU"];
        input.SIMULATION_BOXMODEL = parseBoolString(boxModelSubmenu["Run box model (T/F)"].as<string>(), "Run box model (T/F)");
        input.SIMULATION_BOX_FILENAME = boxModelSubmenu["netCDF filename format (string)"].as<string>();

        if(input.SIMULATION_PARAMETER_SWEEP == input.SIMULATION_MONTECARLO){
            throw std::invalid_argument("In Simulation Menu: Parameter sweep and Monte Carlo cannot have the same value!");
        }
    }
    void readParamMenu(OptInput& input, const YAML::Node& paramNode){

        input.PARAMETER_PARAM_MAP["PLUMEPROCESS"] = parseParamSweepInput(paramNode["Plume Process [hr] (double)"].as<string>(), "Plume Process [hr] (double)");

        YAML::Node metParamSubmenu = paramNode["METEOROLOGICAL PARAMETERS SUBMENU"];
        input.PARAMETER_PARAM_MAP["TEMPERATURE"] = parseParamSweepInput(metParamSubmenu["Temperature [K] (double)"].as<string>(), "Temperature [K] (double)");
        input.PARAMETER_PARAM_MAP["RHW"] = parseParamSweepInput(metParamSubmenu["R.Hum. wrt water [%] (double)"].as<string>(), "R.Hum. wrt water [%] (double)");
        input.PARAMETER_PARAM_MAP["PRESSURE"] = parseParamSweepInput(metParamSubmenu["Pressure [hPa] (double)"].as<string>(), "Pressure [hPa] (double)");
        input.PARAMETER_PARAM_MAP["DH"] = parseParamSweepInput(metParamSubmenu["Horiz. diff. coeff. [m^2/s] (double)"].as<string>(), "Horiz. diff. coeff. [m^2/s] (double)");
        input.PARAMETER_PARAM_MAP["DV"] = parseParamSweepInput(metParamSubmenu["Verti. diff. [m^2/s] (double)"].as<string>(), "Verti. diff. [m^2/s] (double)");
        input.PARAMETER_PARAM_MAP["SHEAR"] = parseParamSweepInput(metParamSubmenu["Wind shear [1/s] (double)"].as<string>(), "Wind shear [1/s] (double)");
        input.PARAMETER_PARAM_MAP["NBV"] = parseParamSweepInput(metParamSubmenu["Brunt-Vaisala Frequency [s^-1] (double)"].as<string>(), "Brunt-Vaisala Frequency [s^-1] (double)");

        YAML::Node locTimeSubmenu = paramNode["LOCATION AND TIME SUBMENU"];
        input.PARAMETER_PARAM_MAP["LONGITUDE"] = parseParamSweepInput(locTimeSubmenu["LON [deg] (double)"].as<string>(), "LAT [deg] (double)");
        input.PARAMETER_PARAM_MAP["LATITUDE"] = parseParamSweepInput(locTimeSubmenu["LAT [deg] (double)"].as<string>(), "LON [deg] (double)");
        input.PARAMETER_PARAM_MAP["EDAY"] = parseParamSweepInput(locTimeSubmenu["Emission day [1-365] (int)"].as<string>(), "Emission day [1-365] (int)");
        input.PARAMETER_PARAM_MAP["ETIME"] = parseParamSweepInput(locTimeSubmenu["Emission time [hr] (double)"].as<string>(), "Emission time [hr] (double)");
       
        YAML::Node backMixRatioSubmenu = paramNode["BACKGROUND MIXING RATIOS SUBMENU"];
        input.PARAMETER_PARAM_MAP["BACKG_NOX"] = parseParamSweepInput(backMixRatioSubmenu["NOx [ppt] (double)"].as<string>(), "NOx [ppt] (double)");
        input.PARAMETER_PARAM_MAP["BACKG_HNO3"] = parseParamSweepInput(backMixRatioSubmenu["HNO3 [ppt] (double)"].as<string>(), "HNO3 [ppt] (double)");
        input.PARAMETER_PARAM_MAP["BACKG_O3"] = parseParamSweepInput(backMixRatioSubmenu["O3 [ppb] (double)"].as<string>(), "O3 [ppb] (double)");
        input.PARAMETER_PARAM_MAP["BACKG_CO"] = parseParamSweepInput(backMixRatioSubmenu["CO [ppb] (double)"].as<string>(), "CO [ppb] (double)");
        input.PARAMETER_PARAM_MAP["BACKG_CH4"] = parseParamSweepInput(backMixRatioSubmenu["CH4 [ppm] (double)"].as<string>(), "CH4 [ppm] (double)");
        input.PARAMETER_PARAM_MAP["BACKG_SO2"] = parseParamSweepInput(backMixRatioSubmenu["SO2 [ppt] (double)"].as<string>(), "SO2 [ppt] (double)");

        YAML::Node eiSubmenu = paramNode["EMISSION INDICES SUBMENU"];
        input.PARAMETER_PARAM_MAP["EI_NOX"] = parseParamSweepInput(eiSubmenu["NOx [g(NO2)/kg_fuel] (double)"].as<string>(), "NOx [g(NO2)/kg_fuel] (double)");
        input.PARAMETER_PARAM_MAP["EI_CO"] = parseParamSweepInput(eiSubmenu["CO [g/kg_fuel] (double)"].as<string>(), "CO [g/kg_fuel] (double)");
        input.PARAMETER_PARAM_MAP["EI_UHC"] = parseParamSweepInput(eiSubmenu["UHC [g/kg_fuel] (double)"].as<string>(), "UHC [g/kg_fuel] (double)");
        input.PARAMETER_PARAM_MAP["EI_SO2"] = parseParamSweepInput(eiSubmenu["SO2 [g/kg_fuel] (double)"].as<string>(), "SO2 [g/kg_fuel] (double)");
        input.PARAMETER_PARAM_MAP["EI_SO2TOSO4"] =  parseParamSweepInput(eiSubmenu["SO2 to SO4 conv [%] (double)"].as<string>(), "SO2 to SO4 conv [%] (double)");
        input.PARAMETER_PARAM_MAP["EI_SOOT"] = parseParamSweepInput(eiSubmenu["Soot [g/kg_fuel] (double)"].as<string>(), "Soot [g/kg_fuel] (double)");
        
        input.PARAMETER_PARAM_MAP["EI_SOOTRAD"] = parseParamSweepInput(paramNode["Soot Radius [m] (double)"].as<string>(), "Soot Radius [m] (double)");
        input.PARAMETER_PARAM_MAP["FF"] = parseParamSweepInput(paramNode["Total fuel flow [kg/s] (double)"].as<string>(), "Total fuel flow [kg/s] (double)");
        input.PARAMETER_PARAM_MAP["AMASS"] = parseParamSweepInput(paramNode["Aircraft mass [kg] (double)"].as<string>(), "Aircraft mass [kg] (double)");
        input.PARAMETER_PARAM_MAP["FSPEED"] = parseParamSweepInput(paramNode["Flight speed [m/s] (double)"].as<string>(), "Flight speed [m/s] (double)");
        input.PARAMETER_PARAM_MAP["NUMENG"] = parseParamSweepInput(paramNode["Num. of engines [2/4] (int)"].as<string>(), " Num. of engines [2/4] (int)"); // Why is this a vector1d in the first place...
        input.PARAMETER_PARAM_MAP["WINGSPAN"] = parseParamSweepInput(paramNode["Wingspan [m] (double)"].as<string>(), "Wingspan [m] (double)");
        input.PARAMETER_PARAM_MAP["COREEXITTEMP"] = parseParamSweepInput(paramNode["Core exit temp. [K] (double)"].as<string>(), "Core exit temp. [K] (double)");
        input.PARAMETER_PARAM_MAP["BYPASSAREA"] = parseParamSweepInput(paramNode["Exit bypass area [m^2] (double)"].as<string>(), "Exit bypass area [m^2] (double)");
        
        //convert hPa to Pa because the solver uses Pa as the default unit
        for(double& i: input.PARAMETER_PARAM_MAP["PRESSURE"]){
            i *= 100;
        }
        //Convert % to ratio
        for(double& i : input.PARAMETER_PARAM_MAP["EI_SO2TOSO4"] ){
            i *= 1.0/100.0;
        }
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
        input.CHEMISTRY_TIMESTEP = parseDoubleString(chemNode["Chemistry Timestep [min] (double)"].as<string>(), "Chemistry Timestep [min] (double)");
        input.CHEMISTRY_JRATE_FOLDER = parseFileSystemPath(chemNode["Photolysis rates folder (string)"].as<string>());
    }
    void readAeroMenu(OptInput& input, const YAML::Node& aeroNode){
        input.AEROSOL_GRAVSETTLING = parseBoolString(aeroNode["Turn on grav. settling (T/F)"].as<string>(), "Turn on grav. settling (T/F)");
        input.AEROSOL_COAGULATION_SOLID = parseBoolString(aeroNode["Turn on solid coagulation (T/F)"].as<string>(), "Turn on solid coagulation (T/F)");
        input.AEROSOL_COAGULATION_LIQUID = parseBoolString(aeroNode["Turn on liquid coagulation (T/F)"].as<string>(), "Turn on liquid coagulation (T/F)");
        input.AEROSOL_COAGULATION_TIMESTEP = parseDoubleString(aeroNode["Coag. timestep [min] (double)"].as<string>(), "Coag. timestep [min] (double)");
        input.AEROSOL_ICE_GROWTH = parseBoolString(aeroNode["Turn on ice growth (T/F)"].as<string>(), "Turn on ice growth (T/F)");
        input.AEROSOL_ICE_GROWTH_TIMESTEP = parseDoubleString(aeroNode["Ice growth timestep [min] (double)"].as<string>(), "Ice growth timestep [min] (double)");
    }
    void readMetMenu(OptInput& input, const YAML::Node& metNode){
        YAML::Node metInputSubmenu = metNode["METEOROLOGICAL INPUT SUBMENU"];
        input.MET_LOADMET = parseBoolString(metInputSubmenu["Use met. input (T/F)"].as<string>(), "Use met. input (T/F)");
        input.MET_FILENAME = parseFileSystemPath(metInputSubmenu["Met input file path (string)"].as<string>());
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
        YAML::Node humidScalingOpts = metInputSubmenu["HUMIDITY SCALING OPTIONS"];
        input.MET_HUMIDSCAL_MODIFICATION_SCHEME = humidScalingOpts["Humidity modification scheme (none / constant / scaling)"].as<string>();
        input.MET_HUMIDSCAL_CONST_RHI = parseDoubleString(humidScalingOpts["Constant RHi [%] (double)"].as<string>(), "Constant RHi [%] (double)");
        input.MET_HUMIDSCAL_SCALING_A = parseDoubleString(humidScalingOpts["Humidity scaling constant a (double)"].as<string>(), "Humidity scaling constant a (double)");
        input.MET_HUMIDSCAL_SCALING_B = parseDoubleString(humidScalingOpts["Humidity scaling constant b (double)"].as<string>(), "Humidity scaling constant b (double)");


        YAML::Node moistLayerSubmenu = metNode["IMPOSE MOIST LAYER DEPTH SUBMENU"];
        input.MET_FIXDEPTH = parseBoolString(moistLayerSubmenu["Impose moist layer depth (T/F)"].as<string>(), "Impose moist layer depth (T/F)");
        input.MET_DEPTH = parseDoubleString(moistLayerSubmenu["Moist layer depth [m] (double)"].as<string>(), "Moist layer depth [m] (double)");

        YAML::Node lapseRateSubmenu = metNode["IMPOSE LAPSE RATE SUBMENU"];
        input.MET_FIXLAPSERATE = parseBoolString(lapseRateSubmenu["Impose lapse rate (T/F)"].as<string>(), "Impose lapse rate (T/F)");
        input.MET_LAPSERATE = parseDoubleString(lapseRateSubmenu["Lapse rate [K/m] (T/F)"].as<string>(), "Lapse rate [K/m] (T/F)");

        input.MET_DIURNAL = parseBoolString(metNode["Add diurnal variations (T/F)"].as<string>(), "Add diurnal variations (T/F)");
        
        YAML::Node tempPerturbMenu = metNode["TEMPERATURE PERTURBATION SUBMENU"];
        input.MET_ENABLE_TEMP_PERTURB = parseBoolString( tempPerturbMenu["Enable Temp. Pert. (T/F)"].as<string>(), "Enable Temp. Pert. (T/F)" );
        input.MET_TEMP_PERTURB_AMPLITUDE = parseDoubleString( tempPerturbMenu["Temp. Perturb. Amplitude (double)"].as<string>(), "Temp. Perturb. Amplitude (double)" );
        input.MET_TEMP_PERTURB_TIMESCALE = parseDoubleString( tempPerturbMenu["Temp. Perturb. Timescale (min)"].as<string>(), "Temp. Perturb. Timescale (min)" );
        
        //Humidity mod scheme must be none, constant, or scaling
        string modSchemeCaps = input.MET_HUMIDSCAL_MODIFICATION_SCHEME;
        for (auto & c: modSchemeCaps) c = toupper(c);

        if(modSchemeCaps != "NONE" && modSchemeCaps != "CONSTANT" && modSchemeCaps != "SCALING") {
            throw std::invalid_argument("Humidity modification scheme must be one of none, constant, or scaling.");
        }
        else {
            for (auto & c: input.MET_HUMIDSCAL_MODIFICATION_SCHEME) c = tolower(c);
        }

        //At least one of load met, impose depth, and fix lapse rate must be on.
        if(input.MET_LOADMET + input.MET_FIXDEPTH + input.MET_FIXLAPSERATE == 0){
            throw std::invalid_argument("At least one of \"Use met. input\", \"Impose moist layer depth\", or \"Impose lapse rate\" must be on!");
        }
        if(input.MET_FIXDEPTH + input.MET_FIXLAPSERATE > 1) {
            throw std::invalid_argument("Cannot fix both moist layer depth and lapse rate");
        }
        if(input.MET_LOADMET){
            if(input.MET_LOADTEMP + input.MET_FIXLAPSERATE + input.MET_FIXDEPTH != 1){
                throw std::invalid_argument("When using met. input, only one of \"Init temp. from met.\", \"Impose moist layer depth\", or \"Impose lapse rate\" must be selected.");
            }
            else if (input.MET_LOADRH + input.MET_FIXLAPSERATE + input.MET_FIXDEPTH != 1){
                throw std::invalid_argument("When using met. input, only one of \"Init RH from met.\" or \"Impose moist layer depth\", or \"Impose lapse rate\" must be selected.");
            }

        }
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
        input.ADV_GRID_NX = parseIntString(gridSubmenu["NX (positive int)"].as<string>(), "NX (positive int)");
        input.ADV_GRID_NY = parseIntString(gridSubmenu["NY (positive int)"].as<string>(), "NY (positive int)");
        input.ADV_GRID_XLIM_RIGHT = parseDoubleString(gridSubmenu["XLIM_RIGHT (positive double)"].as<string>(), "XLIM_RIGHT (positive double)");
        input.ADV_GRID_XLIM_LEFT = parseDoubleString(gridSubmenu["XLIM_LEFT (positive double)"].as<string>(), "XLIM_LEFT (positive double)");
        input.ADV_GRID_YLIM_UP = parseDoubleString(gridSubmenu["YLIM_UP (positive double)"].as<string>(), "YLIM_UP (positive double)");
        input.ADV_GRID_YLIM_DOWN = parseDoubleString(gridSubmenu["YLIM_DOWN (positive double)"].as<string>(), "YLIM_DOWN (positive double)");

        if(input.ADV_GRID_NX < 0 ||
           input.ADV_GRID_NY < 0 ||
           input.ADV_GRID_XLIM_LEFT < 0 || 
           input.ADV_GRID_XLIM_RIGHT < 0 ||
           input.ADV_GRID_YLIM_UP < 0 ||
           input.ADV_GRID_YLIM_DOWN < 0) {
            
            throw std::invalid_argument("No values in GRID SUBMENU can be less than zero!");
        }
    }

    vector<std::unordered_map<string, double>> generateCasesHelper(vector<std::unordered_map<string, double>>& allCases, const vector<std::pair<string, Vector_1D>>& params, const int row){
        if(row ==  params.size()){
            return allCases;
        }

        //For first row, just add maps with one entry for each case, these will be duplicated and filled up.
        string paramName = params[row].first;
        Vector_1D paramCases = params[row].second;
        if(allCases.empty()){
            for(const auto& c: paramCases){
                allCases.push_back(std::unordered_map<string, double>({{paramName, c}}));
            }   
        }
        else{
            //Generate new cases from existing ones as we go down the parameter list
            vector<std::unordered_map<string,double>> newCases;
            for(const auto& i: allCases){
                for(const auto& j: paramCases){
                    std::unordered_map<string,double> tempMap = i;
                    tempMap[paramName] = j;
                    newCases.push_back(tempMap);
                }
            }
            allCases = std::move(newCases);
        }
        //Repeat pattern of adding to the list of cases as we go down the param list until end and all cases have been generated
        return generateCasesHelper(allCases, params, row + 1);
    }

    vector<std::unordered_map<string, double>> generateCases(const OptInput& input){

        //Convert parameter map to vector, each row of the vector represents one parameter
        vector<std::pair<string, Vector_1D>> params;
        for (const auto& p: input.PARAMETER_PARAM_MAP){
            params.push_back(p);
        }
        vector<std::unordered_map<string, double>> allCases;
        return generateCasesHelper(allCases, params, 0);
    }

    Vector_1D parseParamSweepInput(const string paramString, const string paramLocation, bool monteCarlo, int nRuns){
        const string s = trim(paramString);
        const vector<string> colon_split_tokens = split(s, ":");
        if(monteCarlo){
            if(colon_split_tokens.size() != 0 && colon_split_tokens.size() != 2){
                throw std::invalid_argument("Monte Carlo Simulation requires parameter input format of min:max or a singular (constant) value at " + paramLocation + "!");
            }
            setSeed();
            Vector_1D paramVector;
            const double min = parseDoubleString(colon_split_tokens[0], "");
            const double max = parseDoubleString(colon_split_tokens[1], "");
            for(int i = 0; i < nRuns; i++){
                paramVector.push_back(fRand(min, max));
            }
            return paramVector;
        }
        // Non-monte carlo case
        if(colon_split_tokens.size() != 1 && colon_split_tokens.size() != 3){
            throw std::invalid_argument("Parameter sweep step input requires the format min:step:max at " + paramLocation + "!");
        }
        //min:step:max case
        if(colon_split_tokens.size() == 3){
            const double epsilon = 1e-40;
            const double min = parseDoubleString(colon_split_tokens[0], paramLocation);
            const double step = parseDoubleString(colon_split_tokens[1], paramLocation);
            const double max = parseDoubleString(colon_split_tokens[2], paramLocation);
            Vector_1D paramVector;
            for(double d = min; d < max; d+= step){
                paramVector.push_back(d);
            }
            if(abs(paramVector[paramVector.size()-1] - max ) > epsilon) paramVector.push_back(max);
            return paramVector;
        }
        //a b c d ... case
        const vector<string> space_split_tokens = split(s, " ");
        Vector_1D paramVector;
        for (string stri: space_split_tokens){
            paramVector.push_back(parseDoubleString(stri, paramLocation));
        }
        return paramVector;
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
        Vector_1D vec1d = parseParamSweepInput(paramString, paramLocation);
        vector<int> vecint;
        for (double d: vec1d){
            if(std::fmod(d,1.0) > 1e-40) throw (std::invalid_argument("Decimals not allowed in int inputs at " + paramLocation + "!"));
            vecint.push_back((int)(d));
        }
        return vecint;
    }

    std::string parseFileSystemPath(std::string str){
        std::filesystem::path p(str);
        return p.is_absolute() ? str : std::filesystem::weakly_canonical(INPUT_FILE_PATH.parent_path() / str).generic_string();
    }
}
