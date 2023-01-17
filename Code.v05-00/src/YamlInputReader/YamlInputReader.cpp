#include "YamlInputReader/YamlInputReader.hpp"
using std::cout;
using std::endl;
namespace YamlInputReader{
    void readYamlInputFile(OptInput& input, string filename){
        YAML::Node data = YAML::LoadFile(filename);
        readSimMenu(input, data["SIMULATION MENU"]);
        readParamMenu(input, data["PARAMETER MENU"]);
        readTransportMenu(input, data["TRANSPORT MENU"]);
        readChemMenu(input, data["CHEMISTRY MENU"]);
        readAeroMenu(input, data["AEROSOL MENU"]);
        readMetMenu(input, data["METEOROLOGY MENU"]);
        readDiagMenu(input, data["DIAGNOSTIC MENU"]);
    }
    void readSimMenu(OptInput& input, const YAML::Node& simNode){
        YAML::Node paramSweepSubmenu = simNode["PARAM SWEEP SUBMENU"];
        input.SIMULATION_PARAMETER_SWEEP = parseBoolString(paramSweepSubmenu["Parameter sweep (T/F)"].as<string>(), "Parameter sweep (T/F");
        input.SIMULATION_MONTECARLO = parseBoolString(paramSweepSubmenu["Run Monte Carlo (T/F)"].as<string>(), "Run Monte Carlo (T/F)");
        input.SIMULATION_MCRUNS =  parseIntString(paramSweepSubmenu["Num Monte Carlo runs (int)"].as<string>(), "Num Monte Carlo runs (int)");

        YAML::Node outputSubmenu = simNode["OUTPUT SUBMENU"];
        input.SIMULATION_OUTPUT_FOLDER = outputSubmenu["Output folder (string)"].as<string>();
        input.SIMULATION_OVERWRITE = parseBoolString(outputSubmenu["Overwrite if folder exists (T/F)"].as<string>(), "Overwrite if folder exists (T/F)");
        
        input.SIMULATION_RUN_DIRECTORY = simNode["Run directory (string)"].as<string>();
        input.SIMULATION_THREADED_FFT = parseBoolString(simNode["Use threaded FFT (T/F)"].as<string>(), "Use threaded FFT (T/F)");

        YAML::Node fftwWisdomSubmenu = simNode["FFTW WISDOM SUBMENU"];
        input.SIMULATION_USE_FFTW_WISDOM = parseBoolString(fftwWisdomSubmenu["Use FFTW WISDOM (T/F)"].as<string>(), "Use FFTW WISDOM (T/F)");
        input.SIMULATION_DIRECTORY_W_WRITE_PERMISSION = fftwWisdomSubmenu["Dir w/ write permission (string)"].as<string>();
        
        input.SIMULATION_INPUT_BACKG_COND = simNode["Input background condition (string)"].as<string>();

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
        input.PARAMETER_FILEINPUT = parseBoolString(paramNode["Use input file (T/F)"].as<string>(), "Use input file (T/F)");
        input.PARAMETER_FILENAME = paramNode["Input file name (string)"].as<string>();

        input.PARAMETER_PARAM_MAP["PLUMEPROCESS"] = parseParamSweepInput(paramNode["Plume Process [hr] (double)"].as<string>(), "Plume Process [hr] (double)");

        YAML::Node metParamSubmenu = paramNode["METEOROLOGICAL PARAMETERS SUBMENU"];
        input.PARAMETER_PARAM_MAP["TEMPERATURE"] = parseParamSweepInput(metParamSubmenu["Temperature [K] (double)"].as<string>(), "Temperature [K] (double)");
        input.PARAMETER_PARAM_MAP["RHW"] = parseParamSweepInput(metParamSubmenu["R.Hum. wrt water [%] (double)"].as<string>(), "R.Hum. wrt water [%] (double)");
        input.PARAMETER_PARAM_MAP["PRESSURE"] = parseParamSweepInput(metParamSubmenu["Pressure [hPa] (double)"].as<string>(), "Pressure [hPa] (double)");
        input.PARAMETER_PARAM_MAP["DH"] = parseParamSweepInput(metParamSubmenu["Horiz. diff. coeff. [m^2/s] (double)"].as<string>(), "Horiz. diff. coeff. [m^2/s] (double)");
        input.PARAMETER_PARAM_MAP["DV"] = parseParamSweepInput(metParamSubmenu["Verti. diff. [m^2/s] (double)"].as<string>(), "Verti. diff. [m^2/s] (double)");
        input.PARAMETER_PARAM_MAP["SHEAR"] = parseParamSweepInput(metParamSubmenu["Wind shear [1/s] (double)"].as<string>(), "Wind shear [1/s] (double)");

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
        input.TRANSPORT_PART_FLUX = parseBoolString(transportNode["Particle flux correction (T/F)"].as<string>(), "Particle flux correction (T/F)");

        YAML::Node updraftSubmenu = transportNode["PLUME UPDRAFT SUBMENU"];
        input.TRANSPORT_UPDRAFT = parseBoolString(updraftSubmenu["Turn on plume updraft (T/F)"].as<string>(), "Turn on plume updraft (T/F)");
        input.TRANSPORT_UPDRAFT_TIMESCALE = parseDoubleString(updraftSubmenu["Updraft timescale [s] (double)"].as<string>(), "Updraft timescale [s] (double)");
        input.TRANSPORT_UPDRAFT_VELOCITY = parseDoubleString(updraftSubmenu["Updraft veloc. [cm/s] (double)"].as<string>(), "Updraft veloc. [cm/s] (double)");
    }
    void readChemMenu(OptInput& input, const YAML::Node& chemNode){
        input.CHEMISTRY_CHEMISTRY = parseBoolString(chemNode["Turn on Chemistry (T/F)"].as<string>(), "Turn on Chemistry (T/F)");
        input.CHEMISTRY_HETCHEM = parseBoolString(chemNode["Perform hetero. chem. (T/F)"].as<string>(), "Perform hetero. chem. (T/F)");
        input.CHEMISTRY_TIMESTEP = parseDoubleString(chemNode["Chemistry Timestep [min] (double)"].as<string>(), "Chemistry Timestep [min] (double)");
        input.CHEMISTRY_JRATE_FOLDER = chemNode["Photolysis rates folder (string)"].as<string>();
    }
    void readAeroMenu(OptInput& input, const YAML::Node& aeroNode){
        input.AEROSOL_GRAVSETTLING = parseBoolString(aeroNode["Turn on grav. settling (T/F)"].as<string>(), "Turn on grav. settling (T/F)");
        input.AEROSOL_COAGULATION_SOLID = parseBoolString(aeroNode["Turn on solid coagulation (T/F)"].as<string>(), "Turn on solid coagulation (T/F)");
        input.AEROSOL_COAGULATION_LIQUID = parseBoolString(aeroNode["Turn on liquid coagulation (T/F)"].as<string>(), "Turn on liquid coagulation (T/F)");
        input.AEROSOL_COAGULATION_TIMESTEP = parseDoubleString(aeroNode["Coag. timestep [min] (double)"].as<string>(), "Coag. timestep [min] (double)");
        input.AEROSOL_ICE_GROWTH = parseBoolString(aeroNode["Turn on ice growth (T/F)"].as<string>(), "Turn on ice growth (T/F)");
    }
    void readMetMenu(OptInput& input, const YAML::Node& metNode){
        YAML::Node metInputSubmenu = metNode["METEOROLOGICAL INPUT SUBMENU"];
        input.MET_LOADMET = parseBoolString(metInputSubmenu["Use met. input (T/F)"].as<string>(), "Use met. input (T/F)");
        input.MET_FILENAME = metInputSubmenu["Met input file path (string)"].as<string>();
        input.MET_DT = parseDoubleString(metInputSubmenu["Time series data timestep [hr] (double)"].as<string>(), "Time series data timestep [hr] (double)");
        input.MET_LOADTEMP = parseBoolString(metInputSubmenu["Init temp. from met. (T/F)"].as<string>(), "Init temp. from met. (T/F)");
        input.MET_TEMPTIMESERIES = parseBoolString(metInputSubmenu["Temp. time series input (T/F)"].as<string>(), "Temp. time series input (T/F)");
        input.MET_LOADRH = parseBoolString(metInputSubmenu["Init RH from met. (T/F)"].as<string>(), "Init RH from met. (T/F)");
        input.MET_RHTIMESERIES = parseBoolString(metInputSubmenu["RH time series input (T/F)"].as<string>(), "RH time series input (T/F)");
        input.MET_LOADSHEAR = parseBoolString(metInputSubmenu["Init wind shear from met. (T/F)"].as<string>(), "Init wind shear from met. (T/F)");
        input.MET_SHEARTIMESERIES = parseBoolString(metInputSubmenu["Wind shear time series input (T/F)"].as<string>(), "Wind shear time series input (T/F)");

        YAML::Node moistLayerSubmenu = metNode["IMPOSE MOIST LAYER DEPTH SUBMENU"];
        input.MET_FIXDEPTH = parseBoolString(moistLayerSubmenu["Impose moist layer depth (T/F)"].as<string>(), "Impose moist layer depth (T/F)");
        input.MET_DEPTH = parseDoubleString(moistLayerSubmenu["Moist layer depth [m] (double)"].as<string>(), "Moist layer depth [m] (double)");

        YAML::Node lapseRateSubmenu = metNode["IMPOSE LAPSE RATE SUBMENU"];
        input.MET_FIXLAPSERATE = parseBoolString(lapseRateSubmenu["Impose lapse rate (T/F)"].as<string>(), "Impose lapse rate (T/F)");
        input.MET_LAPSERATE = parseDoubleString(lapseRateSubmenu["Lapse rate [K/m] (T/F)"].as<string>(), "Lapse rate [K/m] (T/F)");

        input.MET_DIURNAL = parseBoolString(metNode["Add diurnal variations (T/F)"].as<string>(), "Add diurnal variations (T/F)");

        //At least one of load met, impose depth, and fix lapse rate must be on.
        if(input.MET_LOADMET + input.MET_FIXDEPTH + input.MET_FIXLAPSERATE == 0){
            throw std::invalid_argument("At least one of \"Use met. input\", \"Impose moist layer depth\", or \"Impose lapse rate\" must be on!");
        }
        if(input.MET_LOADMET){
            if(input.MET_LOADTEMP + input.MET_FIXLAPSERATE != 1){
                throw std::invalid_argument("When using met. input, only one of \"Init temp. from met.\" or \"Impose lapse rate\" must be selected.");
            }
            else if (input.MET_LOADRH + input.MET_FIXDEPTH != 1){
                throw std::invalid_argument("When using met. input, only one of \"Init RH from met.\" or \"Impose moist layer depth\" must be selected.");
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
    void performOtherInputValidnessChecks(OptInput& input){

    }
    vector<std::unordered_map<string, double>> generateCasesHelper(vector<std::unordered_map<string, double>>& allCases, const vector<std::pair<string, Vector_1D>>& params, const int row){
        if(row ==  params.size()){
            return allCases;
        }
        string paramName = params[row].first;
        Vector_1D paramCases = params[row].second;
        if(allCases.empty()){
            for(int i = 0; i < paramCases.size(); i++){
                allCases.push_back(std::unordered_map<string, double>({{paramName, paramCases[i]}}));
            }   
        }
        else{
            vector<std::unordered_map<string,double>> newCases;
            for(int i = 0; i < allCases.size(); i++){
                for(int j = 0; j < paramCases.size(); j++){
                    std::unordered_map<string,double> tempMap = allCases[i];
                    tempMap[paramName] = paramCases[j];
                    newCases.push_back(tempMap);
                }
            }
            allCases = newCases;
        }
        return generateCasesHelper(allCases, params, row+1);
    }

    vector<std::unordered_map<string, double>> generateCases(const OptInput& input){
        vector<std::pair<string, Vector_1D>> params;
        for (std::pair<string, Vector_1D> p: input.PARAMETER_PARAM_MAP){
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
        for (RealDouble d: vec1d){
            if(std::fmod(d,1.0) > 1e-40) throw (std::invalid_argument("Decimals not allowed in int inputs at " + paramLocation + "!"));
            vecint.push_back((int)(d));
        }
        return vecint;
    }
}
