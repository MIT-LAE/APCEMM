#ifndef YAMLINPUTREADER_H
#define YAMLINPUTREADER_H
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <yaml-cpp/yaml.h> 
#include <filesystem>
#include "Core/Input_Mod.hpp"
#include "Util/ForwardDecl.hpp"
#include "Util/MC_Rand.hpp"   
using std::string;
using std::vector;

namespace YamlInputReader{
    static std::filesystem::path INPUT_FILE_PATH;
    void readYamlInputFile(OptInput& input, string filename);
    void readSimMenu(OptInput& input, const YAML::Node& inputNode);
    void readParamMenu(OptInput& input, const YAML::Node& paramNode);
    void readTransportMenu(OptInput& input, const YAML::Node& transportNode);
    void readChemMenu(OptInput& input, const YAML::Node& chemNode);
    void readAeroMenu(OptInput& input, const YAML::Node& aeroNode);
    void readMetMenu(OptInput& input, const YAML::Node& metNode);
    void readDiagMenu(OptInput& input, const YAML::Node& diagNode);
    void readAdvancedMenu(OptInput& input, const YAML::Node& advancedNode);
    
    void performOtherInputValidnessChecks(OptInput& input);
    vector<std::unordered_map<string, double>> generateCases(const OptInput& input);
    Vector_1D parseParamSweepInput(const string paramString, const string paramLocation = "", bool monteCarlo = false, int nRuns = 0);
    vector<string> split(const string str, const string delimiter);

    inline string trim(const string str){
        string s = str;
        //trim left
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), \
                                        [](unsigned char ch) {
                                            return !std::isspace(ch);
                                        }));

        //trim right
        s.erase(std::find_if(s.rbegin(), s.rend(), 
                            [](unsigned char ch) {
                                return !std::isspace(ch);
                            }).base(), s.end());
        return s;
    }
    inline double strToDouble(const string str, bool sciNotation = true){
        string s = trim(str);
        bool foundDot = false;
        if (s.length() != 0 && (s[0] == '-' || s[0] == '+')){
            s = s.substr(1);
        }
        if (s.length() == 0) throw (std::invalid_argument("Input Reader tried to process number of length 0! Please double-check your input!"));
        for (std::size_t i = 0; i < s.length(); i++){
            char c = s[i];
            if(std::isdigit(c)) continue;
            //Don't allow something like 10e2. or 3e-2.5
            else if (sciNotation && (c == '.' && !foundDot)){
                foundDot = true;
                continue;
            } 
            else if  (sciNotation && (c == 'e' || c == 'E')){
                return 0*strToDouble(s.substr(i+1), false) + std::stod(str);
            }
            else{
                throw (std::invalid_argument("Something went wrong with processing a number. Please double-check your parameter input format"));
            }
        }
        return std::stod(str);
    }
    inline bool parseBoolString(const string paramString, const string paramLocation = ""){
        string str = trim(paramString);
        std::transform(str.begin(), str.end(), str.begin(), 
                    [](unsigned char c){
                            return std::toupper(c);
                        });
        
        if(str == "T" || str == "TRUE" || str == "YES" || str == "Y" || str == "1"){ return true; }
        else if (str == "F" || str == "FALSE" || str == "NO" || str == "N" || str == "0"){ return false; }
        else{
            throw std::invalid_argument("Unable to read boolean value at parameter " + paramLocation); 
        }
    }
    inline double parseDoubleString(const string paramString, const string paramLocation = ""){
        try{
            return strToDouble(paramString);
        }
        catch (const std::invalid_argument &e){
            throw std::invalid_argument(string(e.what()) + " at parameter " + paramLocation);
        }

    }
    inline int parseIntString(const string paramString, const string paramLocation = ""){
        return (int)(parseDoubleString(paramString, paramLocation));
    }

    vector<int> parseVectorIntString(const string paramString, const string paramLocation = "");
    std::string parseFileSystemPath(std::string str);
}
#endif