#include "YamlInputReader/YamlPathUtils.hpp"
#include "YamlInputReader/YamlInputReader.hpp"  // split
#include "Util/YamlUtils.hpp"                   // getScalarKey, buildKeyPath

#include <cstddef>
#include <filesystem>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

/*
This module exists to resolve relative paths in APCEMM input files.
APCEMM resolves a relative path relative to the directory of the input file
that specifies the path. This is fine for a single input file, but when they
are multiple we need to ensure that each path is resolved against the dir of
the input file it came from. During the parsing of a user input file, all 
relative paths are resolved to absolute paths before they get merged.

To do so we:
1. Keep track of which keys in the input yaml are paths
2. Check that set of path keys against the default yaml at runtime
3. Check that there are no relative paths left after path resolution 
*/
namespace YamlInputReader{
    // The input.yaml keys whose value is a file system path. 
    // The entries are a 'keyPath', built with the same " -> " separator buildKeyPath uses.
    // This set needs to be kept up to date with defaults/input.yaml, every path that
    // appears in the default path must also appear here to be checked.
    // To ensure this set is up to date:
    // checkDefaultPaths() checks that all keys in the PATH_KEYS set
    // are also in default/input.yaml so a key rename cannot bypass the check.
    // readPath() throws on a value left relative, so a path key missing
    // from this set returns an error instead of potentially resolving against the wrong
    // directory.
    static const std::set<std::string> PATH_KEYS = {
        "SIMULATION MENU -> OUTPUT SUBMENU -> Output folder (string)",
        "SIMULATION MENU -> Input background condition (string)",
        "SIMULATION MENU -> Input engine emissions (string)",
        "SIMULATION MENU -> External EPM NetCDF file",
        "METEOROLOGY MENU -> METEOROLOGICAL INPUT SUBMENU -> Met input file path (string)",
    };

    // Sentinel values used in the default/input.yaml
    bool isPathSentinel(const std::string& value) {
        return value == "=MISSING=" || value == "=DEFAULT=";
    }

    // Traverse all nodes of a YAML file, and resolve all keys contained
    // in path keys to absolute paths
    void resolvePaths(YAML::Node node, const std::filesystem::path& baseDir,
                        const std::string& currentPath) {
        if (!node.IsMap()) {
            return;
        }
        for (const auto& it : node) {
            const std::string key = getScalarKey(it.first);
            const std::string keyPath = buildKeyPath(currentPath, key);
            YAML::Node value = it.second;

            if (value.IsMap()) {
                resolvePaths(value, baseDir, keyPath);
                continue;
            }
            // Anything that is not a scalar is either a null value, which means "no
            // override", or a type error that validateYamlKeys has already rejected.
            if (!value.IsScalar() || !PATH_KEYS.contains(keyPath)) {
                continue;
            }

            const std::string pathAsWritten = value.as<std::string>();
            if (isPathSentinel(pathAsWritten)) {
                continue;
            }
            const std::filesystem::path written(pathAsWritten);
            if (written.is_absolute()) {
                continue;
            }
            // baseDir is computed using std::filesystem::parent_path(input_file_path).
            // It is empty when the input file was named without a directory (e.g. passed
            // to APCEMM as './path/to/APCEMM input.yaml'). This means that when baseDir is empty,
            // the input file is in the same dir as the current working directory (assuming no
            // part of APCEMM changes the CWD which is true).
            // weakly_canonical() only makes absolute the leading portion of the path that
            // already exists. A value whose first element does not exist yet, such as an
            // output folder APCEMM creates later (APCEMM_out/ in input.yaml in CWD), stays
            // relative and readPath() then rejects it.
            // Set the base explicitly so the result is absolute either way.
            const std::filesystem::path base = baseDir.empty() ? std::filesystem::current_path() : baseDir;
            // Update the node value to be the resolved absolute path
            // Direct assignment updates the node inplace
            std::filesystem::path resolved = std::filesystem::weakly_canonical(base / written);
            // Normalize trailing / at the end of a directory path
            if (!resolved.has_filename()) {
                resolved = resolved.parent_path();
            }
            value = resolved.generic_string();
        }
        }

        // Get a value from the YAML node tree from a key path but as a vector of keys e.g.
        // ['METEOROLOGY MENU', 'METEOROLOGICAL INPUT SUBMENU' 'Met input file path (string)']
        std::string getValueFromKeyPath(const YAML::Node& node, const std::vector<std::string>& keys, std::size_t depth) {
            if (depth == keys.size()) {
                return node.IsScalar() ? node.as<std::string>() : "";
            }
            if (!node.IsMap()) {
                return "";
            }
            const YAML::Node child = node[keys[depth]];
            if (!child) {
                return "";
            }
            return getValueFromKeyPath(child, keys, depth + 1);
        }

        // Get a value from the YAML node tree from a key path e.g. (splits onto ' -> ')
        // 'METEOROLOGY MENU -> METEOROLOGICAL INPUT SUBMENU -> Met input file path (string)'
        std::string getValueFromKeyPath(const YAML::Node& node, const std::string& keyPath) {
            return getValueFromKeyPath(node, split(keyPath, " -> "), 0);
        }

    // Runs on one input file at a time, before the merge
    void resolvePathsInPlace(YAML::Node node, const std::filesystem::path& baseDir) {
        resolvePaths(node, baseDir, "");
    }

    // The defaults are compiled in, so this only fails when the repository is in a bad
    // state. It checks that PATH_KEYS is current with default/input.yaml
    void checkDefaultPaths(const YAML::Node& defaultData) {
        for (const std::string& keyPath : PATH_KEYS) {
            const std::string value = getValueFromKeyPath(defaultData, keyPath);
            // Here the default/input.yaml has been updated but the set PATH_KEYS has not
            if (value.empty()) {
                throw std::runtime_error("PATH_KEYS names '" + keyPath + "', which is not a scalar leaf of the "
                                         "default input.yaml. Path resolution would silently skip it, update PATH_KEYS.");
            }
            // Here the default/input.yaml has been updated to contain unsupported values
            if (!isPathSentinel(value) && !std::filesystem::path(value).is_absolute()) {
                throw std::runtime_error("The default input.yaml sets '" + keyPath + "' to the relative path '" +
                                         value + "'. The compiled-in defaults have no directory to resolve against, "
                                         "so they must hold absolute paths or the '=MISSING='/'=DEFAULT=' sentinels.");
            }
        }
    }

    // Final check before passing a key that is supposed to be a path to OptInput,
    // if the path is somehow still relative, it means it was not added to PATH_KEYS
    // to be resolved -> throw because otherwise we might resolve the relative path incorrectly
    std::string readPath(const YAML::Node& node, const std::string& key){
        const std::string value = node[key].as<std::string>();
        if (isPathSentinel(value)|| std::filesystem::path(value).is_absolute()) {
            return value;
        }
        throw std::runtime_error("Path '" + key + "' is still relative after reading the input: '" + value +
                                 "'. Add its full key path to PATH_KEYS in YamlPathUtils.cpp so it is "
                                 "resolved against the file that declares it.");
    }
}
