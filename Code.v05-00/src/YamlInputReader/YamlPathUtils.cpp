#include "YamlInputReader/YamlPathUtils.hpp"
#include "YamlInputReader/YamlInputReader.hpp"  // split, parseBoolString, iequals
#include "Util/YamlUtils.hpp"                   // getScalarKey

#include <cstddef>
#include <filesystem>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace YamlInputReader{
    // The keys whose value is a file system path. A relative value under one of these
    // keys resolves against the directory of the file that declares it.
    // The entries are full key paths, built with the same " -> " separator
    // validateYamlKeys uses, because a leaf key alone is not unique: "netCDF filename
    // format (string)" appears four times, and neither it nor "Inst timeseries file
    // (string)" is a path, they are filename patterns.
    // Keep this in step with defaults/input.yaml. checkDefaultPaths() rejects an entry
    // that names no leaf of the default tree, so a key rename cannot silently break
    // resolution, and readPath() throws on a value left relative, so a path key missing
    // from this set fails loudly on the first run instead of resolving against the wrong
    // directory.
    static const std::set<std::string> PATH_KEYS = {
        "SIMULATION MENU -> OUTPUT SUBMENU -> Output folder (string)",
        "SIMULATION MENU -> FFTW WISDOM SUBMENU -> Dir w/ write permission (string)",  // removed by TASK-56
        "SIMULATION MENU -> Input background condition (string)",
        "SIMULATION MENU -> Input engine emissions (string)",
        "SIMULATION MENU -> External EPM NetCDF file",
        "CHEMISTRY MENU -> Photolysis rates folder (string)",                          // removed by TASK-56
        "METEOROLOGY MENU -> METEOROLOGICAL INPUT SUBMENU -> Met input file path (string)",
    };

    bool isPathSentinel(const std::string& value) {
        return value == "=MISSING=" || value == "=DEFAULT=";
    }

    // Internal to path handling. Everything the rest of the reader calls is declared in
    // YamlPathUtils.hpp.
    namespace {
        // The resolution rule, in one place: a relative path resolves against baseDir, the
        // directory of the file that declares it. Nothing else.
        // currentPath is the full key path of the node being walked, which is what
        // PATH_KEYS holds, so the recursion carries it and the public entry point does not
        // have to expose it.
        void resolvePaths(YAML::Node node, const std::filesystem::path& baseDir,
                          const std::string& sourceFile, PathOriginMap& origins,
                          const std::string& currentPath) {
            if (!node.IsMap()) {
                return;
            }
            for (const auto& it : node) {
                const std::string key = getScalarKey(it.first);
                const std::string keyPath = currentPath.empty() ? key : currentPath + " -> " + key;
                YAML::Node value = it.second;

                if (value.IsMap()) {
                    resolvePaths(value, baseDir, sourceFile, origins, keyPath);
                    continue;
                }
                // Anything that is not a scalar is either a null value, which means "no
                // override", or a type error that validateYamlKeys has already rejected.
                if (!value.IsScalar() || !PATH_KEYS.contains(keyPath)) {
                    continue;
                }

                const std::string asWritten = value.as<std::string>();
                if (isPathSentinel(asWritten)) {
                    continue;
                }
                // Recorded whether or not the path is rewritten, so an error about an
                // absolute path can still name the file that declared it.
                origins[keyPath] = PathOrigin{asWritten, sourceFile};

                const std::filesystem::path written(asWritten);
                if (written.is_absolute()) {
                    continue;
                }
                // parent_path() is empty when the input file is named without a directory,
                // and weakly_canonical() would then leave the result relative.
                const std::filesystem::path base = baseDir.empty() ? std::filesystem::current_path() : baseDir;
                node[key] = std::filesystem::weakly_canonical(base / written).generic_string();
            }
        }

        // Name the path the way the user can act on it: the key, the text they wrote,
        // where they wrote it, and what it resolved to. A path that no input file
        // declares comes from the compiled-in defaults.
        std::string describePath(const std::string& keyPath, const std::string& resolved,
                                 const PathOriginMap& origins) {
            const auto it = origins.find(keyPath);
            const std::string declaredIn = it == origins.end()
                ? "the compiled-in default input.yaml"
                : "'" + it->second.declaredIn + "'";
            const std::string asWritten = it == origins.end() ? resolved : it->second.asWritten;
            std::string description = "'" + keyPath + "' is set to '" + asWritten + "' in " + declaredIn;
            if (asWritten != resolved) {
                description += ", which resolves to '" + resolved + "'";
            }
            return description;
        }

        // Read a leaf by its full key path. Returns an empty string when the path names
        // nothing or names a node that is not a scalar.
        // The walk is recursive and the node is taken by const reference on purpose:
        // "node = node[key]" would not rebind the local node, it would assign the child
        // into the tree the node came from, and the non-const Node::operator[] inserts an
        // undefined placeholder pair for a missing key while the const one does not.
        std::string leafAt(const YAML::Node& node, const std::vector<std::string>& keys, std::size_t depth) {
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
            return leafAt(child, keys, depth + 1);
        }

        std::string leafAt(const YAML::Node& node, const std::string& keyPath) {
            return leafAt(node, split(keyPath, " -> "), 0);
        }

        void requirePathExists(const YAML::Node& merged, const std::string& keyPath,
                               const PathOriginMap& origins) {
            const std::string resolved = leafAt(merged, keyPath);
            if (isPathSentinel(resolved) || std::filesystem::exists(resolved)) {
                return;
            }
            throw std::runtime_error("Input file not found: " + describePath(keyPath, resolved, origins) +
                                     ". A relative path is resolved against the directory of the file that "
                                     "sets it, not against the working directory.");
        }
    }

    // Runs on one input file at a time, before the merge, because the merge is what
    // loses track of which file a value came from.
    void resolvePathsInPlace(YAML::Node node, const std::filesystem::path& baseDir,
                             const std::string& sourceFile, PathOriginMap& origins) {
        resolvePaths(node, baseDir, sourceFile, origins, "");
    }

    // One pass over the merged tree, so a path is reported while reading the input
    // instead of much later, when a component fails to open it.
    // A path is only checked when the run actually reads it. The output folder is never
    // checked: Main.cpp creates it.
    void checkPathsExist(const YAML::Node& merged, const PathOriginMap& origins) {
        if (parseBoolString(leafAt(merged, "METEOROLOGY MENU -> METEOROLOGICAL INPUT SUBMENU -> Use met. input (T/F)"),
                            "Use met. input (T/F)")) {
            requirePathExists(merged, "METEOROLOGY MENU -> METEOROLOGICAL INPUT SUBMENU -> Met input file path (string)", origins);
        }
        if (iequals(leafAt(merged, "SIMULATION MENU -> EPM type (original/external/new)"), "external")) {
            requirePathExists(merged, "SIMULATION MENU -> External EPM NetCDF file", origins);
        }
        requirePathExists(merged, "SIMULATION MENU -> Input background condition (string)", origins);
        requirePathExists(merged, "SIMULATION MENU -> Input engine emissions (string)", origins);
    }

    // The defaults are compiled in, so this only fails when the repository is in a bad
    // state. It is what keeps PATH_KEYS honest: an entry that names no leaf would stop
    // resolving without any other sign, and a relative default would need a second
    // resolution rule because the compiled-in defaults have no directory of their own.
    void checkDefaultPaths(const YAML::Node& defaultData) {
        for (const std::string& keyPath : PATH_KEYS) {
            const std::string value = leafAt(defaultData, keyPath);
            if (value.empty()) {
                throw std::runtime_error("PATH_KEYS names '" + keyPath + "', which is not a scalar leaf of the "
                                         "default input.yaml. Path resolution would silently skip it.");
            }
            if (!isPathSentinel(value) && !std::filesystem::path(value).is_absolute()) {
                throw std::runtime_error("The default input.yaml sets '" + keyPath + "' to the relative path '" +
                                         value + "'. The compiled-in defaults have no directory to resolve against, "
                                         "so they must hold absolute paths or the '=MISSING=' sentinel.");
            }
        }
    }

    // Every path reaches OptInput through this accessor. A value left relative means the
    // key is not in PATH_KEYS, so it was never resolved and would be read against the
    // working directory. Refusing it here turns that mistake into a clear error on the
    // first run instead of a file that is silently looked for in the wrong place.
    std::string readPath(const YAML::Node& node, const std::string& key){
        const std::string value = node[key].as<std::string>();
        if (isPathSentinel(value) || std::filesystem::path(value).is_absolute()) {
            return value;
        }
        throw std::runtime_error("Path '" + key + "' is still relative after reading the input: '" + value +
                                 "'. Add its full key path to PATH_KEYS in YamlPathUtils.cpp so it is "
                                 "resolved against the file that declares it.");
    }
}
