#include "Util/YamlUtils.hpp"

#include <stdexcept>
#include <string>

/*
See YamlUtils.hpp.
*/
const std::string& getScalarKey(const YAML::Node& key)
{
    if (key.IsScalar()) return key.Scalar();
    // Keep the position that yaml-cpp records, so the user can find the key.
    // Nodes built in memory instead of parsed from a file have no position.
    const YAML::Mark mark = key.Mark();
    const std::string where = mark.is_null()
        ? ""
        : " at line " + std::to_string(mark.line + 1) + ", column " + std::to_string(mark.column + 1);
    throw std::runtime_error("Found invalid YAML map key" + where + ": map keys must be scalars.");
}

// Adapted from https://stackoverflow.com/a/66205210
const YAML::Node mergeYamlNodes(const YAML::Node& defaults, const YAML::Node& overrides)
{
    // If overrides is not a map, merge result is overrides, unless overrides
    // is null.
    if (!overrides.IsMap()) return overrides.IsNull() ? defaults : overrides;

    // If defaults is not a map, merge result is overrides.
    // This can replace the node type of default by the one in override
    // (e.g. default not a map -> override a map) which does not make sense for an input
    // file: there is one valid structure.
    // In practice this does not happen because we have validated the YAML inputs upstream
    // so that both files have the same structure, and this function is more generic.
    if (!defaults.IsMap()) return overrides;

    // Create a new map 'newNode' with the same mappings as defaults, merged
    // with overrides.
    auto newNode = YAML::Node(YAML::NodeType::Map);
    for (auto node : defaults) {
        const std::string& key = getScalarKey(node.first);
        if (overrides[key]) {
            newNode[node.first] = mergeYamlNodes(node.second, overrides[key]);
        }
        else {
            newNode[node.first] = node.second;
        }
    }

    // Add the mappings from 'overrides' not already in 'newNode'. Keys present
    // in both maps were merged by the loop above, so they must be skipped here:
    // Node::operator[](const Node&) matches by node identity, not by string, so
    // re-assigning them would append a second entry with an equal key instead of
    // overwriting.
    // The lookup goes through a const reference on purpose: the non-const
    // Node::operator[] inserts an undefined placeholder pair for a missing key,
    // the const one does not.
    const YAML::Node& mergedSoFar = newNode;
    for (auto node : overrides) {
        const std::string& key = getScalarKey(node.first);
        if (!mergedSoFar[key]) {
            newNode[node.first] = node.second;
        }
    }

    return YAML::Node(newNode);
}
