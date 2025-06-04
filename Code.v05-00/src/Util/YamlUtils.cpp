#include "Util/YamlUtils.hpp"

#include <string>


// Adapted from https://stackoverflow.com/a/66205210

const YAML::Node mergeYamlNodes(const YAML::Node& defaults, const YAML::Node& overrides)
{
    // If overrides is not a map, merge result is overrides, unless overrides
    // is null.
    if (!overrides.IsMap()) return overrides.IsNull() ? defaults : overrides;

    // If defaults is not a map, merge result is overrides.
    if (!defaults.IsMap()) return overrides;

    if (!defaults.size()) return YAML::Node(overrides);

    // Create a new map 'newNode' with the same mappings as defaults, merged
    // with overrides.
    auto newNode = YAML::Node(YAML::NodeType::Map);
    for (auto node : defaults) {
        if (node.first.IsScalar()) {
            const std::string& key = node.first.Scalar();
            if (overrides[key]) {
                newNode[node.first] = mergeYamlNodes(node.second, overrides[key]);
                continue;
            }
        }
        newNode[node.first] = node.second;
    }

    // Add the mappings from 'overrides' not already in 'newNode'.
    for (auto node : overrides) {
        if (!node.first.IsScalar()) {
            const std::string& key = node.first.Scalar();
            if (defaults[key]) {
                newNode[node.first] = mergeYamlNodes(defaults[key], node.second);
                continue;
            }
        }
        newNode[node.first] = node.second;
    }

    return YAML::Node(newNode);
}
