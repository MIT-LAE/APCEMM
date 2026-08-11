#include "Util/YamlUtils.hpp"

#include <stdexcept>
#include <string>

/*
Return the key as a string. Throw if the key is not a scalar.
YAML allows non-scalar keys. APCEMM inputs never use them, and the merge
below cannot handle them: the merge matches keys by string, and a
non-scalar key does not have a string form: Node::Scalar() returns "".
The lookup with node[key] which the merge relies on then searches for ""
and does not finds the real entry.
The merge concludes the key is missing, inserts it again, and duplicates it.
validateYamlKeys() already rejects non-scalar keys because there are none in the
defaults/input.yaml
This is for completeness so that the merge code does not depend on the default.
*/
const std::string& getScalarKey(const YAML::Node& key)
{
    if (key.IsScalar()) return key.Scalar();
    throw std::runtime_error("Found invalid YAML map key: map keys must be scalars.");
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
