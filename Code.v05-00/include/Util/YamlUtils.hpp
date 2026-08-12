#ifndef YAMLUTILS_H_
#define YAMLUTILS_H_

#include <string>
#include <yaml-cpp/yaml.h>

/*
Return the key as a string. Throw if the key is not a scalar.
YAML allows non-scalar keys (a sequence or a map used as a key). APCEMM inputs
never use them, and every part of the reader matches keys by their string form,
which a non-scalar key does not have: Node::Scalar() returns "".
Every place that turns a map key into a string must go through this function.
Node::as<std::string>() throws a yaml-cpp "bad conversion" instead, which does
not tell the user what is wrong with their input file.
*/
const std::string& getScalarKey(const YAML::Node& key);

const YAML::Node mergeYamlNodes(const YAML::Node& defaults, const YAML::Node& overrides);

#endif // YAMLUTILS_H_
