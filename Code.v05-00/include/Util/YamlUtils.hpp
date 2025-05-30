#ifndef YAMLUTILS_H_
#define YAMLUTILS_H_

#include <yaml-cpp/yaml.h>

const YAML::Node mergeYamlNodes(const YAML::Node& defaults, const YAML::Node& overrides);

#endif // YAMLUTILS_H_
