#ifndef YAMLPATHUTILS_H_
#define YAMLPATHUTILS_H_

#include <filesystem>
#include <string>
#include <yaml-cpp/yaml.h>

/*
See description in YamlPathUtils.cpp
*/
namespace YamlInputReader {
    bool isPathSentinel(const std::string& value);

    void resolvePathsInPlace(YAML::Node node, const std::filesystem::path& baseDir);

    void checkDefaultPaths(const YAML::Node& defaultData);

    std::string readPath(const YAML::Node& node, const std::string& key);
}

#endif // YAMLPATHUTILS_H_
