#ifndef YAMLPATHUTILS_H_
#define YAMLPATHUTILS_H_

#include <filesystem>
#include <map>
#include <string>
#include <yaml-cpp/yaml.h>

/*
Everything APCEMM knows about file system paths in an input file.

ONE rule: a relative path resolves against the directory of the file that
declares it. Not the working directory, and not the last file passed on the
command line. The compiled-in defaults hold no relative path at all, so there is
no second rule for them.

The keys whose value is a path are listed in PATH_KEYS, in YamlPathUtils.cpp.
*/
namespace YamlInputReader {
    // A path leaf as the user wrote it, and the file that declared it. Kept so an
    // unresolvable path can be reported with the text the user typed, not only with
    // the resolved path they never saw.
    struct PathOrigin {
        std::string asWritten;
        std::string declaredIn;
    };
    // Keyed by the full key path, e.g. "SIMULATION MENU -> OUTPUT SUBMENU -> Output folder (string)".
    using PathOriginMap = std::map<std::string, PathOrigin>;

    // "unset" and "use the compiled-in data". Neither is a path, so neither is resolved
    // nor checked, and a menu reader that needs a real path has to reject them itself.
    bool isPathSentinel(const std::string& value);

    // Rewrite every path leaf of one input file to an absolute path, resolving a
    // relative one against baseDir, the directory of the file that declares it.
    // Records what each path looked like as written into origins.
    // yaml-cpp nodes are reference types, so the node is edited in place.
    void resolvePathsInPlace(YAML::Node node, const std::filesystem::path& baseDir,
                             const std::string& sourceFile, PathOriginMap& origins);

    // Report a path that names nothing on disk, once the merge is done.
    void checkPathsExist(const YAML::Node& merged, const PathOriginMap& origins);

    // Reject a PATH_KEYS entry that names no scalar leaf of the default tree, and any
    // default path that is relative. Neither can happen without a broken repository,
    // which is why this runs once per read rather than being a test-only check.
    void checkDefaultPaths(const YAML::Node& defaultData);

    // Read a path leaf that resolvePathsInPlace has already made absolute.
    // Returns the "=MISSING=" and "=DEFAULT=" sentinels unchanged.
    std::string readPath(const YAML::Node& node, const std::string& key);
}

#endif // YAMLPATHUTILS_H_
