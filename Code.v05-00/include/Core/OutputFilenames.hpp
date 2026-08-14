#ifndef OUTPUTFILENAMES_H_INCLUDED
#define OUTPUTFILENAMES_H_INCLUDED
#include <string>

// Default files names that cannot be changed by
// the user through the input.yaml

namespace OutputFiles {
inline const std::string README = "README";
inline const std::string STATUS = "status";
inline const std::string MERGED_YAML = "merged-input.yaml";
inline const std::string EPM_OUTPUT = "epm-output.nc";
inline const std::string MICROPHYSICS = "Micro.out";
} // namespace OutputFiles

#endif
