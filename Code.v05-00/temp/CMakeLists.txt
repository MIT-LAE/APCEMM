cmake_minimum_required(VERSION 3.23.0)

project(APCEMM)

# Version
set(APCEMM_VERSION_NUMBER_MAJOR 5)
set(APCEMM_VERSION_NUMBER_MINOR 0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)



# Options
option(WITH_OPENMP "Build with OpenMP" OFF)


find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDES})

list(APPEND CMAKE_PREFIX_PATH "/home/vmeijer/local/lib/cmake/netCDF/")
list(APPEND CMAKE_MODULE_PATH "../")
find_package(NetCDF REQUIRED)
find_package (netCDFCxx REQUIRED) 
include_directories(${NETCDF_INCLUDES})
include_directories(${netCDFCxx_INCLUDES})

get_cmake_property(_variableNames VARIABLES)
list (SORT _variableNames)
foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
endforeach()

message(WARNING ${netCDFCxx_INCLUDES})
message(WARNING ${netCDFCxx_LIBRARIES})
add_executable(TEMP temp.cpp)
target_link_libraries(TEMP ${netCDFCxx_LIBRARIES} ${NETCDF_LIBRARIES})