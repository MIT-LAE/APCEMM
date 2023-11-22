cmake_minimum_required(VERSION 3.16)

project(APCEMM)

# Options
SET(DEBUG 0)
SET(RINGS 0)
SET(OMP 1)
option(BUILD_TEST ON)

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")	
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
    if(COMPILER_SUPPORTS_MARCH_NATIVE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    endif()
    
    #Enabling unsafe optimizations doesn't seem to affect the runtime much, but could be something to check just in case later.
    #-Michael
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -fno-trapping-math -funsafe-math-optimizations")
endif()

execute_process(COMMAND git rev-parse --short HEAD 
                    OUTPUT_VARIABLE APCEMM_VERSION_BUILD_NUMBER 
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

set(APCEMM_VERSION_BUILD_NUMBER "${APCEMM_VERSION_BUILD_NUMBER}" CACHE STRING "Git Version/Commit Hash")
set(CMAKE_CXX_STANDARD 17)
configure_file(${CMAKE_SOURCE_DIR}/include/APCEMM.h.in ${CMAKE_SOURCE_DIR}/include/APCEMM.h @ONLY)

# Include .cmake files in this directory to find the required packages
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Automated discovert of netCDF using netCDF config
# Find nc-config and ncxx4-config
find_program(NC_CONFIG NAMES "nc-config" DOC "Location of nc-config utility")
find_program(NCXX_CONFIG NAMES "ncxx4-config" DOC "Location of ncxx4-config utility")

# A function to call nx-config with an argument, and append the resulting path to a list
function(inspect_netcdf_config VAR NX_CONFIG ARG)
    execute_process(
        COMMAND ${NX_CONFIG} ${ARG}
        OUTPUT_VARIABLE NX_CONFIG_OUTPUT
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(EXISTS "${NX_CONFIG_OUTPUT}")
        list(APPEND ${VAR} ${NX_CONFIG_OUTPUT})
        set(${VAR} ${${VAR}} PARENT_SCOPE)
    endif()
endfunction()

# Add nc-config and ncxx4-config paths to CMAKE_PREFIX_PATH
inspect_netcdf_config(CMAKE_PREFIX_PATH "${NC_CONFIG}" "--prefix")
inspect_netcdf_config(CMAKE_PREFIX_PATH "${NCXX_CONFIG}" "--prefix")

# Find packages
find_package(NetCDF REQUIRED)
find_package(netCDFCxx REQUIRED)
find_package(FFTW REQUIRED)
find_package(OpenMP REQUIRED)
#find_package(Catch2 3 REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(Eigen3 3.3 REQUIRED)

# Try to get Catch2 automatically
Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.0.1 # or a later release
#  FIND_PACKAGE_ARGS # No specific args - but this ensure find_package runs first
)
#TODO: Figure out why FIND_PACKAGE_ARGS causes chaos

#FetchContent_Declare(
#        yaml-cpp
#        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
#        GIT_TAG 4aad2b1666a4742743b04e765a34742512915674
#)

#set(YAML_CPP_BUILD_TESTS OFF)
# This should attempt to find the package, and will download if not found
FetchContent_MakeAvailable(Catch2)
#FetchContent_MakeAvailable(Catch2 yaml-cpp)

# This ensures the header files of the necessary libraries are included
include_directories(${OpenMP_INCLUDES})
include_directories(${FFTW_INCLUDES})
include_directories(${NETCDF_INCLUDES})
include_directories(${netCDFCxx_INCLUDE_DIRS})
include_directories(${YAML_CPP_INCLUDE_DIR})

# Ensure header files of APCEMM are included
include_directories(${CMAKE_SOURCE_DIR}/include)

# This ensures that the CMakeLists.txt in src/ gets build
add_subdirectory(${CMAKE_SOURCE_DIR}/src)

# Defines the 'main' file we want to create a binary of 
add_executable(${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src/Core/Main.cpp)

# Defines the libraries necessary for compiling the executable
target_link_libraries(${PROJECT_NAME}
                      ${netCDFCxx_LIBRARIES} ${NETCDF_LIBRARIES} ${FFTW_LIBRARIES} ${YAML_CPP_LIBRARIES}
                      OpenMP::OpenMP_CXX Eigen3::Eigen FVM_ANDS Util AIM KPP EPM Core YamlInputReader LAGRID)

# Tests
#if (BUILD_TEST)
include(CTest)
enable_testing()
add_subdirectory(tests)
#endif()

