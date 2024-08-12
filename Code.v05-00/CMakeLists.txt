cmake_minimum_required(VERSION 3.16)

set(CMAKE_TOOLCHAIN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/submodules/vcpkg/scripts/buildsystems/vcpkg.cmake"
  CACHE STRING "Vcpkg toolchain file")

project(APCEMM)

# Options
option(DEBUG "Enable debug mode" OFF)
SET(RINGS 0)
SET(OMP 1)
option(BUILD_TEST ON)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")	
    include(CheckCXXCompilerFlag)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -gdwarf-2 -gstrict-dwarf -O0 -Wall -Wextra -Wconversion")
    #add_compile_options(-fno-omit-frame-pointer -fsanitize=bounds -fsanitize=address -fsanitize=undefined)
    #add_link_options(-fsanitize=bounds -fsanitize=address -fsanitize=undefined)
    add_compile_options(-fno-omit-frame-pointer -fsanitize=bounds -fsanitize=undefined)
    add_link_options(-fsanitize=bounds -fsanitize=undefined)
elseif (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
    if(COMPILER_SUPPORTS_MARCH_NATIVE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    endif()
    
    # Turn on warnings always
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp")


    #Enabling unsafe optimizations doesn't seem to affect the runtime much, but could be something to check just in case later.
    #-Michael
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -fno-trapping-math -funsafe-math-optimizations")
endif()

execute_process(COMMAND git rev-parse --short HEAD 
                    OUTPUT_VARIABLE APCEMM_VERSION_BUILD_NUMBER 
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

set(APCEMM_VERSION_BUILD_NUMBER "${APCEMM_VERSION_BUILD_NUMBER}" CACHE STRING "Git Version/Commit Hash")
set(CMAKE_CXX_STANDARD 20)
configure_file(${CMAKE_SOURCE_DIR}/include/APCEMM.h.in ${CMAKE_SOURCE_DIR}/include/APCEMM.h @ONLY)

# Defines the 'main' file we want to create a binary of
add_executable(${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src/Core/Main.cpp)


# Dependencies managed by vcpkg:

find_package(Boost REQUIRED)
# All the Boost libraries we use are header-only, so no need to link here.

find_package(Catch2 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE Catch2::Catch2 Catch2::Catch2WithMain)

find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)

find_package(fmt CONFIG REQUIRED)
# target_link_libraries(${PROJECT_NAME} PRIVATE fmt::fmt)
target_link_libraries(${PROJECT_NAME} PRIVATE fmt::fmt-header-only)

find_package(FFTW3 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE FFTW3::fftw3)

find_package(netCDFCxx CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE netCDF::netcdf netCDF::netcdf-cxx4)

find_package(yaml-cpp CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE yaml-cpp::yaml-cpp)

# Other dependencies:

find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)


# This ensures the header files of the necessary libraries are included
# SYSTEM flag silences the compiler warnings from external dependencies 
include_directories(SYSTEM ${Boost_INCLUDES})
include_directories(SYSTEM ${OpenMP_INCLUDES})
include_directories(SYSTEM ${FFTW_INCLUDES})
include_directories(SYSTEM ${NETCDF_INCLUDES})
include_directories(SYSTEM ${netCDFCxx_INCLUDE_DIRS})
include_directories(SYSTEM ${YAML_CPP_INCLUDE_DIR})

# Ensure header files of APCEMM are included
include_directories(${CMAKE_SOURCE_DIR}/include)

# This ensures that the CMakeLists.txt in src/ gets build
add_subdirectory(${CMAKE_SOURCE_DIR}/src)

# Defines the libraries necessary for compiling the executable
target_link_libraries(${PROJECT_NAME} PRIVATE
                      FVM_ANDS Util AIM KPP EPM Core YamlInputReader LAGRID)

# Tests
#if (BUILD_TEST)
include(CTest)
enable_testing()
add_subdirectory(tests)
#endif()

if (DEBUG)
  message(STATUS "DEBUG mode is enabled")
endif()
