include(Catch)

set(TARGET_TEST unittest)

set(SRC_TEST
	testmain.cpp
    test_physfunction.cpp
    test_nucleation.cpp
    test_buildkernel.cpp
    test_aerosol.cpp
    test_integrate.cpp
    test_metfunction.cpp
    test_aircraft.cpp
    test_yamlreader.cpp
)
#Add preprocessor def of the tests dir
add_definitions(-DAPCEMM_TESTS_DIR="${CMAKE_SOURCE_DIR}/tests")

add_executable(
	${TARGET_TEST}
	${SRC_TEST}
)
target_link_libraries(${TARGET_TEST}  Catch2::Catch2WithMain yaml-cpp Util AIM EPM YamlInputReader)
catch_discover_tests(${TARGET_TEST})