include(Catch)

set(SRC_TEST
	testmain.cpp
    test_physfunction.cpp
    test_nucleation.cpp
    test_buildkernel.cpp
    test_aerosol.cpp
	#test_meteorology.cpp
    test_integrate.cpp
    test_metfunction.cpp
    test_aircraft.cpp
    test_yamlreader.cpp
)
#Add preprocessor def of the tests dir
add_definitions(-DAPCEMM_TESTS_DIR="${CMAKE_SOURCE_DIR}/tests")

add_executable(unittest ${SRC_TEST})
target_link_libraries(unittest  Catch2::Catch2WithMain Util AIM EPM YamlInputReader)
catch_discover_tests(unittest)

add_executable(test_solver test_adv_diff_solver.cpp)
target_link_libraries(test_solver  Catch2::Catch2WithMain FVM_ANDS Core)
catch_discover_tests(test_solver)

add_executable(test_LAGRID test_LAGRID.cpp)
target_link_libraries(test_LAGRID  Catch2::Catch2WithMain LAGRID)
catch_discover_tests(test_LAGRID)
