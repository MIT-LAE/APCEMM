include(Catch)

set(TARGET_TEST unittest)

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

add_executable(
	${TARGET_TEST}
	${SRC_TEST}
)
target_link_libraries(${TARGET_TEST}  Catch2::Catch2WithMain yaml-cpp Util AIM EPM YamlInputReader)
catch_discover_tests(${TARGET_TEST})

set(TARGET_TEST test_solver)

set(SRC_TEST
	test_adv_diff_solver.cpp
)

add_executable(
	${TARGET_TEST}
	${SRC_TEST}
)
target_link_libraries(${TARGET_TEST}  Catch2::Catch2WithMain FVM_ANDS Core)
catch_discover_tests(${TARGET_TEST})

set(TARGET_TEST test_LAGRID)

set(SRC_TEST
	test_LAGRID.cpp
)

add_executable(
	${TARGET_TEST}
	${SRC_TEST}
)
target_link_libraries(${TARGET_TEST}  Catch2::Catch2WithMain LAGRID)
catch_discover_tests(${TARGET_TEST})