set(SRCS
    YamlInputReader.cpp)
add_library(YamlInputReader STATIC ${SRCS})
target_link_libraries(YamlInputReader yaml-cpp::yaml-cpp)

