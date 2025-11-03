set(PROJECT_NAME_STR GLOBALIZER_tests)
project (${PROJECT_NAME_STR})

set(TESTS_OUTPUT_DIRECTORY ${GLOBALIZER_OUTPUT_DIRECTORY})

include_directories(${GLOBALIZER_COMMON_SRC_DIR})
file(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h)

set(TESTDATA_BIN_PATH ${TESTS_OUTPUT_DIRECTORY})
set(TESTDATA_PATH ${CMAKE_CURRENT_SOURCE_DIR}/data)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/templates/gtest_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config/test_config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/config)

include_directories(${PROJECT_SOURCE_DIR}/third_party/gtest/googletest/include)

add_executable(${PROJECT_NAME_STR} ${SRC_FILES})
set_target_properties(${PROJECT_NAME_STR} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TESTS_OUTPUT_DIRECTORY})
target_link_libraries(${PROJECT_NAME_STR}
  ${GLOBALIZER_LIB_NAME}
  ${GLOBALIZER_MPI_LIB}
  ${GLOBALIZER_MPFR_LIB}
  ${GLOBALIZER_MPIR_LIB}
  gtest)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set_target_properties(${PROJECT_NAME_STR} PROPERTIES LINK_FLAGS "-lm -Wl,--no-as-needed -ldl -lpthread")
endif()

add_test(all_tests ${TESTS_OUTPUT_DIRECTORY}/${PROJECT_NAME_STR})