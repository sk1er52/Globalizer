set(GLOBALIZER_PROJECT_NAME Globalizer)
cmake_minimum_required(VERSION 3.31)

cmake_policy(SET CMP0057 NEW)#adds support for the new IN_LIST operator

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Available build types" FORCE)

set(GLOBALIZER_MPI_LIB_USER_PATH "" CACHE STRING "Path to MPI libraries")
set(GLOBALIZER_MPI_INCLUDE_USER_PATH "" CACHE STRING "Path to MPI includes")
set(GLOBALIZER_MPI "ms" CACHE STRING "Set MPI implementation (ms, mpich, intel)")
set(GLOBALIZER_OUTPUT_LEVEL "L0" CACHE STRING "Set level of output (L0, L1, L2)")
option(GLOBALIZER_USE_MP "Set OFF to disable MPFR (extended precision in evolvent will unavaliable)" ON)
option(GLOBALIZER_USE_MPI "Set OFF to disable MPI.
  In this case a stub implementation of MPI will be used and globalizer won't support multi-process mode" ON)
option(GLOBALIZER_USE_CUDA "Set ON to use CUDA" OFF)
option(GLOBALIZER_BUILD_TESTS "Set ON to build tests" ON)
option(GLOBALIZER_CMAKE_VERBOSE "Set ON to show debug messages in CMake" ON)
option(GLOBALIZER_MAX_DIMENSION "Set MAX DIMENSION Globalizer" OFF)
option(GLOBALIZER_BUILD_PROBLEMS "Set ON to build problems" OFF)


option(GLOBALIZER_ONE_MPI_PROCESS_PER_NODE "ONE PROCESS PER NODE" OFF)

option(GLOBALIZER_PYTHON "Set ON to use Python" OFF)

set(GLOBALIZER_DEBUG_SUFFIX   d)

set(GLOBALIZER_EXECUTABLE_NAME                 ${GLOBALIZER_PROJECT_NAME}${GLOBALIZER_EXECUTABLE_SUFFIX})
set(GLOBALIZER_EXECUTABLE_NAME_DEBUG           ${GLOBALIZER_EXECUTABLE_NAME}_${GLOBALIZER_DEBUG_SUFFIX})
set(GLOBALIZER_EXECUTABLE_NAME_RELEASE         ${GLOBALIZER_EXECUTABLE_NAME})
set(GLOBALIZER_LIB_NAME                        ${GLOBALIZER_EXECUTABLE_NAME}_l)
set(GLOBALIZER_LIB_NAME_RELEASE                ${GLOBALIZER_EXECUTABLE_NAME}_l)
set(GLOBALIZER_LIB_NAME_DEBUG                  ${GLOBALIZER_EXECUTABLE_NAME}_l_${GLOBALIZER_DEBUG_SUFFIX})



set(CMAKE_CXX_STANDARD 17)
project(${GLOBALIZER_PROJECT_NAME})

set(GLOBALIZER_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/_bin CACHE STRING "Output directory")

message("GLOBALIZER_OUTPUT_DIRECTORY" ${GLOBALIZER_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${GLOBALIZER_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${GLOBALIZER_OUTPUT_DIRECTORY})

include(cmake/globalizer_utils.cmake)
include(cmake/globalizer_mpi.cmake)
#include(cmake/globalizer_lib.cmake)
include(cmake/globalizer_doxygen.cmake)
GLOBALIZER_initialize()

if(${GLOBALIZER_MAX_DIMENSION})
	message("true " ${GLOBALIZER_MAX_DIMENSION})
	add_definitions(-DGLOBALIZER_MAX_DIMENSION=${GLOBALIZER_MAX_DIMENSION})
else()
	message("false" ${GLOBALIZER_MAX_DIMENSION})
	add_definitions(-DGLOBALIZER_MAX_DIMENSION=10)
endif()

if (IS_GNU_COMPILER)
  include(cmake/GMPConfig.cmake)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
  find_program(CCACHE_FOUND ccache)
  if(CCACHE_FOUND)
	set(CMAKE_CXX_COMPILER_LAUNCHER ccache)
  endif(CCACHE_FOUND)
elseif (IS_MICROSOFT_COMPILER)
  add_definitions(/EHsc /openmp -D_CRT_SECURE_NO_WARNINGS -DNOMINMAX)
  replace_compiler_flag("W3" "W4")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Ox /wd4100")
  if (GLOBALIZER_TARGET_ARCH MATCHES "x86")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /arch:SSE2")
  endif()
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /wd4100")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /ignore:4099")
else()
  set(OpenMP_C_FLAGS "-fopenmp")
  set(OpenMP_CXX_FLAGS "-fopenmp")
endif()


if(NOT (CMAKE_BUILD_TYPE IN_LIST CMAKE_CONFIGURATION_TYPES))
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "build type" FORCE)
  message(STATUS "Unknown build type has been changed to " ${CMAKE_BUILD_TYPE})
endif()


init_doxygen()
setup_mpi(GLOBALIZER_USE_MPI GLOBALIZER_MPI_LIB GLOBALIZER_MPI_INCLUDES)

foreach(includeFolder ${GLOBALIZER_MPI_INCLUDES})
  if (NOT EXISTS ${includeFolder})
    message(WARNING "MPI include directories is not found: " ${includeFolder})
  endif()
endforeach()

if (GLOBALIZER_MPI_LIB STREQUAL "")
  message(WARNING "MPI library is not found")
endif()




if (${GLOBALIZER_ONE_MPI_PROCESS_PER_NODE})
	add_definitions(-DONE_MPI_PROCESS_PER_NODE)
endif()



if (${GLOBALIZER_USE_MP})
  setup_mpir(GLOBALIZER_MPIR_LIB_DIR GLOBALIZER_MPIR_INCLUDES GLOBALIZER_MPIR_LIB)
  setup_mpfr(GLOBALIZER_MPFR_LIB_DIR GLOBALIZER_MPFR_INCLUDES GLOBALIZER_MPFR_LIB)
  setup_mpreal(GLOBALIZER_MPREAL_INCLUDES)
else()
  add_definitions(-DNO_MPREAL)
endif()

# set include directories for project
set(GLOBALIZER_EXAMPLES_SRC_DIR ${PROJECT_SOURCE_DIR}/examples)
set(GLOBALIZER_COMMON_SRC_DIR ${PROJECT_SOURCE_DIR}/globalizer)
set(GLOBALIZER_PARAMETERS_DIR ${GLOBALIZER_COMMON_SRC_DIR}/parameters)
set(GLOBALIZER_evolvent_DIR ${GLOBALIZER_COMMON_SRC_DIR}/evolvent)
set(GLOBALIZER_extended_DIR ${GLOBALIZER_COMMON_SRC_DIR}/extended)
set(GLOBALIZER_output_system_DIR ${GLOBALIZER_COMMON_SRC_DIR}/output_system)
set(GLOBALIZER_problem_DIR ${GLOBALIZER_COMMON_SRC_DIR}/problem)
set(GLOBALIZER_method_DIR ${GLOBALIZER_COMMON_SRC_DIR}/method)
set(GLOBALIZER_calculation_DIR ${GLOBALIZER_COMMON_SRC_DIR}/method/calculation)
set(GLOBALIZER_local_method_DIR ${GLOBALIZER_COMMON_SRC_DIR}/method/local_method)


include_directories(
  ${GLOBALIZER_MPI_INCLUDES}
  ${GLOBALIZER_MPIR_INCLUDES}
  ${GLOBALIZER_MPFR_INCLUDES}
  ${GLOBALIZER_MPREAL_INCLUDES}
  
  ${GLOBALIZER_COMMON_SRC_DIR}/include
  ${GLOBALIZER_PARAMETERS_DIR}/include
  ${GLOBALIZER_evolvent_DIR}/include
  ${GLOBALIZER_extended_DIR}/include
  ${GLOBALIZER_output_system_DIR}/include
  ${GLOBALIZER_problem_DIR}/include
  ${GLOBALIZER_method_DIR}/include
  ${GLOBALIZER_calculation_DIR}/include
  ${GLOBALIZER_local_method_DIR}/include
)

link_directories(
  ${GLOBALIZER_MPIR_LIB_DIR}
  ${GLOBALIZER_MPFR_LIB_DIR}
)

file(GLOB ROOT_HEADERS ${GLOBALIZER_COMMON_SRC_DIR}/include/*.h ${GLOBALIZER_COMMON_SRC_DIR}/include/*.hpp )
file(GLOB evolvent_HEADERS ${GLOBALIZER_evolvent_DIR}/include/*.h ${GLOBALIZER_evolvent_DIR}/include/*.hpp)
file(GLOB extended_HEADERS ${GLOBALIZER_extended_DIR}/include/*.h ${GLOBALIZER_extended_DIR}/include/*.hpp)
file(GLOB output_system_HEADERS ${GLOBALIZER_output_system_DIR}/include/*.h ${GLOBALIZER_output_system_DIR}/include/*.hpp)
file(GLOB problem_HEADERS ${GLOBALIZER_problem_DIR}/include/*.h ${GLOBALIZER_problem_DIR}/include/*.hpp)
file(GLOB method_HEADERS ${GLOBALIZER_method_DIR}/include/*.h ${GLOBALIZER_method_DIR}/include/*.hpp)
file(GLOB calculation_HEADERS ${GLOBALIZER_calculation_DIR}/include/*.h ${GLOBALIZER_calculation_DIR}/include/*.hpp)
file(GLOB local_method_HEADERS ${GLOBALIZER_local_method_DIR}/include/*.h ${GLOBALIZER_local_method_DIR}/include/*.hpp)						
file(GLOB PARAMETERS_HEADERS ${GLOBALIZER_PARAMETERS_DIR}/include/*.h ${GLOBALIZER_PARAMETERS_DIR}/include/*.hpp)

file(GLOB ROOT_SOURCES ${GLOBALIZER_COMMON_SRC_DIR}/src/*.c ${GLOBALIZER_COMMON_SRC_DIR}/src/*.cpp)
file(GLOB evolvent_SOURCES ${GLOBALIZER_evolvent_DIR}/src/*.c ${GLOBALIZER_evolvent_DIR}/src/*.cpp)
file(GLOB extended_SOURCES ${GLOBALIZER_extended_DIR}/src/*.c ${GLOBALIZER_extended_DIR}/src/*.cpp)
file(GLOB output_system_SOURCES ${GLOBALIZER_output_system_DIR}/src/*.c ${GLOBALIZER_output_system_DIR}/src/*.cpp)
file(GLOB problem_SOURCES ${GLOBALIZER_problem_DIR}/src/*.c ${GLOBALIZER_problem_DIR}/src/*.cpp)
file(GLOB method_SOURCES ${GLOBALIZER_method_DIR}/src/*.c ${GLOBALIZER_method_DIR}/src/*.cpp)
file(GLOB calculation_SOURCES ${GLOBALIZER_calculation_DIR}/src/*.c ${GLOBALIZER_calculation_DIR}/src/*.cpp)
file(GLOB local_method_SOURCES ${GLOBALIZER_local_method_DIR}/src/*.c ${GLOBALIZER_local_method_DIR}/src/*.cpp)
file(GLOB PARAMETERS_SOURCES ${GLOBALIZER_PARAMETERS_DIR}/src/*.c ${GLOBALIZER_PARAMETERS_DIR}/src/*.cpp)


if (GLOBALIZER_USE_CUDA)
  find_package(CUDA REQUIRED)
endif()

# build main library and executable
set(ALL_SRC_FILES ${ROOT_HEADERS} ${PARAMETERS_HEADERS} ${ROOT_SOURCES}
${evolvent_HEADERS} 
${extended_HEADERS}
${output_system_HEADERS}
${problem_HEADERS}
${method_HEADERS}
${calculation_HEADERS}
${local_method_HEADERS}
${evolvent_SOURCES}
${extended_SOURCES}
${output_system_SOURCES}
${problem_SOURCES}
${method_SOURCES}
${calculation_SOURCES}
${local_method_SOURCES}
${PARAMETERS_SOURCES}
)
list(REMOVE_ITEM ALL_SRC_FILES ${GLOBALIZER_COMMON_SRC_DIR}/Main.cpp)

add_library(${GLOBALIZER_LIB_NAME} STATIC ${ALL_SRC_FILES})
add_executable(${GLOBALIZER_PROJECT_NAME} ${GLOBALIZER_EXAMPLES_SRC_DIR}/Main.cpp)


if (${GLOBALIZER_PYTHON})
  find_package (Python COMPONENTS Interpreter Development NumPy)

#set(EXTRA_DEFINES "")
#if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  #set(EXTRA_DEFINES PYLIBNAME="${Python_LIBRARIES}")
  message("Python pname ${Python_LIBRARIES}")
  message("Python end")
  #message(${PYLIBNAME})
  add_definitions(-DUSE_PYTHON)
#endif()

  if (Python_Development_FOUND)
    message("Python Development FOUND")
  endif()
  
  if (Python_NumPy_FOUND)
    message("Python NumPy FOUND")
  endif()

else()
  message("Python not use")
  set(Python_INCLUDE_DIRS "")
  set(Python_NumPy_INCLUDE_DIRS "")
  set(Python_LIBRARIES "")
  
endif()

target_include_directories(${GLOBALIZER_PROJECT_NAME} PRIVATE ${GLOBALIZER_COMMON_SRC_DIR}
												  PRIVATE ${Python_INCLUDE_DIRS} 
												  PRIVATE ${Python_NumPy_INCLUDE_DIRS})
												  
include_directories(${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS})

	
if (${IS_CLANG_COMPILER})	
  target_link_libraries(${GLOBALIZER_PROJECT_NAME} ${GLOBALIZER_LIB_NAME} ${Python_LIBRARIES} OpenMP::OpenMP_CXX  )
  target_link_libraries(${GLOBALIZER_LIB_NAME}
    ${GLOBALIZER_MPI_LIB}
    ${GLOBALIZER_MPFR_LIB}
    ${GLOBALIZER_MPIR_LIB}
    ${DISLIN_LIB}
	OpenMP::OpenMP_CXX  )
else()
  target_link_libraries(${GLOBALIZER_PROJECT_NAME} ${GLOBALIZER_LIB_NAME} ${Python_LIBRARIES})
  target_link_libraries(${GLOBALIZER_LIB_NAME}
    ${GLOBALIZER_MPI_LIB}
    ${GLOBALIZER_MPFR_LIB}
    ${GLOBALIZER_MPIR_LIB}
    ${DISLIN_LIB})

endif()



message("MPI lib: " ${GLOBALIZER_LIB_NAME})

set_target_properties(${GLOBALIZER_LIB_NAME} PROPERTIES OUTPUT_NAME_DEBUG ${GLOBALIZER_LIB_NAME_DEBUG})
set_target_properties(${GLOBALIZER_LIB_NAME} PROPERTIES OUTPUT_NAME_RELEASE ${GLOBALIZER_LIB_NAME_RELEASE})

set_target_properties(${GLOBALIZER_PROJECT_NAME} PROPERTIES OUTPUT_NAME_DEBUG ${GLOBALIZER_EXECUTABLE_NAME_DEBUG})
set_target_properties(${GLOBALIZER_PROJECT_NAME} PROPERTIES OUTPUT_NAME_RELEASE ${GLOBALIZER_EXECUTABLE_NAME_RELEASE})

message("Find OpenMP") 
find_package(OpenMP)
if (OPENMP_FOUND)
  message("OpenMP found")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
  message("OpenMP not found")
endif()

# add linker flags for GCC
if (IS_GNU_COMPILER)
  target_compile_options(${GLOBALIZER_LIB_NAME} INTERFACE -march=native -fopenmp)
  set_target_properties(${GLOBALIZER_PROJECT_NAME} PROPERTIES LINK_FLAGS "-Wl,--no-as-needed -ldl -fopenmp")
endif()

#add problems
if(${GLOBALIZER_BUILD_PROBLEMS})
	add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/Problems)
endif()

#add tests
if(${GLOBALIZER_BUILD_TESTS})
  enable_testing()
  add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/gtest)
  add_subdirectory(${PROJECT_SOURCE_DIR}/tests)
endif()




if(${GLOBALIZER_CMAKE_VERBOSE})
  print_build_config()
endif()


# create source groups for VS IDE projects
if (IS_MICROSOFT_COMPILER)
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
  #set_property(TARGET pugixml PROPERTY FOLDER "3rd-party")
  if(${GLOBALIZER_BUILD_TESTS})
    set_property(TARGET gtest gmock gmock_main gtest_main mpi_stub PROPERTY FOLDER "3rd-party")
  endif()


  source_group("Header Files" FILES ${ROOT_HEADERS})
  source_group("Header Files\\Parameters" FILES ${PARAMETERS_HEADERS})
  source_group("Header Files\\evolvent" FILES ${evolvent_HEADERS}) 
  source_group("Header Files\\extended" FILES ${extended_HEADERS})
  source_group("Header Files\\output_system" FILES ${output_system_HEADERS})
  source_group("Header Files\\problem" FILES ${problem_HEADERS})
  source_group("Header Files\\method" FILES ${method_HEADERS})
  source_group("Header Files\\calculation" FILES ${calculation_HEADERS})
  source_group("Header Files\\local_method" FILES ${local_method_HEADERS})

  source_group("Source Files" FILES ${ROOT_SOURCES})
  source_group("Source Files\\evolvent" FILES ${evolvent_SOURCES})
  source_group("Source Files\\extended" FILES ${extended_SOURCES})
  source_group("Source Files\\output_system" FILES ${output_system_SOURCES})
  source_group("Source Files\\problem" FILES ${problem_SOURCES})
  source_group("Source Files\\method" FILES ${method_SOURCES})
  source_group("Source Files\\calculation" FILES ${calculation_SOURCES})
  source_group("Source Files\\local_method" FILES ${local_method_SOURCES})
  source_group("Source Files\\parameters" FILES ${PARAMETERS_SOURCES})
  
  
  
endif()
