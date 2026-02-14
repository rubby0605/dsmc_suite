cmake_minimum_required(VERSION 3.12)
project(dsmc_suite C)

set(CMAKE_C_STANDARD 11)

# Optimization flags
set(CMAKE_C_FLAGS_RELEASE "-O3 -Wall -Wextra -Wno-unused-parameter")
set(CMAKE_C_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wno-unused-parameter -DDEBUG")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# OpenMP setup
if(APPLE)
    # Apple Clang doesn't ship OpenMP — use Homebrew libomp
    # brew install libomp
    execute_process(
        COMMAND brew --prefix libomp
        OUTPUT_VARIABLE LIBOMP_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(LIBOMP_PREFIX)
        message(STATUS "Found Homebrew libomp at: ${LIBOMP_PREFIX}")
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
        set(OpenMP_C_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")
        include_directories("${LIBOMP_PREFIX}/include")
        link_directories("${LIBOMP_PREFIX}/lib")
    else()
        message(FATAL_ERROR
            "libomp not found. Install with: brew install libomp")
    endif()
    find_package(OpenMP REQUIRED)
else()
    find_package(OpenMP REQUIRED)
endif()

# Embree ray tracing (brew install embree)
if(APPLE)
    execute_process(
        COMMAND brew --prefix embree
        OUTPUT_VARIABLE EMBREE_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(EMBREE_PREFIX)
        message(STATUS "Found Homebrew Embree at: ${EMBREE_PREFIX}")
        set(EMBREE_INCLUDE_DIR "${EMBREE_PREFIX}/include")
        set(EMBREE_LIBRARY "${EMBREE_PREFIX}/lib/libembree4.dylib")
    else()
        message(FATAL_ERROR "Embree not found. Install with: brew install embree")
    endif()
else()
    find_package(embree 4 REQUIRED)
    set(EMBREE_INCLUDE_DIR ${embree_INCLUDE_DIRS})
    set(EMBREE_LIBRARY embree)
endif()

# Include dirs
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${EMBREE_INCLUDE_DIR})

# Math library
find_library(MATH_LIBRARY m)

# Common library (shared modules)
add_library(dsmc_common STATIC
    src/common/grid.c
    src/common/velocity.c
    src/common/integrator.c
    src/common/collision.c
    src/common/polyfit.c
    src/common/fileio.c
)
target_link_libraries(dsmc_common PUBLIC OpenMP::OpenMP_C)
if(MATH_LIBRARY)
    target_link_libraries(dsmc_common PUBLIC ${MATH_LIBRARY})
endif()

# Main executable
add_executable(dsmc_suite
    src/main.c
    src/dsmc/ballistic.c
    src/dsmc/advanced.c
    src/raytracer/shadow_maker.c
    src/thermal/thermal_model.c
)
target_link_libraries(dsmc_suite PRIVATE dsmc_common ${EMBREE_LIBRARY})

# Test executable
add_executable(test_verlet tests/test_verlet.c)
target_link_libraries(test_verlet PRIVATE dsmc_common)

# Install
install(TARGETS dsmc_suite DESTINATION bin)
