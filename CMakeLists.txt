cmake_minimum_required(VERSION 3.20)
project(langevin_bench LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    add_compile_options(/O2 /W3 /MP /EHsc /arch:AVX2)
    add_compile_definitions(_USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS EIGEN_NO_DEBUG)
else()
    add_compile_options(-O3 -march=native -Wall -Wextra)
    add_compile_definitions(EIGEN_NO_DEBUG)
endif()

include(FetchContent)

FetchContent_Declare(eigen3_dl
    URL      https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    URL_HASH SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_GetProperties(eigen3_dl)
if(NOT eigen3_dl_POPULATED)
    FetchContent_Populate(eigen3_dl)
endif()
add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${eigen3_dl_SOURCE_DIR}"
)

FetchContent_Declare(nlopt
    GIT_REPOSITORY https://github.com/stevengj/nlopt.git
    GIT_TAG        v2.7.1
    GIT_SHALLOW    TRUE
)
set(NLOPT_PYTHON      OFF CACHE BOOL "" FORCE)
set(NLOPT_GUILE       OFF CACHE BOOL "" FORCE)
set(NLOPT_SWIG        OFF CACHE BOOL "" FORCE)
set(NLOPT_OCTAVE      OFF CACHE BOOL "" FORCE)
set(NLOPT_MATLAB      OFF CACHE BOOL "" FORCE)
set(NLOPT_CXX         OFF CACHE BOOL "" FORCE)
set(NLOPT_TESTS       OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(nlopt)
if(NOT TARGET NLopt::nlopt)
    add_library(NLopt::nlopt ALIAS nlopt)
endif()

add_library(gsl_mini STATIC gsl_impl.cpp)
target_include_directories(gsl_mini PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(MSVC)
    target_compile_definitions(gsl_mini PRIVATE _USE_MATH_DEFINES)
endif()

add_library(langevin_core STATIC langevin_core.cpp)
target_include_directories(langevin_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(langevin_core PUBLIC LANGEVIN_STATIC)
if(MSVC)
    target_compile_definitions(langevin_core PRIVATE _USE_MATH_DEFINES)
endif()
target_link_libraries(langevin_core
    PUBLIC
        gsl_mini
        NLopt::nlopt
        Eigen3::Eigen
)

