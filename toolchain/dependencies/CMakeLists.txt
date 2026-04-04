# This is a superbuild CMake file
# See: https://cmake.org/cmake/help/latest/module/ExternalProject.html.

cmake_minimum_required(VERSION 3.12)

project(figr-Dependencies LANGUAGES C CXX Fortran)

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/regular")
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
    list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/cce")
endif()

set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)

# Imports
include(GNUInstallDirs)
include(ExternalProject)


# Options
option(FIGR_HDF5    "Build the HDF5  dependency" OFF)
option(FIGR_SILO    "Build the SILO  dependency" OFF)


# HDF5
if (FIGR_HDF5)
    find_package(HDF5)
    if (HDF5_FOUND)
        message(STATUS "HDF5 found.")
        add_custom_target(hdf5)
    else()
        ExternalProject_Add(hdf5
            GIT_REPOSITORY "https://github.com/HDFGroup/hdf5"
            GIT_TAG        hdf5-1_12_2
            GIT_SHALLOW    ON
            GIT_PROGRESS   ON
            CMAKE_ARGS     "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
                            -DBUILD_SHARED_LIBS=OFF
                            -DBUILD_TESTING=OFF
                            -DHDF5_BUILD_UTILS=OFF
                            -DHDF5_BUILD_TOOLS=ON
                            -DHDF5_BUILD_EXAMPLES=OFF
                            -DHDF5_BUILD_CPP_LIB=OFF
                            -DHDF5_BUILD_FORTRAN=OFF
                            -DHDF5_ENABLE_PARALLEL=OFF
        )
    endif()
endif()


# SILO
if (FIGR_SILO)
    find_package(SILO)
    if (SILO_FOUND)
        message(STATUS "SILO found.")
        add_custom_target(silo)
    else()
        find_package(Git REQUIRED)

        ExternalProject_Add(silo
            GIT_REPOSITORY "https://github.com/LLNL/Silo"
            GIT_TAG        4.12.0
            GIT_SHALLOW    ON
            GIT_PROGRESS   ON
            PATCH_COMMAND  "${GIT_EXECUTABLE}" stash
                        && "${GIT_EXECUTABLE}" apply "${CMAKE_SOURCE_DIR}/Silo.patch"
            CMAKE_ARGS     -DSILO_ENABLE_SHARED=OFF
                           -DSILO_ENABLE_SILOCK=OFF
                           -DSILO_ENABLE_BROWSER=OFF
                           -DSILO_ENABLE_ZFP=OFF
                           -DFIND_LIBRARY_USE_LIB64_PATHS=ON
                           "-DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}"
                           "-DCMAKE_FIND_ROOT_PATH=${CMAKE_FIND_ROOT_PATH}"
                           "-DCMAKE_FIND_PACKAGE_REDIRECTS_DIR=${CMAKE_FIND_PACKAGE_REDIRECTS_DIR}"
                           "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
                           "$<$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},Cray>:-DCMAKE_MODULE_PATH=${CMAKE_SOURCE_DIR}/../cmake/cce>"
        )
    endif()
endif()
