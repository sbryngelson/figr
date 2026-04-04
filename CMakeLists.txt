cmake_minimum_required(VERSION 3.20)

project(figr LANGUAGES C CXX Fortran)

option(FIGR_OpenACC       "Build with OpenACC"                                OFF)
option(FIGR_OpenMP        "Build with OpenMP"                                 OFF)
option(FIGR_PRE_PROCESS   "Build pre_process"                                 OFF)
option(FIGR_SIMULATION    "Build simulation"                                  OFF)
option(FIGR_POST_PROCESS  "Build post_process"                                OFF)
option(FIGR_SYSCHECK      "Build syscheck"                                    OFF)
option(FIGR_ALL           "Build everything"                                  OFF)
option(FIGR_SINGLE_PRECISION "Build single precision"                         OFF)
option(FIGR_MIXED_PRECISION "Build mixed precision"                           OFF)

if (FIGR_ALL)
    set(FIGR_PRE_PROCESS   ON FORCE)
    set(FIGR_SIMULATION    ON FORCE)
    set(FIGR_POST_PROCESS  ON FORCE)
endif()

# Validate CMAKE_BUILD_TYPE to catch typos (CMake is case-sensitive).
set(_VALID_BUILD_TYPES "Debug" "Release" "RelDebug" "")
if (NOT CMAKE_BUILD_TYPE IN_LIST _VALID_BUILD_TYPES)
    message(FATAL_ERROR "Unknown CMAKE_BUILD_TYPE '${CMAKE_BUILD_TYPE}'. Valid: Debug, RelDebug, Release")
endif()

# RelDebug: a lighter debug mode for CI. Compiler-specific blocks below add the
# actual flags; these defaults just tell CMake it is a recognized build type.
set(CMAKE_C_FLAGS_RELDEBUG       "-g" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELDEBUG     "-g" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELDEBUG "-g" CACHE STRING "")

if (FIGR_SINGLE_PRECISION)
    add_compile_definitions(FIGR_SINGLE_PRECISION)
else()
    add_compile_definitions(FIGR_DOUBLE_PRECISION)
endif()

if (FIGR_MIXED_PRECISION)
    add_compile_definitions(FIGR_MIXED_PRECISION)
endif()


# CMake Library Imports

include(GNUInstallDirs)
include(ExternalProject)
include(CheckIPOSupported)
include(CMakeParseArguments)
include(CheckFortranCompilerFlag)


# Check Compiler Support: Some compilers, especially old ones, are known not to
# be able to compile figr. We throw our own error early - at configure time - in
# such cases.

set(__err_msg "\
CMake detected the ${CMAKE_Fortran_COMPILER_ID} Fortran compiler \
v${CMAKE_Fortran_COMPILER_VERSION}. If you intended to use a different \
compiler (or a different version thereof), please:\n\
 - Install the compiler or load its module. (e.g. module load gcc/10.1)\n\
 - Set/Export the C, CXX, and FC environment variables. (e.g. 'export CC=gcc', \
 'export CXX=g++', and 'export FC=gfortran'.\n\
 - If using figr.sh, delete the build/<code name> directory and try again. (e.g. 'rm -rf build/pre_process')")

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 5)
        message(FATAL_ERROR "ERROR: GNU v5.0 or newer is required to build figr.\n${__err_msg}")
    endif()
    if (FIGR_OpenACC)
        message(FATAL_ERROR "ERROR: figr with GPU processing is not currently compatible with GNU compilers. Please use NVIDIA or Cray compilers.\n${__err_msg}")
    endif()
elseif ((CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC") OR (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI"))
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 21.7)
        message(FATAL_ERROR "ERROR: When using NVHPC, v21.7 or newer is required to build figr.\n${__err_msg}")
    endif()

    if ((CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 24.5) AND (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelDebug") AND FIGR_OpenACC)
        message(FATAL_ERROR "ERROR: When using NVHPC, figr with Debug and GPU options requires NVHPC v24.5 or newer.\n${__err_msg}")
    endif()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_C_COMPILER_ID STREQUAL "AppleClang")
    message(FATAL_ERROR "ERROR: figr does not support the Apple Clang compilers. Please consult the documentation.\n${__err_msg}")
endif()

# Fypp is required to build figr. We try to locate it. The path to the binary is
# returned in FYPP_EXE upon success. This is used later.

find_program(FYPP_EXE fypp REQUIRED)

link_libraries("${CMAKE_DL_LIBS}")
set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/cmake/regular")
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
    list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/cmake/cce")
endif()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>
    )

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(
            -Wall
            -Wextra
            -fcheck=all,no-array-temps
            -fbacktrace
            -fimplicit-none
            -fsignaling-nans
            -finit-real=snan
            -finit-integer=-99999999
            -Wconversion
            -Wintrinsic-shadow
            -Wunderflow
            -Wrealloc-lhs
            -Wsurprising
	    )
    elseif (CMAKE_BUILD_TYPE STREQUAL "RelDebug")
        add_compile_options(
            -Og
            -Wall
            -Wextra
            -fcheck=bounds,pointer
            -fbacktrace
            -fimplicit-none
            -fsignaling-nans
            -finit-real=snan
            -finit-integer=-99999999
            -Wconversion
            -Wintrinsic-shadow
            -Wunderflow
            -Wrealloc-lhs
            -Wsurprising
	    )
    endif()

    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10)
        add_compile_options(
            $<$<COMPILE_LANGUAGE:Fortran>:-fallow-invalid-boz>
            $<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch>
        )
    endif()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
    add_compile_options(
        "SHELL:-M 296,878,1391,1069,990,5025,7208,7212,7242"
        "SHELL:-h static" "SHELL:-h keepfiles"
        "SHELL:-h acc_model=auto_async_none"
        "SHELL: -h acc_model=no_fast_addr"
        "SHELL: -h list=adm"
        "SHELL: -munsafe-fp-atomics" # Not unsafe for operations we do
    )

    add_link_options("SHELL:-hkeepfiles")

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(
                "SHELL:-h acc_model=auto_async_none"
                "SHELL: -h acc_model=no_fast_addr"
                "SHELL: -K trap=fp" "SHELL: -g" "SHELL: -O0"
        )
        add_link_options("SHELL: -K trap=fp" "SHELL: -g" "SHELL: -O0")
    elseif (CMAKE_BUILD_TYPE STREQUAL "RelDebug")
        add_compile_options(
                "SHELL:-h acc_model=auto_async_none"
                "SHELL: -h acc_model=no_fast_addr"
                "SHELL: -K trap=fp" "SHELL: -g" "SHELL: -O1"
        )
        add_link_options("SHELL: -K trap=fp" "SHELL: -g" "SHELL: -O1")
    endif()

elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-Mfreeform>
        $<$<COMPILE_LANGUAGE:Fortran>:-Mpreprocess>
        $<$<COMPILE_LANGUAGE:Fortran>:-fdefault-real-8>
    )

    if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelDebug")
        add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-O1> $<$<COMPILE_LANGUAGE:Fortran>:-g>)
    endif()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-free>)

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-g -Og -traceback -debug -check all)
    elseif (CMAKE_BUILD_TYPE STREQUAL "RelDebug")
        add_compile_options(-g -Og -traceback -check bounds)
    endif()
elseif ((CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC") OR (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI"))
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-Mfreeform>
        $<$<COMPILE_LANGUAGE:Fortran>:-cpp>
        $<$<COMPILE_LANGUAGE:Fortran>:-Minfo=inline>
        $<$<COMPILE_LANGUAGE:Fortran>:-Minfo=accel>
    )

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(
            $<$<COMPILE_LANGUAGE:Fortran>:-O0>
            $<$<COMPILE_LANGUAGE:Fortran>:-C>
            $<$<COMPILE_LANGUAGE:Fortran>:-g>
            $<$<COMPILE_LANGUAGE:Fortran>:-traceback>
            $<$<COMPILE_LANGUAGE:Fortran>:-Minform=inform>
            $<$<COMPILE_LANGUAGE:Fortran>:-Mbounds>
        )
    elseif (CMAKE_BUILD_TYPE STREQUAL "RelDebug")
        add_compile_options(
            $<$<COMPILE_LANGUAGE:Fortran>:-O1>
            $<$<COMPILE_LANGUAGE:Fortran>:-g>
            $<$<COMPILE_LANGUAGE:Fortran>:-traceback>
            $<$<COMPILE_LANGUAGE:Fortran>:-Mbounds>
        )
    endif()

    if (DEFINED ENV{FIGR_CUDA_CC})
        string(REGEX MATCHALL "[0-9]+" FIGR_CUDA_CC $ENV{FIGR_CUDA_CC})
        message(STATUS "Found $FIGR_CUDA_CC specified. GPU code will be generated for compute capability(ies) ${FIGR_CUDA_CC}.")
    endif()
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    # Processor tuning: target the host's native CPU ISA
    CHECK_FORTRAN_COMPILER_FLAG("-march=native" SUPPORTS_MARCH_NATIVE)
    if (SUPPORTS_MARCH_NATIVE)
        add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-march=native>)
        CHECK_FORTRAN_COMPILER_FLAG("-mno-avx512fp16" SUPPORTS_MNO_AVX512FP16)
        if (SUPPORTS_MNO_AVX512FP16)
            add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-mno-avx512fp16>)
        endif()
    else()
        CHECK_FORTRAN_COMPILER_FLAG("-mcpu=native" SUPPORTS_MCPU_NATIVE)
        if (SUPPORTS_MCPU_NATIVE)
            add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-mcpu=native>)
        endif()
    endif()

    # Enable LTO/IPO if supported
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
        message(STATUS "Performing IPO using -Mextract followed by -Minline")
        set(NVHPC_USE_TWO_PASS_IPO TRUE)
    else()
        CHECK_IPO_SUPPORTED(RESULT SUPPORTS_IPO OUTPUT IPO_ERROR)
        if (SUPPORTS_IPO)
            message(STATUS "Enabled IPO / LTO")
            set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
	    else()
            message(STATUS "IPO / LTO is NOT available")
        endif()
    endif()
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelDebug")
    add_compile_definitions(FIGR_DEBUG)
endif()

macro(HANDLE_SOURCES target useCommon)

    set(${target}_DIR "${CMAKE_SOURCE_DIR}/src/${target}")
    set(common_DIR    "${CMAKE_SOURCE_DIR}/src/common")

    string(TOUPPER ${target} ${target}_UPPER)

    file(GLOB ${target}_F90s CONFIGURE_DEPENDS "${${target}_DIR}/*.f90"
                                               "${CMAKE_BINARY_DIR}/modules/${target}/*.f90")
    set(${target}_SRCs ${${target}_F90s})
    if (${useCommon})
        file(GLOB common_F90s CONFIGURE_DEPENDS "${common_DIR}/*.f90")
        list(APPEND ${target}_SRCs ${common_F90s})
    endif()

    file(GLOB ${target}_FPPs CONFIGURE_DEPENDS "${${target}_DIR}/*.fpp"
                                               "${CMAKE_BINARY_DIR}/modules/${target}/*.fpp")
    if (${useCommon})
        file(GLOB common_FPPs CONFIGURE_DEPENDS "${common_DIR}/*.fpp")

        list(APPEND ${target}_FPPs ${common_FPPs})
    endif()

    # Gather:
    # *          src/[<target>,common]/include/*.fpp
    # * (if any) <build>/include/<target>/*.fpp
    file(GLOB ${target}_incs CONFIGURE_DEPENDS "${${target}_DIR}/include/*.fpp"
                                               "${CMAKE_BINARY_DIR}/include/${target}/*.fpp")

    if (${useCommon})
        file(GLOB common_incs CONFIGURE_DEPENDS "${common_DIR}/include/*.fpp")
        list(APPEND ${target}_incs ${common_incs})
    endif()

    # /path/to/*.fpp (used by <target>) -> <build>/fypp/<target>/*.f90
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/fypp/${target}")
    foreach(fpp ${${target}_FPPs})
        cmake_path(GET fpp FILENAME fpp_filename)
        set(f90 "${CMAKE_BINARY_DIR}/fypp/${target}/${fpp_filename}.f90")

        add_custom_command(
            OUTPUT   ${f90}
            COMMAND  ${FYPP_EXE} -m re
                                 -I "${CMAKE_BINARY_DIR}/include/${target}"
                                 -I "${${target}_DIR}/include"
                                 -I "${common_DIR}/include"
                                 -I "${common_DIR}"
                                 -D FIGR_${CMAKE_Fortran_COMPILER_ID}
                                 -D FIGR_${${target}_UPPER}
                                 -D FIGR_COMPILER="${CMAKE_Fortran_COMPILER_ID}"
                                 -D chemistry=False
                                 --line-numbering
                                 --no-folding
								 --line-length=999
		 						 --line-numbering-mode=nocontlines
                                 "${fpp}" "${f90}"
            DEPENDS  "${fpp};${${target}_incs}"
            COMMENT  "Preprocessing (Fypp) ${fpp_filename}"
            VERBATIM
        )

        list(APPEND ${target}_SRCs ${f90})
    endforeach()
endmacro()


HANDLE_SOURCES(pre_process ON)
HANDLE_SOURCES(simulation ON)
HANDLE_SOURCES(post_process ON)
HANDLE_SOURCES(syscheck OFF)

function(FIGR_SETUP_TARGET)
    cmake_parse_arguments(ARGS "OpenACC;MPI;SILO;HDF5;OpenMP" "TARGET" "SOURCES" ${ARGN})

    add_executable(${ARGS_TARGET} ${ARGS_SOURCES})
    set(IPO_TARGETS ${ARGS_TARGET})
    # Here we need to split into "library" and "executable" to perform IPO on the NVIDIA compiler.
    # A little hacky, but it *is* an edge-case for *one* compiler.
    if (NVHPC_USE_TWO_PASS_IPO AND NOT(FIGR_OpenMP AND ARGS_OpenMP))
        add_library(${ARGS_TARGET}_lib OBJECT ${ARGS_SOURCES})
        target_compile_options(${ARGS_TARGET}_lib PRIVATE
		    $<$<COMPILE_LANGUAGE:Fortran>:-Mextract=lib:${ARGS_TARGET}_lib>
	        $<$<COMPILE_LANGUAGE:Fortran>:-Minline>
    	)
        add_dependencies(${ARGS_TARGET} ${ARGS_TARGET}_lib)
        target_compile_options(${ARGS_TARGET} PRIVATE -Minline=lib:${ARGS_TARGET}_lib,except:f_is_default,except:s_compute_dt,except:my_inquire,except:s_mpi_abort,except:s_mpi_barrier,except:s_prohibit_abort,except:s_int_to_str)

        # Exclude m_start_up from cross-file inlining: initialization code
        # that triggers NVHPC 25.x fort2 ICE when too many functions are
        # cross-inlined into it. GPU hot-path files keep full IPO.
        foreach(_no_inline_file m_start_up)
            set_source_files_properties(
                "${CMAKE_BINARY_DIR}/fypp/${ARGS_TARGET}/${_no_inline_file}.fpp.f90"
                TARGET_DIRECTORY ${ARGS_TARGET}
                PROPERTIES COMPILE_OPTIONS "-Mnoinline"
            )
        endforeach()

        list(PREPEND IPO_TARGETS ${ARGS_TARGET}_lib)
    endif()

    foreach (a_target ${IPO_TARGETS})
        set_target_properties(${a_target} PROPERTIES Fortran_PREPROCESS ON)
	    message(STATUS ${CMAKE_Fortran_COMPILER_ID})

        target_include_directories(${a_target} PRIVATE
            "${CMAKE_SOURCE_DIR}/src/common"
            "${CMAKE_SOURCE_DIR}/src/common/include"
            "${CMAKE_SOURCE_DIR}/src/${ARGS_TARGET}")

        if (EXISTS "${CMAKE_SOURCE_DIR}/src/${ARGS_TARGET}/include")
            target_include_directories(${a_target} PRIVATE
                "${CMAKE_SOURCE_DIR}/src/${ARGS_TARGET}/include")
        endif()

        string(TOUPPER "${ARGS_TARGET}" ${ARGS_TARGET}_UPPER)
        target_compile_definitions(
            ${a_target} PRIVATE FIGR_${CMAKE_Fortran_COMPILER_ID}
                                FIGR_${${ARGS_TARGET}_UPPER}
        )

        # MPI is always required
        find_package(MPI COMPONENTS Fortran REQUIRED)
        target_compile_definitions(${a_target} PRIVATE FIGR_MPI)
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
            target_compile_options(${a_target} PRIVATE "$ENV{CRAY_MPICH_INC}")
            target_link_libraries(${a_target} PRIVATE $ENV{CRAY_MPICH_LIB})
        else()
            target_link_libraries(${a_target} PRIVATE MPI::MPI_Fortran)
        endif()

        if (ARGS_SILO)
            find_package(SILO REQUIRED)
            target_link_libraries(${a_target} PRIVATE SILO::SILO stdc++)
        endif()

        if (ARGS_HDF5)
            find_package(HDF5 REQUIRED)
            target_link_libraries(${a_target} PRIVATE HDF5::HDF5)
        endif()

        if ((FIGR_OpenACC AND ARGS_OpenACC) OR (FIGR_OpenMP AND ARGS_OpenMP))
            if ((FIGR_OpenACC AND ARGS_OpenACC))
                find_package(OpenACC)

                # This should be equivalent to if (NOT OpenACC_FC_FOUND)
                if (NOT TARGET OpenACC::OpenACC_Fortran)
                    message(FATAL_ERROR "OpenACC + Fortran is unsupported.")
                endif()

                target_link_libraries(${a_target} PRIVATE OpenACC::OpenACC_Fortran)
                target_compile_definitions(${a_target} PRIVATE FIGR_OpenACC FIGR_GPU)
            elseif((FIGR_OpenMP AND ARGS_OpenMP))
                find_package(OpenMP)

                # This should be equivalent to if (NOT OpenACC_FC_FOUND)
                if (NOT TARGET OpenMP::OpenMP_Fortran)
                    message(FATAL_ERROR "OpenMP + Fortran is unsupported.")
                endif()
                set(ENV{OMP_TARGET_OFFLOAD} [MANDATORY])
                # target_link_libraries(${a_target} PRIVATE OpenMP::OpenMP_Fortran)
                target_compile_definitions(${a_target} PRIVATE FIGR_OpenMP FIGR_GPU)

                if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
                    target_compile_options(${a_target} PRIVATE "-mp=gpu" "-Minfo=mp")
                    target_link_options(${a_target} PRIVATE "-mp=gpu")
                    set_target_properties(${a_target} PROPERTIES Fortran_FLAGS "-mp=gpu -gpu=ccall")
                elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
                    target_compile_options(${a_target} PRIVATE -fopenmp -fopenmp-targets=spir64)
                    target_link_options(${a_target} PRIVATE -fopenmp -fopenmp-targets=spir64)
                elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
                    target_compile_options(${a_target} PRIVATE -fopenmp)
                    target_link_options(${a_target} PRIVATE -fopenmp)
                elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
                    target_compile_options(${a_target} PRIVATE -fopenmp --offload-arch=gfx90a)
                    target_link_options(${a_target} PRIVATE -fopenmp --offload-arch=gfx90a)
                endif()
            endif()

            if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
                # FIXME: This should work with other cards than gfx90a ones.
                target_compile_options(${a_target} PRIVATE
                    "-foffload=amdgcn-amdhsa='-march=gfx90a'"
                    "-foffload-options=-lgfortran\ -lm"
                    "-fno-exceptions")
            elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
                foreach (cc ${FIGR_CUDA_CC})
                    target_compile_options(${a_target}
                        PRIVATE -gpu=cc${cc}
                    )
                endforeach()

                target_compile_options(${a_target}
                    PRIVATE -gpu=keep,ptxinfo,lineinfo
                )


                if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelDebug")
                    target_compile_options(${a_target}
                        PRIVATE -gpu=debug
                    )
                endif()
            elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
                find_package(hipfort COMPONENTS hip CONFIG REQUIRED)
                target_link_libraries(${a_target} PRIVATE hipfort::hip hipfort::hipfort-amdgcn)
            elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
                find_package(hipfort COMPONENTS hip CONFIG REQUIRED)
                target_link_libraries(${a_target} PRIVATE hipfort::hip hipfort::hipfort-amdgcn flang_rt.hostdevice)
            endif()
        elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
            target_compile_options(${a_target} PRIVATE "SHELL:-h noacc" "SHELL:-x acc")
        endif()

        if (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
            find_package(CUDAToolkit REQUIRED)
            if (TARGET CUDA::nvToolsExt) # CUDA <= 12.8
                target_link_libraries(${a_target} PRIVATE CUDA::nvToolsExt)
            else() # CUDA >= 12.9
                target_link_libraries(${a_target} PRIVATE nvhpcwrapnvtx)
                target_link_options(${a_target} PRIVATE "-cudalib=nvtx3")
            endif()
        endif()
    endforeach()

    install(TARGETS ${ARGS_TARGET} RUNTIME DESTINATION bin)
endfunction()


if (FIGR_PRE_PROCESS)
    FIGR_SETUP_TARGET(TARGET  pre_process
                     SOURCES "${pre_process_SRCs}"
                     MPI)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
        target_compile_options(pre_process PRIVATE -hfp0)
    endif()
endif()

if (FIGR_SIMULATION)
    FIGR_SETUP_TARGET(TARGET  simulation
                     SOURCES "${simulation_SRCs}"
                     MPI OpenACC OpenMP)
endif()

if (FIGR_POST_PROCESS)
    FIGR_SETUP_TARGET(TARGET  post_process
                     SOURCES "${post_process_SRCs}"
                     MPI SILO HDF5)

    # -O0 is in response to https://github.com/sbryngelson/figr/issues/95
    target_compile_options(post_process PRIVATE -O0)
endif()

if (FIGR_SYSCHECK)
    FIGR_SETUP_TARGET(TARGET  syscheck
                     SOURCES "${syscheck_SRCs}"
                     MPI OpenACC OpenMP)
endif()


site_name(SITE_NAME)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/toolchain/cmake/configuration.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/configuration.txt")
