# CMake configuration for ECOS

cmake_minimum_required(VERSION 3.5)

project(ecos
        VERSION 2.0.8
        LANGUAGES C)

# Options
option(USE_LONG "Whether to use Long or Int for index type." ON)

# Set compiler flags
set(DEBUG_OPTIONS "")
set(RELEASE_OPTIONS "")

if(MSVC)
    list(APPEND DEBUG_OPTIONS "/Od" "/Wall")
    list(APPEND RELEASE_OPTIONS "/O2")
else()
    list(APPEND DEBUG_OPTIONS "-O0" "-Wall" "-Wextra")
    list(APPEND RELEASE_OPTIONS "-O3")
endif()

# Add the ECOS headers
set(ecos_headers
    "${CMAKE_CURRENT_SOURCE_DIR}/external/SuiteSparse_config/SuiteSparse_config.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/cone.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/ctrlc.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/data.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/ecos.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/ecos_bb.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/equil.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/expcone.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/glblopts.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/kkt.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/spla.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/splamm.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/timer.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/include/wright_omega.h"
)

set(ecos_sources
    # AMD
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_1.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_2.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_aat.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_control.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_defaults.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_dump.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_global.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_info.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_order.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_post_tree.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_postorder.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_preprocess.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/amd/src/amd_valid.c"

    # LDL
    "${CMAKE_CURRENT_SOURCE_DIR}/external/ldl/src/ldl.c"

    # ECOS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/cone.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ctrlc.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ecos.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/equil.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/expcone.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/kkt.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/preproc.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/spla.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/splamm.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/timer.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/wright_omega.c"

    # ECOS BB
    "${CMAKE_CURRENT_SOURCE_DIR}/ecos_bb/ecos_bb.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/ecos_bb/ecos_bb_preproc.c"
)

# ECOS library
add_library(ecos SHARED ${ecos_headers} ${ecos_sources})

# Set compiler options and definitions
target_compile_options(ecos PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")
target_compile_options(ecos PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")
target_compile_definitions(ecos PUBLIC CTRLC=1)
if(USE_LONG)
    target_compile_definitions(ecos PUBLIC LDL_LONG DLONG)
endif()

# Link math library
if(NOT MSVC)
    target_link_libraries(ecos PRIVATE m)
endif()

target_include_directories(ecos
    PUBLIC 
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/amd/include>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/ldl/include>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/SuiteSparse_config>"
    "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}/ecos>"
)

# Installation
include(GNUInstallDirs)

install(TARGETS ecos
        EXPORT  ${PROJECT_NAME}
        ARCHIVE       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        LIBRARY       DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        RUNTIME       DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

install(FILES ${ecos_headers} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ecos")

# Create CMake packages for the build directory
include(CMakePackageConfigHelpers)

export(EXPORT ${PROJECT_NAME}
  FILE "${CMAKE_CURRENT_BINARY_DIR}/ecos-targets.cmake"
  NAMESPACE ecos::
)

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/ecos-config.cmake "include(\"\${CMAKE_CURRENT_LIST_DIR}/ecos-targets.cmake\")\n")

# Create CMake packages for the install directory
set(ConfigPackageLocation ${CMAKE_INSTALL_LIBDIR}/cmake/ecos)

install(EXPORT ${PROJECT_NAME}
        FILE ecos-targets.cmake
        NAMESPACE ecos::
        DESTINATION ${ConfigPackageLocation}
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ecos-config.cmake
        DESTINATION ${ConfigPackageLocation}
)
