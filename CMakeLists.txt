# CMake configuration for ECOS

cmake_minimum_required(VERSION 3.5)
project(ecos C)

# Options
option(USE_LONG "Whether to use Long or Int for index type." TRUE)
option(ECOS_ENABLE_TESTING "Whether to compile tests." FALSE)

# AMD
add_library(amd SHARED)
target_include_directories(amd SYSTEM PUBLIC
    external/amd/include
    external/SuiteSparse_config
)
target_sources(amd PRIVATE
    external/amd/src/amd_1.c
    external/amd/src/amd_2.c
    external/amd/src/amd_aat.c
    external/amd/src/amd_control.c
    external/amd/src/amd_defaults.c
    external/amd/src/amd_dump.c
    external/amd/src/amd_global.c
    external/amd/src/amd_info.c
    external/amd/src/amd_order.c
    external/amd/src/amd_post_tree.c
    external/amd/src/amd_postorder.c
    external/amd/src/amd_preprocess.c
    external/amd/src/amd_valid.c
)

# LDL
add_library(ldl SHARED)
target_include_directories(ldl SYSTEM PUBLIC
    include
    external/ldl/include
    external/SuiteSparse_config
)
target_sources(ldl PRIVATE
    external/ldl/src/ldl.c
)

# ECOS
add_library(ecos SHARED)
target_include_directories(ecos SYSTEM PUBLIC
    include
    external/SuiteSparse_config
)
target_sources(ecos PRIVATE
    src/cone.c
    src/ctrlc.c
    src/ecos.c
    src/equil.c
    src/expcone.c
    src/kkt.c
    src/preproc.c
    src/spla.c
    src/splamm.c
    src/timer.c
    src/wright_omega.c

    ecos_bb/ecos_bb.c
    ecos_bb/ecos_bb_preproc.c
)

# Set compiler flags
set(DEBUG_OPTIONS -Wall -Wextra)
set(RELEASE_OPTIONS -O2)
target_compile_options(ecos PUBLIC "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")
target_compile_options(ecos PUBLIC "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")

# Set compiler definitions
target_compile_definitions(ecos PUBLIC fPIC CTRLC=1)
if(USE_LONG)
    target_compile_definitions(ldl PUBLIC LDL_LONG DLONG)
    target_compile_definitions(amd PUBLIC LDL_LONG DLONG)
endif()

target_link_libraries(ecos amd ldl m)

# Tests
if(ECOS_ENABLE_TESTING)
    add_library(ecos_test)

    target_include_directories(ecos_test SYSTEM PUBLIC
        test
        test/generated
    )

    target_sources(ecos_test PRIVATE
        test/generated/inv_pos/inv_pos.c
        test/generated/norm/norm.c
        test/generated/quad_over_lin/quad_over_lin.c
        test/generated/sq_norm/sq_norm.c
        test/generated/sum_sq/sum_sq.c
        test/generated/qcml_utils.c
    )

    target_link_libraries(ecos_test ecos)

    add_executable(runecos src/runecos.c)
    target_link_libraries(runecos ecos_test)

    add_executable(runecos_exp src/runecos_exp.c)
    target_link_libraries(runecos_exp ecos_test)

    add_executable(ecostester test/ecostester.c)
    target_link_libraries(ecostester ecos_test)

    add_executable(bb_test test/bb_test.c)
    target_link_libraries(bb_test ecos_test)
endif()