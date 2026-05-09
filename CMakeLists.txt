cmake_minimum_required(VERSION 3.15)
project(liscio VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(BUILD_TESTS "Build test programs" ON)
option(BUILD_BENCHMARK "Build benchmark program" ON)
option(BUILD_SHARED "Build shared library in addition to static" ON)

# ============================================================================
# Library sources
# ============================================================================
set(LISCIO_SOURCES
    src/liscio.c
    src/liscio_arc9d.c
    src/liscio_helix9.c
    src/liscio_bezier9.c
    src/liscio_rbezier.c
    src/liscio_bspline9.c
)

# ============================================================================
# Static library
# ============================================================================
add_library(liscio STATIC ${LISCIO_SOURCES})
target_include_directories(liscio
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(liscio PRIVATE m)
target_compile_options(liscio PRIVATE -Wall -Wextra -Wno-unused-parameter -O2)
set_target_properties(liscio PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    C_VISIBILITY_PRESET hidden
)

# ============================================================================
# Shared library (optional)
# ============================================================================
if(BUILD_SHARED)
    add_library(liscio_shared SHARED ${LISCIO_SOURCES})
    set_target_properties(liscio_shared PROPERTIES
        OUTPUT_NAME liscio
        VERSION     ${PROJECT_VERSION}
        SOVERSION   ${PROJECT_VERSION_MAJOR}
    )
    target_include_directories(liscio_shared
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src
    )
    target_link_libraries(liscio_shared PRIVATE m)
    target_compile_options(liscio_shared PRIVATE -Wall -Wextra -Wno-unused-parameter -O2)
endif()

# ============================================================================
# Tests
# ============================================================================
if(BUILD_TESTS)
    enable_testing()

    # Shared NGC parser (handles G17/G18/G19, G20/G21, G90/G91, O-word while).
    set(NGC_PARSER_SOURCES test/ngc_parser.c)

    foreach(T IN ITEMS test_arc_fit test_helix_fit test_bezier_fit test_rbezier test_bspline test_fuzz test_arc_merge)
        add_executable(${T} test/${T}.c)
        target_link_libraries(${T} PRIVATE liscio m)
        target_include_directories(${T} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
        target_compile_definitions(${T} PRIVATE _POSIX_C_SOURCE=199309L)
    endforeach()
    foreach(T IN ITEMS test_ngc_replay test_dump_csv test_metrics test_compare test_verify)
        add_executable(${T} test/${T}.c ${NGC_PARSER_SOURCES})
        target_link_libraries(${T} PRIVATE liscio m)
        target_include_directories(${T} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/test)
        target_compile_definitions(${T} PRIVATE _POSIX_C_SOURCE=199309L)
    endforeach()

    # ------------------------------------------------------------------
    # Cross-library benchmark (test_xcompare): vendors Schneider FitCurves
    # 3D port; conditionally links SISL (Debian libsisl-dev) and uses a
    # Python venv for geomdl (Piegl & Tiller §9.4 LSQ B-spline reference).
    # ------------------------------------------------------------------
    set(XCOMPARE_SOURCES
        test/test_xcompare.c
        test/third_party/schneider3d.c
        ${NGC_PARSER_SOURCES})
    add_executable(test_xcompare ${XCOMPARE_SOURCES})
    target_link_libraries(test_xcompare PRIVATE liscio m)
    target_include_directories(test_xcompare PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/test)
    target_compile_definitions(test_xcompare PRIVATE _POSIX_C_SOURCE=199309L)

    # SISL (optional — Debian package libsisl-dev provides /usr/include/sisl.h
    # and -lsisl).  Detect via find_path / find_library, fall back gracefully.
    find_path(SISL_INCLUDE_DIR sisl.h)
    find_library(SISL_LIBRARY  NAMES sisl)
    if(SISL_INCLUDE_DIR AND SISL_LIBRARY)
        target_compile_definitions(test_xcompare PRIVATE HAVE_SISL=1)
        target_include_directories(test_xcompare PRIVATE ${SISL_INCLUDE_DIR})
        target_link_libraries(test_xcompare PRIVATE ${SISL_LIBRARY})
        message(STATUS "test_xcompare: SISL enabled (${SISL_LIBRARY})")
    else()
        message(STATUS "test_xcompare: SISL not found, sisl strategy disabled")
    endif()

    add_test(NAME arc_fit     COMMAND test_arc_fit)
    add_test(NAME helix_fit   COMMAND test_helix_fit)
    add_test(NAME bezier_fit  COMMAND test_bezier_fit)
    add_test(NAME rbezier     COMMAND test_rbezier)
    add_test(NAME bspline     COMMAND test_bspline)
    add_test(NAME fuzz        COMMAND test_fuzz 1 30 200)
    add_test(NAME arc_merge   COMMAND test_arc_merge)

    # Geometric regression gate: no fit may exceed tol_xyz against any raw
    # waypoint, AND no fit polyline may drift to > 1.5× its absorbed-chord
    # length.  Catches both pointwise miss and inter-waypoint drift.
    file(GLOB LISCIO_NGC_REGRESSION CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/test/data/*.ngc"
        "${CMAKE_CURRENT_SOURCE_DIR}/test/data/*.NGC")
    if(LISCIO_NGC_REGRESSION)
        add_test(NAME verify_geometry
            COMMAND test_verify --tol_xyz=0.05 --arc_to_helix=16 ${LISCIO_NGC_REGRESSION})
    endif()
endif()

if(BUILD_BENCHMARK)
    add_executable(test_benchmark test/test_benchmark.c)
    target_link_libraries(test_benchmark PRIVATE liscio m rt)
    target_include_directories(test_benchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
endif()

# ============================================================================
# Tools — user-facing utilities (not regression tests).
# ============================================================================
option(BUILD_TOOLS "Build user-facing tools (liscio_compress)" ON)
if(BUILD_TOOLS)
    add_executable(liscio_compress
        tool/liscio_compress.c
        test/ngc_parser.c)
    target_link_libraries(liscio_compress PRIVATE liscio m)
    target_include_directories(liscio_compress PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/test)
    target_compile_definitions(liscio_compress PRIVATE _POSIX_C_SOURCE=199309L)
endif()

# ============================================================================
# Helper targets: dump + plot bundled CAM NGCs
#   make dump_ngcs   → CSV files alongside each NGC in ngc_out/
#   make plot_ngcs   → PNG via matplotlib (needs python3 + matplotlib)
# ============================================================================
file(GLOB LISCIO_NGC_DATA CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/test/data/*.ngc"
    "${CMAKE_CURRENT_SOURCE_DIR}/test/data/*.NGC")
if(LISCIO_NGC_DATA)
    add_custom_command(
        OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/ngc_out/.stamp"
        COMMAND ${CMAKE_COMMAND} -E make_directory
            "${CMAKE_CURRENT_BINARY_DIR}/ngc_out"
        COMMAND ${CMAKE_COMMAND} -E copy ${LISCIO_NGC_DATA}
            "${CMAKE_CURRENT_BINARY_DIR}/ngc_out/"
        COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/ngc_out && \
            shopt -s nullglob 2>/dev/null; \
            ${CMAKE_CURRENT_BINARY_DIR}/test_dump_csv \
            --tol_xyz=0.05 --samples=30 --arc_to_helix=16 *.ngc *.NGC"
        COMMAND ${CMAKE_COMMAND} -E touch
            "${CMAKE_CURRENT_BINARY_DIR}/ngc_out/.stamp"
        DEPENDS test_dump_csv ${LISCIO_NGC_DATA}
        COMMENT "liscio: dump bundled NGCs to ngc_out/"
        VERBATIM
    )
    add_custom_target(dump_ngcs
        DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/ngc_out/.stamp")

    add_custom_target(plot_ngcs
        COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/ngc_out && \
            shopt -s nullglob 2>/dev/null; \
            python3 ${CMAKE_CURRENT_SOURCE_DIR}/test/plot_fit.py *.ngc *.NGC"
        DEPENDS dump_ngcs
        COMMENT "liscio: plot dumped NGCs (needs matplotlib)"
        VERBATIM
    )

    # Auto-regenerate quantitative metrics tables when NGC corpus changes.
    # Default — tol_xyz=0.05 matches G64 P0.05 (industry mainstream
     # metal cutting), with helix auto-detection (default cfg).
    add_custom_command(
        OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/metrics.csv"
        COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_metrics
            --tol_xyz=0.05
            --csv=${CMAKE_CURRENT_BINARY_DIR}/metrics.csv
            ${LISCIO_NGC_DATA}
        DEPENDS test_metrics ${LISCIO_NGC_DATA}
        COMMENT "liscio: regenerate metrics.csv (tol_xyz=0.05, default cfg)"
        VERBATIM
    )
    # Tight precision table — tol_xyz=0.01 (10 µm precision machining).
    add_custom_command(
        OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/metrics_g1.csv"
        COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_metrics
            --tol_xyz=0.01
            --arc_subseg=8
            --csv=${CMAKE_CURRENT_BINARY_DIR}/metrics_g1.csv
            ${LISCIO_NGC_DATA}
        DEPENDS test_metrics ${LISCIO_NGC_DATA}
        COMMENT "liscio: regenerate metrics_g1.csv (tol_xyz=0.01 arc_subseg=8)"
        VERBATIM
    )
    add_custom_target(metrics_csv ALL
        DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/metrics.csv"
                "${CMAKE_CURRENT_BINARY_DIR}/metrics_g1.csv")

    # ----------------------------------------------------------------
    # `liscio_results` target — single command that syncs everything
    # the regression archive needs:
    #   1. dump CSVs   (test_dump_csv per ngc)
    #   2. plot PNGs   (matplotlib)
    #   3. metrics CSVs (test_metrics, both default and arc_subseg=8)
    #   4. copy PNGs + CSVs into results/{images,data}/
    #
    # Triggered manually (`cmake --build build --target liscio_results`).
    # Not in ALL because it needs matplotlib at runtime; metrics_csv
    # alone is in ALL for build-time freshness.
    # ----------------------------------------------------------------
    add_custom_target(liscio_results
        COMMAND bash
            "${CMAKE_CURRENT_SOURCE_DIR}/cmake/sync_results.sh"
            "${CMAKE_CURRENT_BINARY_DIR}"
            "${CMAKE_CURRENT_SOURCE_DIR}"
        DEPENDS plot_ngcs metrics_csv
        COMMENT "liscio: sync ngc_out PNGs + metrics CSVs into results/"
    )
endif()

# ============================================================================
# Install
# ============================================================================
install(TARGETS liscio
    EXPORT liscio-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT development
)
if(BUILD_SHARED)
    install(TARGETS liscio_shared
        EXPORT liscio-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT runtime
    )
endif()

install(DIRECTORY include/liscio
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT development
    FILES_MATCHING PATTERN "*.h"
)

# CMake package config
set(LISCIO_INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

install(EXPORT liscio-targets
    DESTINATION ${LISCIO_INSTALL_CONFIGDIR}
    NAMESPACE liscio::
    COMPONENT development
)

configure_package_config_file(
    cmake/liscio-config.cmake.in
    liscio-config.cmake
    INSTALL_DESTINATION ${LISCIO_INSTALL_CONFIGDIR}
)
write_basic_package_version_file(
    liscio-config-version.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/liscio-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/liscio-config-version.cmake"
    DESTINATION ${LISCIO_INSTALL_CONFIGDIR}
    COMPONENT development
)

message(STATUS "liscio v${PROJECT_VERSION} — pure C geometry preprocessor")
