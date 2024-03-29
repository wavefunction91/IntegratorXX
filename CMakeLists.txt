cmake_minimum_required( VERSION 3.17 ) # Require CMake 3.17+

# Set up project definition + version information
project( IntegratorXX VERSION 1.2.0 LANGUAGES CXX )

option(INTEGRATORXX_HEADER_ONLY "Force header-only build" OFF)

if(INTEGRATORXX_HEADER_ONLY)
  add_library( integratorxx INTERFACE )
  set(INTEGRATORXX_TARGET_TYPE INTERFACE)
else()
  add_subdirectory(src)
  set(INTEGRATORXX_TARGET_TYPE PUBLIC)
endif()

# Target features
target_compile_features( integratorxx ${INTEGRATORXX_TARGET_TYPE} cxx_std_17 )
target_include_directories( integratorxx
  ${INTEGRATORXX_TARGET_TYPE} 
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-Wno-missing-braces" INTEGRATORXX_HAS_NO_MISSING_BRACES )
if( INTEGRATORXX_HAS_NO_MISSING_BRACES )
	target_compile_options( integratorxx INTERFACE $<$<COMPILE_LANGUAGE:CXX>: -Wno-missing-braces> )
endif()




# INSTALL rules
add_library( IntegratorXX::IntegratorXX ALIAS integratorxx )
include( GNUInstallDirs )
set( INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/IntegratorXX )

# Targets
install(TARGETS integratorxx
    EXPORT integratorxx-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

set_target_properties( integratorxx PROPERTIES EXPORT_NAME IntegratorXX )

# Install Headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )

# Export target to script
install(EXPORT integratorxx-targets
  FILE
    IntegratorXXTargets.cmake
  NAMESPACE
    IntegratorXX::
  DESTINATION
    ${INSTALL_CONFIGDIR}
)

# Export build-tree targets also (e.g. to be usable by FetchContent)
export(EXPORT integratorxx-targets
    NAMESPACE IntegratorXX::
    FILE "${PROJECT_BINARY_DIR}/IntegratorXXTargets.cmake")

# Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/IntegratorXXConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

# Setup IntegratorXXConfig.cmake
configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/IntegratorXXConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/IntegratorXXConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

#Install the config, configversion and custom find modules
install(DIRECTORY 
  ${PROJECT_SOURCE_DIR}/cmake/ 
  DESTINATION ${INSTALL_CONFIGDIR} 
  FILES_MATCHING PATTERN "*.cmake"
)


install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/IntegratorXXConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/IntegratorXXConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

# Testing
if( NOT DEFINED INTEGRATORXX_ENABLE_TESTS )
  set( INTEGRATORXX_ENABLE_TESTS ON )
endif()


if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
endif()

if( CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND INTEGRATORXX_ENABLE_TESTS AND BUILD_TESTING )
  add_subdirectory( test )
endif()



