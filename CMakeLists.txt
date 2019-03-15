cmake_minimum_required( VERSION 3.11 ) # Require CMake 3.11+

# Set up project definition + version information
project( IntegratorXX VERSION 1.2.0 LANGUAGES CXX )


add_library( integratorxx INTERFACE )
add_library( IntegratorXX::integratorxx ALIAS integratorxx )


# Ranges-V3
find_package( range-v3 QUIET )

if( range-v3_FOUND )

  target_link_libraries( integratorxx INTERFACE range-v3 )

else()

  
  set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
  include( DownloadProject )

  download_project(
    PROJ                range-v3
    GIT_REPOSITORY      https://github.com/ericniebler/range-v3
    GIT_TAG             master
    UPDATE_DISCONNECTED 1
  )

  set(RANGE_V3_TESTS    OFF) 
  set(RANGE_V3_EXAMPLES OFF)
  set(RANGES_CXX_STD    17)

  add_subdirectory(${range-v3_SOURCE_DIR} ${range-v3_BINARY_DIR})
  target_link_libraries( integratorxx INTERFACE range-v3 )
endif()



add_subdirectory( test )

# Target features
target_compile_features( integratorxx INTERFACE cxx_std_17 )
target_include_directories( integratorxx
  INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

# INSTALL rules
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

