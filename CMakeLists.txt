cmake_minimum_required( VERSION 3.11 ) # Require CMake 3.11+

# Set up project definition + version information
project( IntegratorXX VERSION 1.2.0 LANGUAGES CXX )


add_library( integratorxx INTERFACE )
add_library( IntegratorXX::integratorxx ALIAS integratorxx )

option( INTEGRATORXX_ENABLE_TESTS "Enable Tets" OFF )

if( INTEGRATORXX_ENABLE_TESTS )
  find_package( Catch2 QUIET )
  if( Catch2-FOUND )
    target_link_libraries( integratorxx INTERFACE Catch2::Catch2)
  else()
  
    set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
    include( DownloadProject )
  
    download_project(
      PROJ                Catch2
      GIT_REPOSITORY      https://github.com/catchorg/Catch2.git
      GIT_TAG             master
      UPDATE_DISCONNECTED 1
    )
  
    add_subdirectory(${Catch2_SOURCE_DIR} ${Catch2_BINARY_DIR})
  
  endif()

  add_subdirectory( test )
endif()


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

