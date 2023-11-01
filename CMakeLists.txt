cmake_minimum_required( VERSION 3.17 ) # Require CMake 3.17+

# Set up project definition + version information
project( IntegratorXX VERSION 1.2.0 LANGUAGES CXX )

#add_library( integratorxx INTERFACE )
add_subdirectory(src)

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



