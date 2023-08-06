find_package( Catch2 3.0.1 CONFIG)
if( NOT Catch2_FOUND )
  message( STATUS "Could NOT Find Catch2 (Building v3.3.2 Locally)" )
  include( FetchContent )
  FetchContent_Declare(
    catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.3.2
  )
  if( NOT catch2_POPULATED )
    FetchContent_Populate( catch2 )
    set(CATCH_BUILD_TESTING OFF CACHE BOOL "Build SelfTest project" FORCE)
    set(CATCH_INSTALL_DOCS OFF CACHE BOOL "Install documentation alongside library" FORCE)
    set(CATCH_INSTALL_HELPERS OFF CACHE BOOL "Install contrib alongside library" FORCE)
    add_subdirectory( ${catch2_SOURCE_DIR} ${catch2_BINARY_DIR} )
  endif()
else()
  message( STATUS "Found Catch2 VERSION ${Catch2_VERSION} DIR ${Catch2_DIR}" )
endif()

add_library(integratorxx_common_ut quad_matcher.cxx)
target_link_libraries( integratorxx_common_ut PUBLIC Catch2::Catch2WithMain integratorxx )

add_executable( quadrature_manipulation quadrature_manipulation.cxx )
target_link_libraries( quadrature_manipulation PUBLIC integratorxx_common_ut )

add_executable( 1d_quadratures 1d_quadratures.cxx )
target_link_libraries( 1d_quadratures PUBLIC integratorxx_common_ut )

add_executable( gausslegendre gausslegendre.cxx )
target_link_libraries( gausslegendre PUBLIC integratorxx_common_ut )

add_executable( gausslobatto gausslobatto.cxx )
target_link_libraries( gausslobatto PUBLIC integratorxx_common_ut )
add_executable( muraknowles muraknowles.cxx )
target_link_libraries( muraknowles PUBLIC integratorxx_common_ut )

add_executable( gausschebyshev1 gausschebyshev1.cxx )
target_link_libraries( gausschebyshev1 PUBLIC Catch2::Catch2WithMain integratorxx )

add_executable( gausschebyshev2 gausschebyshev2.cxx )
target_link_libraries( gausschebyshev2 PUBLIC Catch2::Catch2WithMain integratorxx )

add_executable( composite_quadratures composite_quadratures.cxx )
target_link_libraries( composite_quadratures PUBLIC integratorxx_common_ut )

add_test( NAME QUADRATURE_MANIP         COMMAND quadrature_manipulation )
add_test( NAME QUADRATURES_1D           COMMAND 1d_quadratures          )
add_test( NAME QUADRATURES_COMPOSITE    COMMAND composite_quadratures   )
add_test( NAME QUADRATURES_LEGENDRE     COMMAND gausslegendre           )
add_test( NAME QUADRATURES_LOBATTO      COMMAND gausslobatto            )
add_test( NAME QUADRATURES_CHEBYSHEV1   COMMAND gausschebyshev1         )
add_test( NAME QUADRATURES_CHEBYSHEV2   COMMAND gausschebyshev2         )
