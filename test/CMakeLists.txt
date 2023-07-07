find_package( Catch2 CONFIG QUIET )
if( NOT Catch2_FOUND )
  include( FetchContent )
  FetchContent_Declare(
    catch2_download
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.3.2
  )
  if( NOT catch2_download_POPULATED )
    FetchContent_Populate( catch2_download )
    set(CATCH_BUILD_TESTING OFF CACHE BOOL "Build SelfTest project" FORCE)
    set(CATCH_INSTALL_DOCS OFF CACHE BOOL "Install documentation alongside library" FORCE)
    set(CATCH_INSTALL_HELPERS OFF CACHE BOOL "Install contrib alongside library" FORCE)
    add_subdirectory( ${catch2_download_SOURCE_DIR} ${catch2_download_BINARY_DIR} )
  endif()
endif()


#add_executable( my_test test.cxx )
#target_link_libraries( my_test       PUBLIC integratorxx )


#add_library( catch2_main STATIC ut_main.cxx )
#target_link_libraries( catch2_main PUBLIC Catch2::Catch2 )


add_executable( quadrature_manipulation quadrature_manipulation.cxx )
target_link_libraries( quadrature_manipulation PUBLIC Catch2Main Catch2 integratorxx )

add_executable( 1d_quadratures 1d_quadratures.cxx )
target_link_libraries( 1d_quadratures PUBLIC Catch2Main Catch2 integratorxx )

add_executable( composite_quadratures composite_quadratures.cxx )
target_link_libraries( composite_quadratures PUBLIC Catch2Main Catch2 integratorxx )

add_test( NAME QUADRATURE_MANIP      COMMAND quadrature_manipulation )
add_test( NAME QUADRATURES_1D        COMMAND 1d_quadratures          )
add_test( NAME QUADRATURES_COMPOSITE COMMAND composite_quadratures   )
