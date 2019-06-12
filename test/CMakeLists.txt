#add_executable( my_test test.cxx )
#target_link_libraries( my_test       PUBLIC integratorxx )


add_library( catch2_main STATIC ut_main.cxx )
target_compile_definitions( catch2_main PRIVATE "CATCH_CONFIG_MAIN" )
target_link_libraries( catch2_main PUBLIC Catch2::Catch2 )


add_executable( 1d_quadratures 1d_quadratures.cxx )
target_link_libraries( 1d_quadratures PUBLIC catch2_main integratorxx )

add_executable( batch batch.cxx )
target_link_libraries( batch PUBLIC catch2_main integratorxx )
