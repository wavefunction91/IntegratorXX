cmake_minimum_required( VERSION 3.17 FATAL_ERROR )
project( integratorxx_cmake_discovery LANGUAGES C CXX )

find_package( IntegratorXX REQUIRED )
add_executable( integratorxx_link_tester integratorxx_link_tester.cxx )
target_link_libraries( integratorxx_link_tester PUBLIC IntegratorXX::IntegratorXX )
