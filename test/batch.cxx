#include "catch2/catch.hpp"
#include <integratorxx/quadrature.hpp>
#include <integratorxx/batch.hpp>
#include <cmath>
#include <complex>

TEST_CASE( "Spherical Micro Batches", "[sph-batch]" ) {

  using namespace IntegratorXX;
  EulerMaclaurin<double> rquad(150);
  Lebedev<double> aquad(770);

  auto sph_micro_batch = SphericalMicroBatch( rquad, aquad );

}
