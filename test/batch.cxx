#include "catch2/catch.hpp"
#include <integratorxx/deprecated/quadrature.hpp>
#include <integratorxx/deprecated/batch.hpp>
#include <cmath>
#include <complex>

TEST_CASE( "Spherical Micro Batches", "[sph-batch]" ) {

  using namespace IntegratorXX;
  IntegratorXX::deprecated::EulerMaclaurin<double> rquad(150);
  IntegratorXX::deprecated::Lebedev<double> aquad(770);

  auto sph_micro_batch = IntegratorXX::deprecated::SphericalMicroBatch( rquad, aquad );

}
