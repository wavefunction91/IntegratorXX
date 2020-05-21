#include "catch2/catch.hpp"

#include <integratorxx/quadratures/muraknowles.hpp>
#include <integratorxx/quadratures/lebedev_laikov.hpp>
#include <integratorxx/composite_quadratures/spherical_quadrature.hpp>
#include <integratorxx/batch/spherical_micro_batcher.hpp>

//#include <integratorxx/deprecated/batch.hpp>
//#include <integratorxx/deprecated/quadrature.hpp>

TEST_CASE( "Spherical Quadratures", "[sph-quad]" ) {

  size_t nrad = 150;
  size_t nang = 770;

  SECTION("Correctness") {
    IntegratorXX::MuraKnowles<double,double> r(nrad);
    IntegratorXX::LebedevLaikov<double> q(nang);

    IntegratorXX::SphericalQuadrature s( r, q );
    size_t npts = s.npts();

    REQUIRE( npts == nrad * nang );

    double res = 0;
    for( auto i = 0; i < npts; ++i ) {
      const auto p = s.points()[i];
      const auto rsq = p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
      res += s.weights()[i] * std::exp(-rsq);
    }

    CHECK( res == Approx( M_PI * std::sqrt(M_PI) ) );
  }

#if 0
  SECTION("Consistency") {

    std::array<double,3> cen = {0., 0., 1.};
    double scale_factor = 0.5;
    IntegratorXX::MuraKnowles<double,double> r(nrad, scale_factor);
    IntegratorXX::LebedevLaikov<double> q(nang);
    
    IntegratorXX::deprecated::Knowles<double>       r_old(nrad);
    IntegratorXX::deprecated::Lebedev<double> q_old(nang);

#if 0
    IntegratorXX::SphericalQuadrature s( r, q, cen );
#else
    IntegratorXX::SphericalQuadrature s( r, q );
    s.recenter( cen );
#endif

    auto s_old = IntegratorXX::deprecated::SphericalBatch( r_old, q_old,
      cen, scale_factor, nrad, nang );

    auto [points_old, weights_old] = *s_old.begin();
    

    REQUIRE( points_old.size() == s.npts() );

    for( auto i = 0; i < s.npts(); ++i ) {
      auto p_old = points_old[i];
      auto w_old = weights_old[i];

      auto p_new = s.points()[i];
      auto w_new = s.weights()[i];

      CHECK( w_old == Approx( w_new ) );
      CHECK( p_old[0] == Approx( p_new[0] ) );
      CHECK( p_old[1] == Approx( p_new[1] ) );
      CHECK( p_old[2] == Approx( p_new[2] ) );
    }
  }
#endif

  SECTION("Batching") {

    size_t max_batch_sz = 512;
    IntegratorXX::MuraKnowles<double,double> r(nrad);
    IntegratorXX::LebedevLaikov<double> q(nang);
    IntegratorXX::SphericalQuadrature s( r, q );

    size_t npts = s.npts();
    IntegratorXX::SphericalMicroBatcher batcher( max_batch_sz, s );

    size_t npts_c = 0;
    for( auto&& [box_lo, box_up, points_b, weights_b] : batcher ) {

      auto npts_b = points_b.size();
      CHECK( npts_b != 0 );

      npts_c += npts_b;
      //for( auto& p : points_b ) {
      //  auto pc = std::count( s.points().begin(), s.points().end(), p );
      //  CHECK( pc == 1 );
      //}

    }

    CHECK( npts_c == npts );

    auto batcher_clone = batcher.clone();
    CHECK( &batcher.quadrature() != &batcher_clone.quadrature() );
    CHECK( batcher.npts() == batcher_clone.npts() );
    CHECK( batcher.nbatches() == batcher_clone.nbatches() );
  
    CHECK( batcher.points() == batcher_clone.points() );
    CHECK( batcher.weights() == batcher_clone.weights() );
  }
};
