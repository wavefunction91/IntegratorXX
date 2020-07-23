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


  SECTION("Centering") {

    std::array<double,3> cen = {0., 1., 0.};
    std::array<double,3> org = {0., 0., 0.};

    IntegratorXX::LebedevLaikov<double>      ang( nang );
    IntegratorXX::MuraKnowles<double,double> rad( nrad );
    
    IntegratorXX::SphericalQuadrature s_origin( rad, ang );
    IntegratorXX::SphericalQuadrature s_cen(rad, ang, cen);
    IntegratorXX::SphericalQuadrature so_cpy( s_origin );
    IntegratorXX::SphericalQuadrature sc_cpy( s_cen    );

    // Center on construction
    CHECK( s_origin.center() == org );
    CHECK( s_cen.center()    == cen );

    auto npts = s_origin.npts();
    CHECK( s_cen.npts() == npts );

    for(size_t i = 0; i < npts; ++i) {
      auto po = s_origin.points()[i];
      auto pc = s_cen.points()[i];

      CHECK( pc[0] == Approx( po[0] + cen[0] ) );
      CHECK( pc[1] == Approx( po[1] + cen[1] ) );
      CHECK( pc[2] == Approx( po[2] + cen[2] ) );
    }


    // Single recenter 
    s_origin.recenter( cen );
    s_cen.recenter( org );

    CHECK( s_origin.center() == cen );
    CHECK( s_cen.center()    == org );

    for(size_t i = 0; i < npts; ++i) {
      auto po = s_origin.points()[i];
      auto pc = s_cen.points()[i];
      auto po_ref = so_cpy.points()[i];
      auto pc_ref = sc_cpy.points()[i];

      for(int x = 0; x < 3; ++x) {
        CHECK( po[x] == Approx(pc_ref[x]) );
        CHECK( pc[x] == Approx(po_ref[x]) );
      }

    }
    
    // Multiple recenter
    s_origin.recenter( org );
    s_cen.recenter( cen );

    CHECK( s_origin.center() == org );
    CHECK( s_cen.center()    == cen );

    for(size_t i = 0; i < npts; ++i) {
      auto po = s_origin.points()[i];
      auto pc = s_cen.points()[i];
      auto po_ref = so_cpy.points()[i];
      auto pc_ref = sc_cpy.points()[i];

      for(int x = 0; x < 3; ++x) {
        CHECK( po[x] == Approx(po_ref[x]) );
        CHECK( pc[x] == Approx(pc_ref[x]) );
      }
    }


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

#if 1
  SECTION("Batching") {

    size_t max_batch_sz = 512;
    IntegratorXX::MuraKnowles<double,double> r(nrad);
    IntegratorXX::LebedevLaikov<double> q(nang);
    IntegratorXX::SphericalQuadrature s( r, q );

    size_t npts = s.npts();
    //IntegratorXX::SphericalMicroBatcher batcher = 
    //  make_batcher( max_batch_sz, s );

    IntegratorXX::SphericalMicroBatcher batcher = 
      make_batcher( max_batch_sz, s );
  
    size_t npts_c = 0;
    for( auto&& [box_lo, box_up, points_b, weights_b] : batcher ) {

      auto npts_b = points_b.size();
      CHECK( npts_b != 0 );
      npts_c += npts_b;

    }

    npts_c = 0;
    for( size_t i = 0; i < batcher.nbatches(); ++i ) {

      auto&& [box_lo, box_up, points_b, weights_b] = batcher.at(i);

      auto npts_b = points_b.size();
      CHECK( npts_b != 0 );
      npts_c += npts_b;

    }

    CHECK( npts_c == npts );

    const auto& cbatcher = batcher;
    npts_c = 0;
    for( size_t i = 0; i < cbatcher.nbatches(); ++i ) {

      auto&& [box_lo, box_up, points_b, weights_b] = cbatcher.at(i);

      auto npts_b = points_b.size();
      CHECK( npts_b != 0 );
      npts_c += npts_b;

    }

    CHECK( npts_c == npts );

    auto batcher_clone = batcher.clone();
    CHECK( &batcher.quadrature() != &batcher_clone.quadrature() );
    CHECK( batcher.npts() == batcher_clone.npts() );
    CHECK( batcher.nbatches() == batcher_clone.nbatches() );
  
    CHECK( batcher.points() == batcher_clone.points() );
    CHECK( batcher.weights() == batcher_clone.weights() );
  
    batcher.quadrature().recenter({0.,0.,0.}); // just check if this compiles....
  }
#endif
};
