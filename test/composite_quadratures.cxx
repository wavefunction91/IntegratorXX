#include "catch2/catch.hpp"

#include <integratorxx/quadratures/muraknowles.hpp>
#include <integratorxx/quadratures/lebedev_laikov.hpp>
#include <integratorxx/composite_quadratures/spherical_quadrature.hpp>
#include <integratorxx/composite_quadratures/pruned_spherical_quadrature.hpp>
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
      
}



TEST_CASE( "Pruned Spherical Quadratures", "[sph-quad]" ) {

     
  SECTION("Consistency With Spherical Quadrature") {
    size_t nrad = 4;
    size_t nang = 6;

    IntegratorXX::MuraKnowles<double,double> r(nrad);
    IntegratorXX::LebedevLaikov<double> q(nang);
    IntegratorXX::RadialGridPartition<decltype(q)> rgp;

    SECTION("Single Radial Partition") {
      rgp = IntegratorXX::RadialGridPartition( r, 0, q );
    }
    SECTION("Multiple Radial Partition") {
      rgp = IntegratorXX::RadialGridPartition( r, 0, q, r.npts()/2, q );
    }

    IntegratorXX::SphericalQuadrature       s( r, q );
    IntegratorXX::PrunedSphericalQuadrature ps( r, rgp );
    REQUIRE( s.npts() == ps.npts() );

    // Impose lexiographical ordering for check
    // b/c pruning doesn't have radial points as globally fastest running
    std::vector<size_t> idx_s( s.npts() ), idx_ps( ps.npts() );
    std::iota( idx_s.begin(), idx_s.end(), 0 );
    std::iota( idx_ps.begin(), idx_ps.end(), 0 );

    auto comparator_pt = [](const auto& _q, auto i, auto j) {
      return _q.points()[i] < _q.points()[j];
    };
    auto comparator_s = [&](auto i, auto j) {
      return comparator_pt(s,i,j);
    };
    auto comparator_ps = [&](auto i, auto j) {
      return comparator_pt(ps,i,j);
    };

    std::sort( idx_s.begin(), idx_s.end(), comparator_s );
    std::sort( idx_ps.begin(), idx_ps.end(), comparator_ps );

    for( int i = 0; i < s.npts(); ++i ) {
      CHECK( s.points()[idx_s[i]] == ps.points()[idx_ps[i]] );
      CHECK( s.weights()[idx_s[i]] == ps.weights()[idx_ps[i]] );
    }

  }
      

  SECTION("Validity") {

    size_t nrad = 12;
    IntegratorXX::MuraKnowles<double,double> r(nrad);
    IntegratorXX::LebedevLaikov<double> 
      q_6(6), q_14(14), q_302(302), q_770(770);

    IntegratorXX::RadialGridPartition 
      rgp( r, 0ul, q_6, 4ul, q_302, 6ul, q_14, 9ul, q_770 );
    IntegratorXX::PrunedSphericalQuadrature ps( r, rgp );

    size_t npts = 0;
    for( const auto& [r_range, q] : rgp ) {
      auto [r1,r2] = r_range;
      for( auto j = 0;  j < q.npts(); ++j )
      for( auto i = r1; i < r2;       ++i ) {
        auto pt = ps.points()[npts];
        CHECK( pt[0] == r.points()[i] * q.points()[j][0] );
        CHECK( pt[1] == r.points()[i] * q.points()[j][1] );
        CHECK( pt[2] == r.points()[i] * q.points()[j][2] );

        const auto _r = r.points()[i];
        const auto _rw = r.weights()[i];
        const auto _aw = q.weights()[j];
        CHECK( ps.weights()[npts] == Approx( 4 * M_PI * _r * _r * _rw * _aw ) );
        npts++;
      }
    }
    CHECK( npts == ps.npts() );
      
  }

  SECTION("With Batcher") {
    IntegratorXX::MuraKnowles<double,double> r(4);
    IntegratorXX::LebedevLaikov<double>      q(302);
    IntegratorXX::RadialGridPartition        rgp( r, 0ul, q );

    IntegratorXX::SphericalQuadrature       s(r,q);
    IntegratorXX::PrunedSphericalQuadrature ps(r,rgp);

    IntegratorXX::SphericalMicroBatcher batcher_s =
      make_batcher( 1, s );
    IntegratorXX::SphericalMicroBatcher batcher_ps =
      make_batcher( 1, ps );

    REQUIRE( batcher_s.nbatches() == batcher_ps.nbatches() );

    for( auto i = 0ul; i < batcher_s.nbatches(); ++ i ) {
      auto [s_box_lo,  s_box_up,  s_points,  s_weights ] = batcher_s.at(i);
      auto [ps_box_lo, ps_box_up, ps_points, ps_weights] = batcher_ps.at(i);

      CHECK( s_box_lo == ps_box_lo );
      CHECK( s_box_up == ps_box_up );
      CHECK( s_points == ps_points );
      CHECK( s_weights == ps_weights );
    }
  }

}








