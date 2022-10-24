#include "catch2/catch.hpp"

#include <integratorxx/quadratures/muraknowles.hpp>
#include <integratorxx/quadratures/lebedev_laikov.hpp>
#include <integratorxx/composite_quadratures/spherical_quadrature.hpp>
#include <integratorxx/composite_quadratures/pruned_spherical_quadrature.hpp>
#include <integratorxx/batch/spherical_micro_batcher.hpp>
#include <integratorxx/batch/hilbert_partition.hpp>

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
}

TEST_CASE("Hilbert") {
  SECTION("Exhaustive 3x2 -> uint64_t") {
    std::vector<std::array<uint32_t,2>> Xs = {
      {0,0}, {0,1}, {1,1}, {1,0}, {2,0}, {3,0}, {3,1}, {2,1}, {2,2}, {3,2},
      {3,3}, {2,3}, {1,3}, {1,2}, {0,2}, {0,3}, {0,4}, {1,4}, {1,5}, {0,5},
      {0,6}, {0,7}, {1,7}, {1,6}, {2,6}, {2,7}, {3,7}, {3,6}, {3,5}, {2,5},
      {2,4}, {3,4}, {4,4}, {5,4}, {5,5}, {4,5}, {4,6}, {4,7}, {5,7}, {5,6},
      {6,6}, {6,7}, {7,7}, {7,6}, {7,5}, {6,5}, {6,4}, {7,4}, {7,3}, {7,2},
      {6,2}, {6,3}, {5,3}, {4,3}, {4,2}, {5,2}, {5,1}, {4,1}, {4,0}, {5,0},
      {6,0}, {6,1}, {7,1}, {7,0}
    };
    size_t i = 0;
    for( auto X : Xs ) {
      auto h = IntegratorXX::hilbert_encode<uint64_t,3>(X);
      REQUIRE(h == i++);
    }
  }

  SECTION("Spot check 3D -> uint64_t") {
    // This comes from the Skilling paper
    std::array<uint32_t,3> X = {5, 10, 20};
    auto hh = IntegratorXX::hilbert_encode<uint64_t,5>(X);
    REQUIRE(hh == 7865ul);
  }
}


TEMPLATE_TEST_CASE("SphericalMicroBatching", "[sph-quad]",
  IntegratorXX::OctreeGridPartitioner, IntegratorXX::HilbertGridPartitioner<21>) {

  using grid_partitioner_type = TestType;

  const size_t nrad = 35;
  const size_t nang = 110;

  size_t max_batch_sz = 512;
  IntegratorXX::MuraKnowles<double,double> r(nrad);
  IntegratorXX::LebedevLaikov<double> q(nang);
  IntegratorXX::SphericalQuadrature s( r, q );

  size_t npts = s.npts();

  IntegratorXX::SphericalMicroBatcher batcher = 
    IntegratorXX::make_batcher<grid_partitioner_type>( max_batch_sz, s );

  size_t npts_c = 0;
  for( size_t i = 0; i < batcher.nbatches(); ++i ) {

    auto&& [box_lo, box_up, points_b, weights_b] = batcher.at(i);

    auto npts_b = points_b.size();
    REQUIRE( npts_b != 0 );
    npts_c += npts_b;

    const double vol = 
      (box_up[0] - box_lo[0]) *
      (box_up[1] - box_lo[1]) *
      (box_up[2] - box_lo[2]); 

    //std::cout << "b" << i << " = [" << std::endl;
    //for( auto pt : points_b ) std::cout << pt[0] << ", " << pt[1] << ", " << pt[2] << std::endl;
    //std::cout << "];" << std::endl;
  }
  

  REQUIRE( npts_c == npts );

  auto batcher_clone = batcher.clone();
  CHECK( &batcher.quadrature() != &batcher_clone.quadrature() );
  CHECK( batcher.npts() == batcher_clone.npts() );
  CHECK( batcher.nbatches() == batcher_clone.nbatches() );
  
  CHECK( batcher.points() == batcher_clone.points() );
  CHECK( batcher.weights() == batcher_clone.weights() );
  
  batcher.quadrature().recenter({0.,0.,0.}); // just check if this compiles....

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

    using partitioner_type = IntegratorXX::OctreeGridPartitioner;
    IntegratorXX::SphericalMicroBatcher batcher_s =
      IntegratorXX::make_batcher<partitioner_type>( 1, s );
    IntegratorXX::SphericalMicroBatcher batcher_ps =
      IntegratorXX::make_batcher<partitioner_type>( 1, ps );

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








