#include <functional>
#include <type_traits>
#include <integratorxx/quadrature.hpp>
#include <integratorxx/batch.hpp>
#include <integratorxx/integrator.hpp>

#include <algorithm>
#include <numeric>
#include <type_traits>
#include <iostream>

template <typename T> class type_printer;

int main() {

  using namespace IntegratorXX;

  GaussLegendre<double> quad(10,-1,1);
  GaussChebyshev1<double> quad2(10,-1,1);
  GaussChebyshev2<double> quad3(10,-1,1);
  EulerMaclaurin<double> quad4(100);
  Lebedev<double> quad5(302);

  std::cout << quad4.points().front() << ", " << quad4.points().back() << std::endl;
  std::cout << quad4.weights().front() << ", " << quad4.weights().back() << std::endl;


  std::cout << (quad4.points() | ranges::view::take(100)) << std::endl;

  // Functions / reference integrands
  constexpr auto gaussian = 
    []( double x ){ return x*x*std::exp(-x*x); };
  constexpr double ref_gaussian = std::sqrt( M_PI ) / 4.;

  // 1-D Quadratue
  auto quad_range = ranges::view::zip( quad4.points(), quad4.weights() );


  double ret_1 = 0.;
  for( auto&& [pt,wgt] : quad_range ) // Loop over quadrature point-wise
    ret_1 += wgt * gaussian( pt );

  std::cout << std::abs( ret_1 - ref_gaussian ) << std::endl;



  // Chunk over zipped range
  auto batched_1d_range = quad_range | ranges::view::chunk(4);


  ret_1 = 0.;
  for( auto&& q : batched_1d_range ) 
  for( auto&& [pt, wgt] : q ) 
    ret_1 += wgt * gaussian( pt );

  std::cout << std::abs( ret_1 - ref_gaussian ) << std::endl;




  // Zip over chunked range
  auto batched_1d_pts   = quad4.points()  | ranges::view::chunk(4);
  auto batched_1d_wgt   = quad4.weights() | ranges::view::chunk(4);
  auto batched_range_1d = ranges::view::zip( batched_1d_pts, batched_1d_wgt );

  ret_1 = 0.;
  for( auto&& [pts, wgts] : batched_range_1d ){ 
    auto x = pts[0]; // XXX: Is it lazy up to this point? Or does this force evaluation of pts?
  for( auto&& [pt, wgt]   : ranges::view::zip( pts, wgts ) )
    ret_1 += wgt * gaussian( pt );
  }
  std::cout << std::abs( ret_1 - ref_gaussian ) << std::endl;



  using eval_type = std::tuple< std::vector<double>, std::vector<double> >;

  ret_1 = 0.;
  for( eval_type&& batch : batched_range_1d ) {
    auto&& [pts, wgts] = std::move(batch); // Force evaluation of range into normal containers
  for( auto&& [pt, wgt]   : ranges::view::zip( pts, wgts ) )
    ret_1 += wgt * gaussian( pt );
  }

  std::cout << std::abs( ret_1 - ref_gaussian ) << std::endl;
  


  // 2-D Quadrature

  // Transform 2D quadrature into something more natural ( R^2 x R )
  constexpr auto flatten_2d_grid = []( const auto& y ) {
    const auto& [q1, q2] = y;
    const auto& [p1, w1] = q1;
    const auto& [p2, w2] = q2;
    return std::make_tuple( p1, p2, w1*w2 );
  };

  // Pointwise cartesian product
  auto quad_range_2d = 
    ranges::view::cartesian_product( quad_range, quad_range );

  double ret_2 = 0.;
  size_t eval_cnt = 0;
  for( auto&& [q1,q2] : quad_range_2d ) {
    const auto& [pt1, wgt1] = q1;
    const auto& [pt2, wgt2] = q2;

    const double wgt  = wgt1 * wgt2;
    const double eval = gaussian( pt1 ) * gaussian( pt2 );

    ret_2 += wgt * eval;
    ++eval_cnt;
  }

  std::cout << std::abs( ret_2 - std::pow(ref_gaussian,2) )  << ", " << eval_cnt
            << std::endl;







  // Flatten from constructed range
  auto flat_quad_range_2d = quad_range_2d |
    ranges::view::transform( flatten_2d_grid );

  ret_2 = 0.;
  eval_cnt = 0;
  for( auto&& [p1, p2, w] : flat_quad_range_2d ) {
    ret_2 += w * gaussian(p1) * gaussian(p2);
    ++eval_cnt;
  }
  
  std::cout << std::abs( ret_2 - std::pow(ref_gaussian,2) )  << ", " << eval_cnt
            << std::endl;





  // Flatten directly
  auto flat_quad_range_2d_direct = 
    ranges::view::cartesian_product( quad_range, quad_range ) |
    ranges::view::transform( flatten_2d_grid );

  ret_2 = 0.;
  eval_cnt = 0;
  for( auto&& [p1, p2, w] : flat_quad_range_2d_direct ){
    ret_2 += w * gaussian(p1) * gaussian(p2);
    ++eval_cnt;
  }
  
  std::cout << std::abs( ret_2 - std::pow(ref_gaussian,2) )  << ", " << eval_cnt
            << std::endl;




  // Chunk of flattened 2d quadrature
  auto batched_range_2d = flat_quad_range_2d_direct | ranges::view::chunk(4);
  ret_2 = 0.;
  eval_cnt = 0;
  for( auto&& batch : batched_range_2d )
  for( auto&& [p1, p2, w] : batch ){
    ret_2 += w * gaussian(p1) * gaussian(p2);
    ++eval_cnt;
  }

  std::cout << std::abs( ret_2 - std::pow(ref_gaussian,2) )  << ", " << eval_cnt
            << std::endl;



/*
  // Compile error
  auto batched_range_2d = 
    ranges::view::cartesian_product( batched_1d_range, batched_1d_range );

  ret_2 = 0.;
  for( auto&& [q1,q2] : batched_range_2d ) {
    const auto& [pt1_batch, wgt1_batch] = q1;
    const auto& [pt2_batch, wgt2_batch] = q2;

    for( auto [pt1, pt2, wgt1, wgt2] : ranges::view::zip(pt1_batch, pt2_batch, wgt1_batch, wgt2_batch) ) {
      const double wgt  = wgt1 * wgt2;
      const double eval = gaussian( pt1 ) * gaussian( pt2 );

      ret_2 += wgt * eval;
    }
  }

  std::cout << std::abs( ret_2 - std::pow(ref_gaussian,2) ) 
            << std::endl;
*/


  // This is really gross...
  auto batched_range_2d_2 = ranges::view::cartesian_product( batched_range_1d, batched_range_1d );

  ret_2 = 0.;
  eval_cnt = 0;
 
  for( auto&& [q1, q2] : batched_range_2d_2 ){

    auto&& [pt1_batch, wgt1_batch] = q1;
    auto&& [pt2_batch, wgt2_batch] = q2;
  

    auto pts_batch = ranges::view::cartesian_product( pt1_batch, pt2_batch );
    auto wgts_batch = ranges::view::cartesian_product( wgt1_batch, wgt2_batch ) | 
                      ranges::view::transform([](auto z){ return std::get<0>(z) * std::get<1>(z); });

    for( auto&& [qq1, wgt] : ranges::view::zip(pts_batch, wgts_batch) ) {
      auto&& [pt1, pt2]   = qq1;
      const double eval = gaussian( pt1 ) * gaussian( pt2 );

      ret_2 += wgt * eval;
      ++eval_cnt;
    }

  }

  std::cout << std::abs( ret_2 - std::pow(ref_gaussian,2) )  << ", " << eval_cnt
            << std::endl;




  auto flat_batched_range_2d = batched_range_2d_2 |
    ranges::view::transform( [](auto&& z) {
      auto&& [q1, q2] = std::move(z);
      auto&& [pts1, wgts1] = std::move( q1 );
      auto&& [pts2, wgts2] = std::move( q2 );

      auto pts_batch = ranges::view::cartesian_product( pts1, pts2 ) |
                       ranges::view::transform( []( auto&& t ) { 
                         return std::array{std::get<0>(t), std::get<1>(t)}; 
                       } );

      auto wgts_batch = ranges::view::cartesian_product( wgts1, wgts2 ) |
                        ranges::view::transform( []( auto&& t ) { 
                          return std::get<0>(t) * std::get<1>(t); 
                        } );

      return std::make_tuple( pts_batch, wgts_batch );
    });

  ret_2 = 0.;
  eval_cnt = 0;
  for( auto&& [pts, wgts] : flat_batched_range_2d ) {
//  std::vector< std::array<double,2> > z = pts;
//  std::cout << z.size() << ", " << z.data() << std::endl;
  for( auto&& [pt, wgt]   : ranges::view::zip(pts, wgts) ) {
    ret_2 += wgt * gaussian(pt[0]) * gaussian(pt[1]);
    eval_cnt++;
  }
  }
  std::cout << std::abs( ret_2 - std::pow(ref_gaussian,2) )  << ", " << eval_cnt
            << std::endl;

  return 0;
}
