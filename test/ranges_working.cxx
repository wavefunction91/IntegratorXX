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


  EulerMaclaurin<double> quad(100); // Generate quadrature

  // Functions / reference integrands
  constexpr auto gaussian = 
    []( double x ){ return x*x*std::exp(-x*x); };

  constexpr auto gaussian_2d = 
    [=]( double x, double y ){ return gaussian(x) * gaussian(y); };

  constexpr double ref_gaussian = std::sqrt( M_PI ) / 4.;
  constexpr double ref_gaussian_2d =  M_PI  / 16.;


  // Zip over chunked range
  auto batched_1d_pts   = quad.points()  | ranges::view::chunk(4);
  auto batched_1d_wgt   = quad.weights() | ranges::view::chunk(4);
  auto batched_range_1d = ranges::view::zip( batched_1d_pts, batched_1d_wgt );


  
  // Loop over batches and accumulate 1D integrand
  double ret_1 = 0.;
  // pts / wgts are also convertable to std::vector<double>
  // internal loop not necessary
  for( auto&& [pts, wgts] : batched_range_1d )
  for( auto&& [pt, wgt]   : ranges::view::zip( pts, wgts ) )
    ret_1 += wgt * gaussian( pt );

  std::cout << std::abs( ret_1 - ref_gaussian ) << std::endl;


  // 2D...
  auto batched_range_2d = 
    ranges::view::cartesian_product( batched_range_1d, 
                                     batched_range_1d );

  auto flat_batched_range_2d = batched_range_2d |
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


  double ret_2 = 0.;
  for( auto&& [pts, wgts] : flat_batched_range_2d ) 
  for( auto&& [pt, wgt]   : ranges::view::zip(pts, wgts) ) 
    ret_2 += wgt * gaussian_2d(pt[0], pt[1]);

  std::cout << std::abs( ret_2 - ref_gaussian_2d ) << std::endl;

  return 0;
}
