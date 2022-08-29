#pragma once

#include <integratorxx/quadrature.hpp>
#include <integratorxx/util/bound_transform.hpp>

namespace IntegratorXX {


template <typename PointType, typename WeightType>
class GaussLegendre : 
  public Quadrature<GaussLegendre<PointType,WeightType>> {

  using base_type = Quadrature<GaussLegendre<PointType,WeightType>>;

public:

  using point_type       = typename base_type::point_type;
  using weight_type      = typename base_type::weight_type;
  using point_container  = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;
  
  GaussLegendre(size_t npts, point_type lo, point_type up):
    base_type( npts, lo, up ) { }

  GaussLegendre( const GaussLegendre& ) = default;
  GaussLegendre( GaussLegendre&& ) noexcept = default;
};







template <typename PointType, typename WeightType>
struct quadrature_traits<
  GaussLegendre<PointType,WeightType>
> {

  using point_type  = PointType;
  using weight_type = WeightType;

  using point_container  = std::vector< point_type >;
  using weight_container = std::vector< weight_type >;

  inline static std::tuple<point_container,weight_container>
    generate( size_t npts, point_type lo, point_type up ) {

    assert( npts % 2 == 0 );
    assert( lo != -std::numeric_limits<double>::infinity() );
    assert( up != std::numeric_limits<double>::infinity() );

    point_container  points( npts );
    weight_container weights( npts );

    const size_t mid = (npts + 1) / 2;
    const double eps = 3.e-11; // Convergence tolerance

    for( size_t i = 1; i <= mid; ++i ) {

      const point_type ip = i;

      point_type z = std::cos( M_PI * (ip - 0.25) / (npts + 0.5));
      point_type pp(0), z1(0);

      // Iteratively determine the i-th root 
      while(std::abs(z-z1) > eps) {
        point_type p1(1.), p2(0.);
  
        // Loop over the recurrence relation to evaluate the
        // Legendre polynomial at position z
        for( size_t j = 1; j <= npts; ++j ){
          const point_type jp = j;

          point_type p3 = p2;
          p2 = p1; 
          p1 = ( (2. * jp - 1.) * z * p2 - (jp - 1.) * p3) / jp;
        } // end j for
  
        // p1 is now the desired Legrendre polynomial. We next compute 
        // pp, its derivative, by a standard relation involving also p2,
        // the polynomial of one lower order
        pp = npts * (z * p1 - p2) / (z * z - 1.);
        z1 = z;
        z  = z1 - p1 / pp;
      } // end while
  
      point_type  pt = z;
      weight_type wgt = 2. / (1. - z*z) / pp / pp;

      // Transform points and populate arrays
      std::tie(pt,wgt) = transform_minus_one_to_one(lo,up,pt,wgt);
      points[i-1] = pt;
      weights[i-1] = wgt;

      // Reflect the points
      points[npts-i]  = (lo + up) - pt;
      weights[npts-i] = wgt;

    }; // Loop over points

    return std::make_tuple( points, weights );

  }

};

}
