#pragma once

#include <integratorxx/quadrature.hpp>
#include <integratorxx/quadratures/uniform.hpp>

namespace IntegratorXX {


template <typename PointType, typename WeightType>
class MurrayHandyLaming : 
  public Quadrature<MurrayHandyLaming<PointType,WeightType>> {

  using base_type = Quadrature<MurrayHandyLaming<PointType,WeightType>>;

public:

  using point_type       = typename base_type::point_type;
  using weight_type      = typename base_type::weight_type;
  using point_container  = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;
  
  MurrayHandyLaming(size_t npts, weight_type R = 1., int m = 2): base_type( npts, R, m ) { }

  MurrayHandyLaming( const MurrayHandyLaming& ) = default;
  MurrayHandyLaming( MurrayHandyLaming&& ) noexcept = default;
};







template <typename PointType, typename WeightType>
struct quadrature_traits<
  MurrayHandyLaming<PointType,WeightType>
> {

  using point_type  = PointType;
  using weight_type = WeightType;

  using point_container  = std::vector< point_type >;
  using weight_container = std::vector< weight_type >;

  inline static std::tuple<point_container,weight_container>
    generate( size_t npts, weight_type R, int m ) {

    point_container  points( npts );
    weight_container weights( npts );


    using base_quad_traits = 
      quadrature_traits<UniformTrapezoid<PointType,WeightType>>;

    // Generate uniform trapezoid points on [0,1]
    // ux(j) = j/(m-1)
    // uw(j) = 1/(m-1)
    // m = npts + 2, j in [0, npts+2)
    auto [ux, uw] = base_quad_traits::generate( npts+2, 0., 1. );

    // Perform Murray, Handy and Laming transformation
    // Molecular Physics, 78:4, 997-1014,
    // DOI: https://doi.org/10.1080/00268979300100651
    // 
    // x(i) = ux(i+1)
    // r(i) = [x(i) / (1 - x(i))]^m
    // w(i) = uw(i+1) * m * x(i)^(m-1) * [1 - x(i)]^(-m-1)
    // i in [0, npts)
    //
    // XXX: i+1 offset on trapezoid points ignores enpoints of trapezoid 
    //      quadrature (f(r) = 0 with r in {0,inf})
    for( size_t i = 0; i < npts; ++i ) {
      const auto xi       = ux[i+1];
      const auto one_m_xi = 1. - xi;
      points[i]  = R * std::pow( xi / one_m_xi, m );
      weights[i] = R * uw[i+1] * m * std::pow( xi, m-1 ) / std::pow( one_m_xi, m+1);
    }
    


    return std::tuple( points, weights );

  }

};

}
