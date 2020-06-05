#pragma once

#include <integratorxx/quadrature.hpp>
#include <integratorxx/quadratures/uniform.hpp>

namespace IntegratorXX {

/**
 *  @brief Implementation of the Murray-Handy-Laming radial quadrature.
 *
 *  Generates a quadrature on the bounds (0, inf). Suitable for integrands
 *  which tend to zero as their argument tends to 0 and inf. Tailored for
 *  radial integrands, i.e. r^2 * f(r), with lim_{r->inf} f(r) = 0.
 *
 *  Reference:
 *  Molecular Physics, 78:4, 997-1014,
 *  DOI: https://doi.org/10.1080/00268979300100651
 *
 *  @tparam PointType  Type describing the quadrature points  
 *  @tparam WeightType Type describing the quadrature weights 
 */
template <typename PointType, typename WeightType>
class MurrayHandyLaming : 
  public Quadrature<MurrayHandyLaming<PointType,WeightType>> {

  using base_type = Quadrature<MurrayHandyLaming<PointType,WeightType>>;

public:

  using point_type       = typename base_type::point_type;
  using weight_type      = typename base_type::weight_type;
  using point_container  = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;
  
  /**
   *  @brief Construct the Murray-Handy-Laming radial quadrature
   *
   *  @param[in] npts Number of quadrature points to generate
   *  @param[in] R    Radial scaling factor. Suggested to be set to 
   *                  correspond to the Bragg-Slater radius of the
   *                  appropriate nucleus 
   *  @param[in] m    Exponent factor for the quadrature. 2 was the
   *                  default suggested in the reference. 
   */
  MurrayHandyLaming(size_t npts, weight_type R = 1., int m = 2): 
    base_type( npts, R, m ) { }

  MurrayHandyLaming( const MurrayHandyLaming& )     = default;
  MurrayHandyLaming( MurrayHandyLaming&& ) noexcept = default;
};






/**
 *  @brief Quadrature traits for the Murray-Handy-Laming quadrature
 *
 *  @tparam PointType  Type describing the quadrature points  
 *  @tparam WeightType Type describing the quadrature weights 
 */
template <typename PointType, typename WeightType>
struct quadrature_traits<
  MurrayHandyLaming<PointType,WeightType>,
  std::enable_if_t<
    std::is_floating_point_v<PointType> and
    std::is_floating_point_v<WeightType>
  >
> {

  using point_type  = PointType;
  using weight_type = WeightType;

  using point_container  = std::vector< point_type >;
  using weight_container = std::vector< weight_type >;


  /**
   *  @brief Generator for the Murray-Handy-Laming quadrature
   *
   *  @param[in] npts Number of quadrature points to generate
   *  @param[in] R    Radial scaling factor. Suggested to be set to 
   *                  correspond to the Bragg-Slater radius of the
   *                  appropriate nucleus 
   *  @param[in] m    Exponent fector for the quadrature. 2 was the
   *                  default suggested in the original reference. 
   *
   *  @returns Tuple of quadrature points and weights
   */
  inline static std::tuple<point_container,weight_container>
    generate( size_t npts, weight_type R, int m ) {

    point_container  points( npts );
    weight_container weights( npts );


    using base_quad_traits = 
      quadrature_traits<UniformTrapezoid<PointType,WeightType>>;

    /*
     * Generate uniform trapezoid points on [0,1]
     * ux(j) = j/(m-1)
     * uw(j) = 1/(m-1)
     * m = npts + 2, j in [0, npts+2)
     */
    auto [ux, uw] = base_quad_traits::generate( npts+2, 0., 1. );

    /*
     * Perform Murray, Handy and Laming transformation
     *
     * Original reference:
     * Molecular Physics, 78:4, 997-1014,
     * DOI: https://doi.org/10.1080/00268979300100651
     *
     * Closed form for points / weights obtained from:
     * Journal of Computational Chemistry, 24: 732–740, 2003
     * DOI: https://doi.org/10.1002/jcc.10211
     * 
     * x(i) = ux(i+1)
     * r(i) = [x(i) / (1 - x(i))]^m
     * w(i) = uw(i+1) * m * x(i)^(m-1) * [1 - x(i)]^(-m-1)
     * i in [0, npts)
     *
     * XXX: i+1 offset on trapezoid points ignores enpoints of trapezoid 
     *      quadrature (f(r) = 0 with r in {0,inf})
     */
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
