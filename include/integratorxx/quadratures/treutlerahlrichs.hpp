#pragma once

#include <integratorxx/quadrature.hpp>
#include <integratorxx/quadratures/gausscheby2.hpp>

namespace IntegratorXX {


/**
 *  @brief Implementation of the Treutler-Ahlrichs radial quadrature.
 *
 *  Generates a quadrature on the bounds (0, inf). Suitable for integrands
 *  which tend to zero as their argument tends to 0 and inf. Tailored for
 *  radial integrands, i.e. r^2 * f(r), with lim_{r->inf} f(r) = 0.
 *
 *  Reference:
 *  J. Chem. Phys. 102, 346 (1995)
 *  DOI: https://doi.org/10.1063/1.469408
 *
 *  @tparam PointType  Type describing the quadrature points
 *  @tparam WeightType Type describing the quadrature weights
 */
template <typename PointType, typename WeightType>
class TreutlerAhlrichs :
  public Quadrature<TreutlerAhlrichs<PointType,WeightType>> {

  using base_type = Quadrature<TreutlerAhlrichs<PointType,WeightType>>;

public:

  using point_type       = typename base_type::point_type;
  using weight_type      = typename base_type::weight_type;
  using point_container  = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;

  /**
   *  @brief Construct the Treutler-Ahlrichs radial quadrature
   *
   *  @param[in] npts   Number of quadrature points to generate
   *  @param[in] R      Radial scaling factor. Table for suggested
   *                    values is given in the original reference.
   *  @param[in] alpha  Exponent factor for the quadrature. 0.6 was the
   *                    default suggested in the reference.
   */
  TreutlerAhlrichs(size_t npts, weight_type R = 1., weight_type alpha = 0.6):
    base_type( npts, R, alpha ) { }

  TreutlerAhlrichs( const TreutlerAhlrichs& )     = default;
  TreutlerAhlrichs( TreutlerAhlrichs&& ) noexcept = default;
};






/**
 *  @brief Quadrature traits for the Treutler-Ahlrichs quadrature
 *
 *  @tparam PointType  Type describing the quadrature points
 *  @tparam WeightType Type describing the quadrature weights
 */

template <typename PointType, typename WeightType>
struct quadrature_traits<
  TreutlerAhlrichs<PointType,WeightType>
> {

  using point_type  = PointType;
  using weight_type = WeightType;

  using point_container  = std::vector< point_type >;
  using weight_container = std::vector< weight_type >;

  /**
   *  @brief Generator for the Treutler-Ahlrichs quadrature
   *
   *  @param[in] npts   Number of quadrature points to generate
   *  @param[in] R      Radial scaling factor. Table for suggested
   *                    values is given in the original reference.
   *  @param[in] alpha  Exponent factor for the quadrature. 0.6 was the
   *
   *  @returns Tuple of quadrature points and weights
   */
  inline static std::tuple<point_container,weight_container>
    generate( size_t npts, weight_type R, weight_type alpha ) {


    const point_type ln_2 = std::log(2.);
    const point_type pi_ov_npts_p1 = M_PI / (npts + 1);

    point_container  points( npts );
    weight_container weights( npts );

    /*
     * Treutler-Ahlrichs quadrature
     *
     * Original reference:
     *  J. Chem. Phys. 102, 346 (1995)
     *  DOI: https://doi.org/10.1063/1.469408
     *
     * Closed form for points / weights obtained from:
     * Journal of Computational Chemistry, 24: 732–740, 2003
     * DOI: https://doi.org/10.1002/jcc.10211
     *
     */
    for( size_t i = 0; i < npts; ++i ) {
      const auto xi = std::cos( (i+1) * pi_ov_npts_p1 );

      const auto pow_term  = std::pow( 1. + xi, alpha ) / ln_2;
      const auto log_term  = std::log( (1. - xi)/2. );
      const auto sqrt_term = std::sqrt( (1.+xi)/(1.-xi) );

      points[i]  = -R * pow_term * log_term;
      weights[i] = R * pi_ov_npts_p1 * pow_term *
        ( sqrt_term - alpha * log_term / sqrt_term );
    }

    return std::make_tuple( points, weights );

  }

};





struct TreutlerAhlrichsRadialTraits {

  inline static constexpr double alpha = 0.6;
  inline static constexpr double ln_2  = 0.693147180559945309417232;

  template <typename PointType>
  static auto radial_transform(PointType x) {
    const auto pow_term = std::pow(1.0 + x, alpha);
    const auto log_term = std::log(0.5 * (1.0 - x));
    return  pow_term * log_term / ln_2; 
  };


  template <typename PointType>
  static auto radial_jacobian(PointType x) {
    const auto pow_term = std::pow(1.0 + x, alpha);
    const auto log_term = std::log(0.5 * (1.0 - x));
    return pow_term / ln_2 * (-alpha * log_term / (1.0 + x) + 1.0 / (1.0 - x));
  }

};



}
