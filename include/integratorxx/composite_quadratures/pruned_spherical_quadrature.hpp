#pragma once

#include <integratorxx/composite_quadratures/spherical_quadrature.hpp>
#include <integratorxx/composite_quadratures/sub_quadrature.hpp>
#include <integratorxx/types.hpp>
#include <integratorxx/type_traits.hpp>

namespace IntegratorXX {

template <typename RadialQuad, typename AngularQuad>
class PrunedSphericalQuadrature : 
  public Quadrature<PrunedSphericalQuadrature<RadialQuad,AngularQuad>>,
  public SphericalQuadratureBase<
    typename AngularQuad::point_container,
    typename AngularQuad::weight_container
  > {


  using self_type = PrunedSphericalQuadrature<RadialQuad,AngularQuad>;
  using traits    = quadrature_traits<self_type>;

  using sph_base = SphericalQuadratureBase<
    typename AngularQuad::point_container,
    typename AngularQuad::weight_container
  >;

public:

  using quad_base_type = Quadrature<self_type>;

  using point_type       = typename traits::point_type;
  using weight_type      = typename traits::weight_type;
  using point_container  = typename traits::point_container;
  using weight_container = typename traits::weight_container;

protected:

  const point_container& sph_points_adaptor() const override { 
    return quad_base_type::points();
  }
  point_container& sph_points_adaptor() override {;
    return quad_base_type::points();
  }
  const weight_container& sph_weights_adaptor() const override {;
    return quad_base_type::weights();
  }
  weight_container& sph_weights_adaptor() override {;
    return quad_base_type::weights();
  }

public:

  using radial_index_range = std::pair<size_t,size_t>;
  template <typename _Q>
  using angular_subquadrature = 
    std::pair<radial_index_range, _Q>;


  PrunedSphericalQuadrature( 
    const Quadrature<RadialQuad>& rq, 
    const std::vector<angular_subquadrature<AngularQuad>>& aqs,
    const point_type cen = point_type({0., 0., 0.}) 
  ) :
    quad_base_type( rq, aqs, cen ),
    sph_base(cen) { }


  PrunedSphericalQuadrature( const PrunedSphericalQuadrature& ) = default;
  PrunedSphericalQuadrature( PrunedSphericalQuadrature&& ) noexcept = default;


  const auto& points() const { return quad_base_type::points(); }
        auto& points()       { return quad_base_type::points(); }
  const auto& weights() const { return quad_base_type::weights(); }
        auto& weights()       { return quad_base_type::weights(); }

  size_t npts() const { return quad_base_type::npts(); }


  std::shared_ptr<sph_base> clone() const override {
    return std::make_shared<self_type>( *this );
  }

  //inline void recenter( point_type new_center ) {
  //  sph_base::recenter( new_center, this->points_ );
  //}


}; 



template <typename T,size_t N>
std::ostream& operator<<( std::ostream& os, const std::array<T,N>& arr ) {
  os << "{";
  for( auto i = 0; i < (N-1); ++i) os << arr[i] << ", ";
  os << arr[N-1] << "}";
  return os;
}


template <typename RadialQuad, typename AngularQuad>
struct quadrature_traits<
  PrunedSphericalQuadrature<RadialQuad, AngularQuad>
> {

  using radial_subquad_type = SubQuadrature<RadialQuad>;
  using prod_type  = ProductQuadrature<
    detail::spherical_combine_op<typename RadialQuad::point_type>,
    radial_subquad_type,
    AngularQuad
  >; 

  using radial_subquad_traits = quadrature_traits< radial_subquad_type>;
  using prod_traits           = quadrature_traits< prod_type >;

  using radial_traits = quadrature_traits< RadialQuad >;
  using radial_point_container = typename radial_traits::point_container;

  using point_type       = typename prod_traits::point_type;
  using weight_type      = typename prod_traits::weight_type;
  using point_container  = typename prod_traits::point_container;
  using weight_container = typename prod_traits::weight_container;


  inline static void shift_grid( point_container& points, point_type vector ) {

    size_t npts = points.size();
    for( size_t i = 0; i < npts; ++i ) {
      auto& p = points[i];
      p[0] += vector[0];
      p[1] += vector[1];
      p[2] += vector[2];
    }

  }


  using radial_index_range = std::pair<size_t,size_t>;
  using angular_subquadrature = 
    std::pair< radial_index_range, AngularQuad >;

  inline static std::tuple<point_container,weight_container>
    generate( 
      const Quadrature<RadialQuad>&             r, 
      const std::vector<angular_subquadrature>& as,
      const point_type&                         cen
    ) {

    point_container  points;
    weight_container weights;
#if 1
    // Loop over angular quadratures + associated ranges
    for( const auto& [r_range, a] : as ) {

      // Generate radial subquadrature
      radial_subquad_type r_subquad( r_range, r );

      // Generate product subquadrature
      auto [sub_points, sub_weights] = prod_traits::generate( r_subquad, a );

      // Include Spherical Jacobian
      const auto npr  = r_subquad.npts();
      const auto npa  = a.npts();

      const auto& rp = r_subquad.points();

      // Include Spherical Jacobian
      for( size_t j = 0; j < npa; ++j )
      for( size_t i = 0; i < npr; ++i ) {

        const auto ij = i + j * npr;
        sub_weights[ij] *= 4 * M_PI * rp[i] * rp[i];

      }

      // Insert points/weights into final storage
      points.insert( points.end(), sub_points.begin(), sub_points.end() );
      weights.insert( weights.end(), sub_weights.begin(), sub_weights.end() );
    }
#endif


    // Recenter
    shift_grid( points, cen );

    return std::make_tuple( points, weights );

  }


};

}
