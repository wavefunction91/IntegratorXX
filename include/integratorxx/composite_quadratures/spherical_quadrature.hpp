#pragma once

#include <integratorxx/composite_quadratures/product_quadrature.hpp>
#include <integratorxx/types.hpp>
#include <integratorxx/type_traits.hpp>

namespace IntegratorXX {

namespace detail {

template <typename RealType>
struct spherical_combine_op {

  cartesian_pt_t<RealType> 
    operator()( const RealType& r, const cartesian_pt_t<RealType>& a ) {
    return { r*a[0], r*a[1], r*a[2] };    
  }

  RealType operator()( const RealType& x, const RealType& y ) {
    return x*y;
  }

};

}


template <typename PointContainer, typename WeightContainer >
class SphericalQuadratureBase { 

public:

  using point_container   = PointContainer;
  using weight_container  = WeightContainer;

  using weight_type       = typename weight_container::value_type;
  using point_type        = typename point_container::value_type;

protected:

  point_type center_;

  virtual const point_container& sph_points_adaptor() const = 0;
  virtual       point_container& sph_points_adaptor()       = 0;
  virtual const weight_container& sph_weights_adaptor() const = 0;
  virtual       weight_container& sph_weights_adaptor()       = 0;

public:

  SphericalQuadratureBase( point_type center ) : center_(center) { }
  virtual ~SphericalQuadratureBase() noexcept = default;
  
  auto center() const { return center_; }

  inline void recenter( point_type new_center ) {

    if( new_center != center_ ) {
      point_type shift = new_center;
      shift[0] -= center_[0];
      shift[1] -= center_[1];
      shift[2] -= center_[2];

      auto& points_ = points();
      size_t npts = points_.size();
      for( size_t i = 0; i < npts; ++i ) {
        auto& p = points_[i];
        p[0] += shift[0];
        p[1] += shift[1];
        p[2] += shift[2];
      }

      center_ = new_center;
    }

  }

  const point_container& points() const { return sph_points_adaptor(); }
        point_container& points()       { return sph_points_adaptor(); }
  const weight_container& weights() const { return sph_weights_adaptor(); }
        weight_container& weights()       { return sph_weights_adaptor(); }

  
  size_t npts() const { return points().size(); }

  virtual std::shared_ptr<SphericalQuadratureBase> clone() const = 0;

};















template <typename RadialQuad, typename AngularQuad>
class SphericalQuadrature : 
  public Quadrature<SphericalQuadrature<RadialQuad,AngularQuad>>,
  public SphericalQuadratureBase<
    typename AngularQuad::point_container,
    typename AngularQuad::weight_container
  > {


  using self_type = SphericalQuadrature<RadialQuad,AngularQuad>;
  using traits    = quadrature_traits<self_type>;

  using sph_base = SphericalQuadratureBase<
    typename AngularQuad::point_container,
    typename AngularQuad::weight_container
  >;


public:

  using quad_base_type = Quadrature<SphericalQuadrature<RadialQuad,AngularQuad>>;

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

  SphericalQuadrature( 
    const Quadrature<RadialQuad>& rq, 
    const Quadrature<AngularQuad>& aq, 
    const point_type cen = point_type({0., 0., 0.}) 
  ) :
    quad_base_type( rq, aq, cen ),
    sph_base(cen) { }


  SphericalQuadrature( const SphericalQuadrature& ) = default;
  SphericalQuadrature( SphericalQuadrature&& ) noexcept = default;


  const auto& points() const { return quad_base_type::points(); }
        auto& points()       { return quad_base_type::points(); }
  const auto& weights() const { return quad_base_type::weights(); }
        auto& weights()       { return quad_base_type::weights(); }

  size_t npts() const { return quad_base_type::npts(); }


  std::shared_ptr<sph_base> clone() const override {
    return std::make_shared<self_type>( *this );
  }

}; 





template <typename RadialQuad, typename AngularQuad>
struct quadrature_traits<
  SphericalQuadrature<RadialQuad, AngularQuad>
> {

  using prod_type  = ProductQuadrature<
    detail::spherical_combine_op<typename RadialQuad::point_type>,
    RadialQuad,
    AngularQuad
  >; 
  using prod_traits = quadrature_traits< prod_type >;

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



  inline static std::tuple<point_container,weight_container>
    generate( 
      const Quadrature<RadialQuad>&  r, 
      const Quadrature<AngularQuad>& a,
      const point_type&              cen
    ) {

    auto [points, weights] = prod_traits::generate(r, a); 

    const auto npr  = r.npts();
    const auto npa  = a.npts();

    const auto& rp = r.points();

    // Include Spherical Jacobian
    for( size_t j = 0; j < npa; ++j )
    for( size_t i = 0; i < npr; ++i ) {

      const auto ij = i + j * npr;
      weights[ij] *= 4 * M_PI * rp[i] * rp[i];

    }

    // Recenter
    shift_grid( points, cen );

    return std::make_tuple( points, weights );

  }


};

}
