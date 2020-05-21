#pragma once

#include <integratorxx_new/quadrature.hpp>
#include <integratorxx_new/quadratures/lebedev_laikov/lebedev_laikov_grids.hpp>

namespace IntegratorXX::redeux {

namespace detail {
  template <class StaticGrid, typename PointContainer, typename WeightsContainer>
  void copy_grid( PointContainer& points, WeightsContainer& weights ){

    const auto& static_points  = StaticGrid::points;
    const auto& static_weights = StaticGrid::weights;

    std::copy( static_points.begin(), static_points.end(), points.begin() );
    std::copy( static_weights.begin(), static_weights.end(), weights.begin() );

  }
}

template <typename RealType>
class LebedevLaikov :
  public Quadrature<LebedevLaikov<RealType>> {

  using base_type = Quadrature<LebedevLaikov<RealType>>;

public:

  using point_type       = typename base_type::point_type;
  using weight_type      = typename base_type::weight_type;
  using point_container  = typename base_type::point_container;
  using weight_container = typename base_type::weight_container;
  
  LebedevLaikov(size_t npts): base_type( npts ) { }

  LebedevLaikov( const LebedevLaikov& ) = default;
  LebedevLaikov( LebedevLaikov&& ) noexcept = default;
};


template <typename RealType>
struct quadrature_traits< LebedevLaikov<RealType> > {

  using point_type  = cartesian_pt_t<RealType>;
  using weight_type = RealType;

  using point_container  = std::vector< point_type >;
  using weight_container = std::vector< weight_type >;

  inline static std::tuple<point_container,weight_container>
    generate( size_t npts ) {

    point_container  points( npts );
    weight_container weights( npts );

    using namespace LebedevLaikovGrids;

    if( npts == 6 ) 
      detail::copy_grid<lebedev_laikov_6<RealType>>( points, weights );
    else if( npts == 14 ) 
      detail::copy_grid<lebedev_laikov_14<RealType>>( points, weights );
    else if( npts == 26 ) 
      detail::copy_grid<lebedev_laikov_26<RealType>>( points, weights );
    else if( npts == 38 ) 
      detail::copy_grid<lebedev_laikov_38<RealType>>( points, weights );
    else if( npts == 50 ) 
      detail::copy_grid<lebedev_laikov_50<RealType>>( points, weights );
    else if( npts == 74 ) 
      detail::copy_grid<lebedev_laikov_74<RealType>>( points, weights );
    else if( npts == 86 ) 
      detail::copy_grid<lebedev_laikov_86<RealType>>( points, weights );
    else if( npts == 110 ) 
      detail::copy_grid<lebedev_laikov_110<RealType>>( points, weights );
    else if( npts == 146 ) 
      detail::copy_grid<lebedev_laikov_146<RealType>>( points, weights );
    else if( npts == 170 ) 
      detail::copy_grid<lebedev_laikov_170<RealType>>( points, weights );
    else if( npts == 194 ) 
      detail::copy_grid<lebedev_laikov_194<RealType>>( points, weights );
    else if( npts == 230 ) 
      detail::copy_grid<lebedev_laikov_230<RealType>>( points, weights );
    else if( npts == 266 ) 
      detail::copy_grid<lebedev_laikov_266<RealType>>( points, weights );
    else if( npts == 302 ) 
      detail::copy_grid<lebedev_laikov_302<RealType>>( points, weights );
    else if( npts == 350 ) 
      detail::copy_grid<lebedev_laikov_350<RealType>>( points, weights );
    else if( npts == 434 ) 
      detail::copy_grid<lebedev_laikov_434<RealType>>( points, weights );
    else if( npts == 590 ) 
      detail::copy_grid<lebedev_laikov_590<RealType>>( points, weights );
    else if( npts == 770 ) 
      detail::copy_grid<lebedev_laikov_770<RealType>>( points, weights );
    else if( npts == 974 ) 
      detail::copy_grid<lebedev_laikov_974<RealType>>( points, weights );
    else if( npts == 1202 ) 
      detail::copy_grid<lebedev_laikov_1202<RealType>>( points, weights );
    else if( npts == 1454 ) 
      detail::copy_grid<lebedev_laikov_1454<RealType>>( points, weights );
    else if( npts == 1730 ) 
      detail::copy_grid<lebedev_laikov_1730<RealType>>( points, weights );
    else if( npts == 2030 ) 
      detail::copy_grid<lebedev_laikov_2030<RealType>>( points, weights );
    else if( npts == 2354 ) 
      detail::copy_grid<lebedev_laikov_2354<RealType>>( points, weights );
    else if( npts == 2702 ) 
      detail::copy_grid<lebedev_laikov_2702<RealType>>( points, weights );
    else if( npts == 3074 ) 
      detail::copy_grid<lebedev_laikov_3074<RealType>>( points, weights );
    else if( npts == 3470 ) 
      detail::copy_grid<lebedev_laikov_3470<RealType>>( points, weights );
    else if( npts == 3890 ) 
      detail::copy_grid<lebedev_laikov_3890<RealType>>( points, weights );
    else if( npts == 4334 ) 
      detail::copy_grid<lebedev_laikov_4334<RealType>>( points, weights );
    else if( npts == 4802 ) 
      detail::copy_grid<lebedev_laikov_4802<RealType>>( points, weights );
    else if( npts == 5810 ) 
      detail::copy_grid<lebedev_laikov_5810<RealType>>( points, weights );
        
    return std::tuple( points, weights );

  }
};

}
