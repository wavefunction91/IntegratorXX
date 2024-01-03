#pragma once

#include <integratorxx/quadrature.hpp>
#include <integratorxx/quadratures/s2/lebedev_laikov/lebedev_laikov_grids.hpp>
#include <integratorxx/util/copy_grid.hpp>
#include <vector>

namespace IntegratorXX {


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

    // Pretabulated weights are missing 4 pi
    for(auto i=0; i < npts; i++)
      weights[i] *= 4.0*M_PI;
        
    return std::make_tuple( points, weights );

  }

  inline static int64_t npts_by_algebraic_order(int64_t order) {
  
    switch(order) {
      case 3  : return 6    ; 
      case 5  : return 14   ; 
      case 7  : return 26   ; 
      case 9  : return 38   ; 
      case 11 : return 50   ; 
      case 13 : return 74   ; 
      case 15 : return 86   ; 
      case 17 : return 110  ; 
      case 19 : return 146  ; 
      case 21 : return 170  ; 
      case 23 : return 194  ; 
      case 25 : return 230  ; 
      case 27 : return 266  ; 
      case 29 : return 302  ; 
      case 31 : return 350  ; 
      case 35 : return 434  ; 
      case 41 : return 590  ; 
      case 47 : return 770  ; 
      case 53 : return 974  ; 
      case 59 : return 1202 ; 
      case 65 : return 1454 ; 
      case 71 : return 1730 ; 
      case 77 : return 2030 ; 
      case 83 : return 2354 ; 
      case 89 : return 2702 ; 
      case 95 : return 3074 ; 
      case 101: return 3470 ; 
      case 107: return 3890 ; 
      case 113: return 4334 ; 
      case 119: return 4802 ; 
      case 131: return 5810 ; 
      default:  return -1;
    }

  }

  inline static int64_t algebraic_order_by_npts(int64_t npts) {
  
    switch(npts) {
      case 6   : return 3   ; 
      case 14  : return 5   ; 
      case 26  : return 7   ; 
      case 38  : return 9   ; 
      case 50  : return 11  ; 
      case 74  : return 13  ; 
      case 86  : return 15  ; 
      case 110 : return 17  ; 
      case 146 : return 19  ; 
      case 170 : return 21  ; 
      case 194 : return 23  ; 
      case 230 : return 25  ; 
      case 266 : return 27  ; 
      case 302 : return 29  ; 
      case 350 : return 31  ; 
      case 434 : return 35  ; 
      case 590 : return 41  ; 
      case 770 : return 47  ; 
      case 974 : return 53  ; 
      case 1202: return 59  ; 
      case 1454: return 65  ; 
      case 1730: return 71  ; 
      case 2030: return 77  ; 
      case 2354: return 83  ; 
      case 2702: return 89  ; 
      case 3074: return 95  ; 
      case 3470: return 101 ; 
      case 3890: return 107 ; 
      case 4334: return 113 ; 
      case 4802: return 119 ; 
      case 5810: return 131 ; 
      default:  return -1;
    }

  }


  inline static int64_t next_algebraic_order(int64_t order) {
  
     if( order <= 3  )      return 3  ; 
     else if( order <= 5  ) return 5  ; 
     else if( order <= 7  ) return 7  ; 
     else if( order <= 9  ) return 9  ; 
     else if( order <= 11 ) return 11 ; 
     else if( order <= 13 ) return 13 ; 
     else if( order <= 15 ) return 15 ; 
     else if( order <= 17 ) return 17 ; 
     else if( order <= 19 ) return 19 ; 
     else if( order <= 21 ) return 21 ; 
     else if( order <= 23 ) return 23 ; 
     else if( order <= 25 ) return 25 ; 
     else if( order <= 27 ) return 27 ; 
     else if( order <= 29 ) return 29 ; 
     else if( order <= 31 ) return 31 ; 
     else if( order <= 35 ) return 35 ; 
     else if( order <= 41 ) return 41 ; 
     else if( order <= 47 ) return 47 ; 
     else if( order <= 53 ) return 53 ; 
     else if( order <= 59 ) return 59 ; 
     else if( order <= 65 ) return 65 ; 
     else if( order <= 71 ) return 71 ; 
     else if( order <= 77 ) return 77 ; 
     else if( order <= 83 ) return 83 ; 
     else if( order <= 89 ) return 89 ; 
     else if( order <= 95 ) return 95 ; 
     else if( order <= 101) return 101; 
     else if( order <= 107) return 107; 
     else if( order <= 113) return 113; 
     else if( order <= 119) return 119; 
     else                   return 131; 

  }

};


namespace detail {

template <typename QuadType>
static constexpr bool is_lebedev_laikov_v = std::is_same_v<
  QuadType, 
  LebedevLaikov<typename QuadType::weight_type>
>;

}

}

