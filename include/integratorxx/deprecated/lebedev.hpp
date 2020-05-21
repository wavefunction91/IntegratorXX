#ifndef __INCLUDED_INTEGRATORXX_QUADRATURE_LEBEDEV_HPP__
#define __INCLUDED_INTEGRATORXX_QUADRATURE_LEBEDEV_HPP__

#include "lebedev/lebedev_grids.hpp"

namespace IntegratorXX::deprecated {

  /**
   *  \brief Lebedev-Laikov quadrature factory
   */ 
  template <
    typename PointType, 
    typename wght_t,
    template<typename...> class ContiguousContainer
  >
  class GenerateQuadrature< Lebedev<PointType,wght_t,ContiguousContainer> > {

    using point_container  = typename Lebedev<PointType,wght_t,ContiguousContainer>::point_container;
    using weight_container = typename Lebedev<PointType,wght_t,ContiguousContainer>::weight_container;

    /**
     *  \brief Generate the Lebedev-Laikov quadrature rule of a specific order (impl)
     *
     *  \param[in] nPts Number of quadrature points
     *  \returns [points,weights] tuple of quadrature points and weights
     *
     */ 
    static auto generate_impl( const size_t nPts );

  public:

    /**
     *  \brief Generate the Lebedev-Laikov quadrature rule of a specific order (interface)
     */ 
    template <typename... Args>
    inline static auto generate(Args&&... args){
      return generate_impl( std::forward<Args>(args)... );
    }

  };

  // Implementation of Lebedev-Laikov quadrature factory
  template <
    typename PointType, 
    typename wght_t,
    template<typename...> class ContiguousContainer
  >
  auto GenerateQuadrature< Lebedev<PointType,wght_t,ContiguousContainer> >::
    generate_impl(const size_t nPts) {

    point_container   pts(nPts);
    weight_container  wghts(nPts);

    auto copy_from_iterable = []( const auto& from, auto& to ) {
      std::copy( from.begin(), from.end(), to.begin() );
    };

    if( nPts == 6 ) {
      copy_from_iterable( LebedevGrids::lebedev_6<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_6<PointType>::weights, wghts );
    } else if( nPts == 14 ) {
      copy_from_iterable( LebedevGrids::lebedev_14<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_14<PointType>::weights, wghts );
    } else if( nPts == 26 ) {
      copy_from_iterable( LebedevGrids::lebedev_26<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_26<PointType>::weights, wghts );
    } else if( nPts == 38 ) {
      copy_from_iterable( LebedevGrids::lebedev_38<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_38<PointType>::weights, wghts );
    } else if( nPts == 50 ) {
      copy_from_iterable( LebedevGrids::lebedev_50<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_50<PointType>::weights, wghts );
    } else if( nPts == 74 ) {
      copy_from_iterable( LebedevGrids::lebedev_74<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_74<PointType>::weights, wghts );
    } else if( nPts == 86 ) {
      copy_from_iterable( LebedevGrids::lebedev_86<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_86<PointType>::weights, wghts );
    } else if( nPts == 110 ) {
      copy_from_iterable( LebedevGrids::lebedev_110<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_110<PointType>::weights, wghts );
    } else if( nPts == 146 ) {
      copy_from_iterable( LebedevGrids::lebedev_146<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_146<PointType>::weights, wghts );
    } else if( nPts == 14 ) {
      copy_from_iterable( LebedevGrids::lebedev_14<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_14<PointType>::weights, wghts );
    } else if( nPts == 170 ) {
      copy_from_iterable( LebedevGrids::lebedev_170<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_170<PointType>::weights, wghts );
    } else if( nPts == 194 ) {
      copy_from_iterable( LebedevGrids::lebedev_194<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_194<PointType>::weights, wghts );
    } else if( nPts == 230 ) {
      copy_from_iterable( LebedevGrids::lebedev_230<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_230<PointType>::weights, wghts );
    } else if( nPts == 266 ) {
      copy_from_iterable( LebedevGrids::lebedev_266<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_266<PointType>::weights, wghts );
    } else if( nPts == 302 ) {
      copy_from_iterable( LebedevGrids::lebedev_302<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_302<PointType>::weights, wghts );
    } else if( nPts == 350 ) {
      copy_from_iterable( LebedevGrids::lebedev_350<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_350<PointType>::weights, wghts );
    } else if( nPts == 434 ) {
      copy_from_iterable( LebedevGrids::lebedev_434<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_434<PointType>::weights, wghts );
    } else if( nPts == 590 ) {
      copy_from_iterable( LebedevGrids::lebedev_590<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_590<PointType>::weights, wghts );
    } else if( nPts == 770 ) {
      copy_from_iterable( LebedevGrids::lebedev_770<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_770<PointType>::weights, wghts );
    } else if( nPts == 974 ) {
      copy_from_iterable( LebedevGrids::lebedev_974<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_974<PointType>::weights, wghts );
    } else if( nPts == 1202 ) {
      copy_from_iterable( LebedevGrids::lebedev_1202<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_1202<PointType>::weights, wghts );
    } else if( nPts == 1454 ) {
      copy_from_iterable( LebedevGrids::lebedev_1454<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_1454<PointType>::weights, wghts );
    } else if( nPts == 1730 ) {
      copy_from_iterable( LebedevGrids::lebedev_1730<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_1730<PointType>::weights, wghts );
    } else if( nPts == 2030 ) {
      copy_from_iterable( LebedevGrids::lebedev_2030<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_2030<PointType>::weights, wghts );
    } else if( nPts == 2354 ) {
      copy_from_iterable( LebedevGrids::lebedev_2354<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_2354<PointType>::weights, wghts );
    } else if( nPts == 2702 ) {
      copy_from_iterable( LebedevGrids::lebedev_2702<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_2702<PointType>::weights, wghts );
    } else if( nPts == 3074 ) {
      copy_from_iterable( LebedevGrids::lebedev_3074<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_3074<PointType>::weights, wghts );
    } else if( nPts == 3470 ) {
      copy_from_iterable( LebedevGrids::lebedev_3470<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_3470<PointType>::weights, wghts );
    } else if( nPts == 3890 ) {
      copy_from_iterable( LebedevGrids::lebedev_3890<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_3890<PointType>::weights, wghts );
    } else if( nPts == 4334 ) {
      copy_from_iterable( LebedevGrids::lebedev_4334<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_4334<PointType>::weights, wghts );
    } else if( nPts == 4802 ) {
      copy_from_iterable( LebedevGrids::lebedev_4802<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_4802<PointType>::weights, wghts );
    } else if( nPts == 5294 ) {
      copy_from_iterable( LebedevGrids::lebedev_5294<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_5294<PointType>::weights, wghts );
    } else if( nPts == 5810 ) {
      copy_from_iterable( LebedevGrids::lebedev_5810<PointType>::points,  pts   );
      copy_from_iterable( LebedevGrids::lebedev_5810<PointType>::weights, wghts );
    }

    return std::tuple( std::move(pts), std::move(wghts) );

  }

}

#endif
