#pragma once

#include <tuple>

namespace IntegratorXX {

// Forward decl of quadrature traits
template <typename Quadrature>
struct quadrature_traits;


/**
 *  @brief Quadrature base class: Provides minimal storage and user interface 
 *  for quadrature manipulation
 *
 *  @tparam DerivedQuadrature Implemented quadrature type, must admit template
 *  specialization for quadrature_traits
 */
template <typename DerivedQuadrature>
class Quadrature {

private:

  using derived_traits = quadrature_traits<DerivedQuadrature>;

public:

  using point_type  = typename derived_traits::point_type;
  using weight_type = typename derived_traits::weight_type;

  using point_container  = typename derived_traits::point_container;
  using weight_container = typename derived_traits::weight_container;

private:


  template <typename... Args>
  inline static constexpr auto 
    generate( Args&&... args ) {
    return derived_traits::generate( std::forward<Args>(args)... );
  }


protected:

  point_container  points_;
  weight_container weights_;

  Quadrature( const point_container& p, const weight_container& w ):
    points_(p), weights_(w) { }
  Quadrature( point_container&& p, weight_container&& w ):
    points_(std::move(p)), weights_(std::move(w)) { }

  using quadrature_return_type =
    std::tuple<point_container,weight_container>;

  Quadrature( const quadrature_return_type& q ) :
    Quadrature( std::get<0>(q), std::get<1>(q) ) { }

  Quadrature( quadrature_return_type&& q ) :
    Quadrature( std::move(std::get<0>(q)), std::move(std::get<1>(q)) ) { }

public:

  template <typename... Args>
  Quadrature( Args&&... args ) :
    Quadrature( generate( std::forward<Args>(args)... ) ) { } 

  Quadrature( const Quadrature& ) = default;
  Quadrature( Quadrature&& ) noexcept = default;

  const auto& points()  const { return points_; }
  auto&       points()        { return points_; }
  const auto& weights() const { return weights_; }
  auto&       weights()       { return weights_; }

  size_t npts() const { return points_.size(); }

};


}
