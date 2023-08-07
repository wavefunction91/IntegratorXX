#pragma once

#include <tuple>
#include <cstdint>
#include <cstddef>
#include <cassert>
#include <integratorxx/type_traits.hpp>

namespace IntegratorXX {

// Forward decl of quadrature traits
template <typename Quadrature, typename = std::void_t<>>
struct quadrature_traits;


/**
 *  \brief Base class for Quadrature.
 */
template <typename PointContainer, typename WeightContainer>
class QuadratureBase {

public:

  using point_container  = PointContainer;
  using weight_container = WeightContainer;

  using point_type  = typename point_container::value_type;
  using weight_type = typename weight_container::value_type;

protected:

  point_container  points_;
  weight_container weights_;

  void sanity_check() {
#ifndef INTEGRATORXX_DISABLE_SANITY_CHECK
    assert(points_.size() == weights_.size());
#endif
  }

  QuadratureBase( const point_container& p, const weight_container& w ):
    points_(p), weights_(w) { sanity_check(); }
  QuadratureBase( point_container&& p, weight_container&& w ):
    points_(std::move(p)), weights_(std::move(w)) { sanity_check(); }

public:


  /**
   *  @brief Getter for internal storage of quadrature points
   *  
   *  Const variant.
   *
   *  @returns A const reference to internal point storage
   */
  const auto& points()  const { return points_; }

  /**
   *  @brief Getter for internal storage of quadrature points
   *  
   *  Non-const variant.
   *
   *  @returns A non-const reference to internal point storage
   */
  auto&       points()        { return points_; }

  /**
   *  @brief Getter for internal storage of quadrature weights
   *  
   *  Const variant.
   *
   *  @returns A const reference to internal weight storage
   */
  const auto& weights() const { return weights_; }

  /**
   *  @brief Getter for internal storage of quadrature weights
   *  
   *  Non-const variant.
   *
   *  @returns A non-const reference to internal weight storage
   */
  auto&       weights()       { return weights_; }



  /**
   *  @brief Obtain a pointer for underlying quadrature point storage
   *  
   *  Const variant.
   *
   *  @returns A const pointer to quadrature point storage
   */
  const auto* points_data()  const { return points_.data();  }

  /**
   *  @brief Obtain a pointer for underlying quadrature point storage
   *  
   *  Non-const variant.
   *
   *  @returns A non-const pointer to quadrature point storage
   */
  auto*       points_data()        { return points_.data();  }

  /**
   *  @brief Obtain a pointer for underlying quadrature weight storage
   *  
   *  Const variant.
   *
   *  @returns A const pointer to quadrature weight storage
   */
  const auto* weights_data() const { return weights_.data(); }

  /**
   *  @brief Obtain a pointer for underlying quadrature weight storage
   *  
   *  Non-const variant.
   *
   *  @returns A non-const pointer to quadrature weight storage
   */
  auto*       weights_data()       { return weights_.data(); }








  /**
   *  @brief Obtain number of points in implemented quadrature
   *  @returns Number of points in implemented quadrature
   */
  size_t npts() const { return points_.size(); }


  /**
   *  @brief Obtain a copy of a particular quadrature point
   *
   *  Performs a bounds check, should not use in performance
   *  critial code
   *
   *  @param[in] i Point index of desired quadrature point
   *  @returns   Quadrature point at index i
   */
  point_type  points(size_t i)  const { return points_.at(i);  }

  /**
   *  @brief Obtain a copy of a particular quadrature weight
   *
   *  Performs a bounds check, should not use in performance
   *  critial code
   *
   *  @param[in] i Point index of desired quadrature weight
   *  @returns   Quadrature weight at index i
   */
  weight_type weights(size_t i) const { return weights_.at(i); }

};

/**
 *  @brief Quadrature base class: Provides minimal storage and user interface 
 *  for quadrature manipulation
 *
 *  @tparam DerivedQuadrature Implemented quadrature type, must admit template
 *  specialization for quadrature_traits
 */
template <typename Derived>
class Quadrature : public
  QuadratureBase< 
    typename quadrature_traits<Derived>::point_container,
    typename quadrature_traits<Derived>::weight_container
  > {

private:

  using derived_traits = quadrature_traits<Derived>;

public:

  using point_type  = typename derived_traits::point_type;
  using weight_type = typename derived_traits::weight_type;

  using point_container  = typename derived_traits::point_container;
  using weight_container = typename derived_traits::weight_container;

private:

  using base_type = QuadratureBase< point_container, weight_container >;

  template <typename... Args>
  inline static constexpr auto 
    generate( Args&&... args ) {
    return derived_traits::generate( std::forward<Args>(args)... );
  }


protected:

  using quadrature_return_type =
    std::tuple<point_container,weight_container>;

  Quadrature( const quadrature_return_type& q ) :
    base_type( std::get<0>(q), std::get<1>(q) ) { }

  Quadrature( quadrature_return_type&& q ) :
    base_type( std::move(std::get<0>(q)), std::move(std::get<1>(q)) ) { }

public:

  /**
   *  @brief Construct a Quadrature object from options.
   *
   *  Generates the implemented quadrature and populates internal
   *  storage.
   *
   *  @tparam Args Parameter pack to be forwarded to quadrature generator
   *
   *  @param[in] args Parameter pack to be forwarded to quadrature generator
   *
   */
  template <typename... Args,
    typename = std::enable_if_t<
      // Disable generic construction for copy/move to and from base type
      detail::all_are_not< Quadrature<Derived>, std::decay_t<Args>... >::value and
      detail::all_are_not<            Derived , std::decay_t<Args>... >::value
    >
  >
  Quadrature( Args&&... args ) :
    Quadrature( generate( std::forward<Args>(args)... ) ) { } 


  // Default copy and move ctors

  Quadrature( const Quadrature<Derived>& )     = default;
  Quadrature( Quadrature<Derived>&& ) noexcept = default;

  Quadrature( const Derived& d ) :
    Quadrature( dynamic_cast<const Quadrature&>(d) ) { }
  Quadrature( Derived&& d ) noexcept :
    Quadrature( dynamic_cast<Quadrature&&>(std::move(d)) ) { }

};


}
