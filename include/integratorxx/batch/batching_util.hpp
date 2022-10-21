#pragma once
#include <integratorxx/type_traits.hpp>


namespace IntegratorXX::detail {

template <typename T>
using pw_type = std::pair< cartesian_pt_t<T>, T >;

template <typename T>
using pw_container = typename std::vector< pw_type<T> >;

template <typename T>
using pw_iterator  = typename pw_container<T>::iterator;





template <typename T>
bool point_in_box( 
  const cartesian_pt_t<T>& lo, 
  const cartesian_pt_t<T>& up,
  const cartesian_pt_t<T>& pt
) {

  if( pt[0] > up[0] or pt[0] < lo[0] )      return false;
  else if( pt[1] > up[1] or pt[1] < lo[1] ) return false;
  else if( pt[2] > up[2] or pt[2] < lo[2] ) return false;
  else                                      return true;

}


template <size_t I, typename T>
inline constexpr T get_coordinate( const cartesian_pt_t<T>& pt ) {
  return pt[I];
}

template <size_t I, typename T, typename BinaryPred>
inline constexpr bool comp_coordinate( 
  const cartesian_pt_t<T>& a, 
  const cartesian_pt_t<T>& b, 
  const BinaryPred& pred 
) {

  return pred( get_coordinate<I>(a), get_coordinate<I>(b) );

}


template <typename PWIterator>
auto get_box_bounds(
  PWIterator  pw_batch_begin,
  PWIterator  pw_batch_end
) {

  using traits = std::iterator_traits<PWIterator>;
  using cart_type = typename traits::value_type::first_type::value_type;
  using less = std::less<cart_type>;

  auto comp_x = [](const auto& a, const auto& b) {
    return comp_coordinate<0>( a.first, b.first, less());
  };
  auto comp_y = [](const auto& a, const auto& b) {
    return comp_coordinate<1>( a.first, b.first, less());
  };
  auto comp_z = [](const auto& a, const auto& b) {
    return comp_coordinate<2>( a.first, b.first, less());
  };

  auto [min_x_it, max_x_it] = 
    std::minmax_element( pw_batch_begin, pw_batch_end, comp_x );
  auto [min_y_it, max_y_it] = 
    std::minmax_element( pw_batch_begin, pw_batch_end, comp_y );
  auto [min_z_it, max_z_it] = 
    std::minmax_element( pw_batch_begin, pw_batch_end, comp_z );

  const auto min_x = min_x_it->first[0];
  const auto max_x = max_x_it->first[0];
  const auto min_y = min_y_it->first[1];
  const auto max_y = max_y_it->first[1];
  const auto min_z = min_z_it->first[2];
  const auto max_z = max_z_it->first[2];

  return std::make_tuple(
    cartesian_pt_t<cart_type>{ min_x, min_y, min_z },
    cartesian_pt_t<cart_type>{ max_x, max_y, max_z }
  );
}

template <typename PIterator>
auto get_box_bounds_points(
  PIterator  p_batch_begin,
  PIterator  p_batch_end
) {

  using traits = std::iterator_traits<PIterator>;
  using cart_type = typename traits::value_type::value_type;
  using less = std::less<cart_type>;

  auto comp_x = [](const auto& a, const auto& b) {
    return comp_coordinate<0>( a, b, less());
  };
  auto comp_y = [](const auto& a, const auto& b) {
    return comp_coordinate<1>( a, b, less());
  };
  auto comp_z = [](const auto& a, const auto& b) {
    return comp_coordinate<2>( a, b, less());
  };

  auto [min_x_it, max_x_it] = 
    std::minmax_element( p_batch_begin, p_batch_end, comp_x );
  auto [min_y_it, max_y_it] = 
    std::minmax_element( p_batch_begin, p_batch_end, comp_y );
  auto [min_z_it, max_z_it] = 
    std::minmax_element( p_batch_begin, p_batch_end, comp_z );

  const auto min_x = (*min_x_it)[0];
  const auto max_x = (*max_x_it)[0];
  const auto min_y = (*min_y_it)[1];
  const auto max_y = (*max_y_it)[1];
  const auto min_z = (*min_z_it)[2];
  const auto max_z = (*max_z_it)[2];

  return std::make_tuple(
    cartesian_pt_t<cart_type>{ min_x, min_y, min_z },
    cartesian_pt_t<cart_type>{ max_x, max_y, max_z }
  );
}


}
