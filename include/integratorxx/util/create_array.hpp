#pragma once

namespace IntegratorXX {

namespace detail {
template <typename T, std::size_t... Is>
constexpr std::array<T, sizeof...(Is)> create_array(
    T value, std::index_sequence<Is...>) {
  // cast Is to void to remove the warning: unused value
  return {{(static_cast<void>(Is), value)...}};
}

template <std::size_t N, typename T>
constexpr std::array<T, N> create_array(const T& value) {
  return detail::create_array(value, std::make_index_sequence<N>());
}
}  // namespace detail
}  // namespace IntegratorXX
