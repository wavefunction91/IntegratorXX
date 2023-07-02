#pragma once

namespace IntegratorXX {

namespace detail {
template <class StaticGrid, typename PointContainer, typename WeightsContainer>
void copy_grid(PointContainer &points, WeightsContainer &weights) {

  const auto &static_points = StaticGrid::points;
  const auto &static_weights = StaticGrid::weights;

  std::copy(static_points.begin(), static_points.end(), points.begin());
  std::copy(static_weights.begin(), static_weights.end(), weights.begin());
}
} // namespace detail
} // namespace IntegratorXX
