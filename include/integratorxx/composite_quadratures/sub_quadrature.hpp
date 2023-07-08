#pragma once

#include <cassert>
#include <integratorxx/quadrature.hpp>

namespace IntegratorXX {

template <typename DerivedQuad>
class SubQuadrature : public Quadrature<SubQuadrature<DerivedQuad>> {
 private:
  using self_type = SubQuadrature<DerivedQuad>;
  using base_type = Quadrature<self_type>;
  using traits = quadrature_traits<self_type>;

 public:
  using point_type = typename traits::point_type;
  using weight_type = typename traits::weight_type;
  using point_container = typename traits::point_container;
  using weight_container = typename traits::weight_container;
  using range_idx_type = typename traits::range_idx_type;

  SubQuadrature(range_idx_type idx_range, const Quadrature<DerivedQuad>& quad)
      : base_type(idx_range, quad) {}
};

template <typename DerivedQuad>
struct quadrature_traits<SubQuadrature<DerivedQuad>> {
  using derived_traits = quadrature_traits<DerivedQuad>;

  using range_idx_type = std::pair<size_t, size_t>;

  using point_type = typename derived_traits::point_type;
  using weight_type = typename derived_traits::weight_type;
  using point_container = typename derived_traits::point_container;
  using weight_container = typename derived_traits::weight_container;

  inline static std::tuple<point_container, weight_container> generate(
      range_idx_type idx_range, const Quadrature<DerivedQuad>& quad) {
    const auto [begin_idx, end_idx] = idx_range;
    assert(begin_idx >= 0);
    assert(begin_idx <= quad.npts());
    assert(end_idx >= 0);
    assert(end_idx <= quad.npts());
    assert(begin_idx < end_idx);

    point_container points(quad.points().begin() + begin_idx,
                           quad.points().begin() + end_idx);
    point_container weights(quad.weights().begin() + begin_idx,
                            quad.weights().begin() + end_idx);

    return std::make_tuple(points, weights);
  }
};
}  // namespace IntegratorXX
