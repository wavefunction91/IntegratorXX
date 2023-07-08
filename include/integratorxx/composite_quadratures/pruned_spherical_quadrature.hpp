#pragma once

#include <integratorxx/composite_quadratures/spherical_quadrature.hpp>
#include <integratorxx/composite_quadratures/sub_quadrature.hpp>
#include <integratorxx/type_traits.hpp>
#include <integratorxx/types.hpp>
#include <utility>

namespace IntegratorXX {

template <typename AngularQuad>
class RadialGridPartition {
  std::vector<size_t> partition_idx_;
  std::vector<Quadrature<AngularQuad>> quads_;

  template <typename RadialQuad>
  inline void add_quads(const RadialQuad& rq) {
    finalize(rq);
  }

  template <typename RadialQuad, typename... Args>
  void add_quads(const RadialQuad& rq, size_t idx,
                 const Quadrature<AngularQuad>& q, Args&&... args) {
    add_quad(rq, idx, q);
    add_quads(rq, std::forward<Args>(args)...);
  }

 public:
  template <typename RadialQuad>
  void add_quad(const RadialQuad& rq, size_t idx,
                const Quadrature<AngularQuad>& q) {
    assert(partition_idx_.size()
               ? (idx > partition_idx_.back() && idx < (rq.npts() - 1))
               : (idx == 0));
    (void)(rq);

    partition_idx_.emplace_back(idx);
    quads_.emplace_back(q);
  }

  template <typename RadialQuad>
  void finalize(const RadialQuad& rq) {
    if(partition_idx_.back() != rq.npts())
      partition_idx_.emplace_back(rq.npts());
  }

  using angular_type = AngularQuad;

  RadialGridPartition() = default;

  RadialGridPartition(const RadialGridPartition& other)
      : partition_idx_(other.partition_idx_), quads_(other.quads_){};
  RadialGridPartition(RadialGridPartition&& other) noexcept
      : partition_idx_(std::move(other.partition_idx_)),
        quads_(std::move(other.quads_)){};

  RadialGridPartition& operator=(const RadialGridPartition& other) {
    partition_idx_ = other.partition_idx_;
    quads_ = other.quads_;
    return *this;
  }
  RadialGridPartition& operator=(RadialGridPartition&& other) {
    partition_idx_ = std::move(other.partition_idx_);
    quads_ = std::move(other.quads_);
    return *this;
  }

  template <typename RadialQuad, typename... Args>
  RadialGridPartition(const RadialQuad& rq, size_t idx,
                      const Quadrature<AngularQuad>& q, Args&&... args) {
    add_quads(rq, idx, q, std::forward<Args>(args)...);
  }

  auto& partition_idx() { return partition_idx_; }

  template <typename index_iterator, typename quad_iterator>
  struct rgp_iterator {
    using index_range_type = std::pair<size_t, size_t>;
    using quad_type = Quadrature<AngularQuad>;
    using value_type = std::pair<index_range_type, quad_type>;
    using difference_type = size_t;
    using pointer = value_type*;
    using reference = value_type&;
    using iterator_catagory = std::input_iterator_tag;

    index_iterator idx_it;
    quad_iterator quad_it;

    rgp_iterator() = delete;
    rgp_iterator(index_iterator ii, quad_iterator pi)
        : idx_it(ii), quad_it(pi) {}

    rgp_iterator& operator++() {
      idx_it++;
      quad_it++;
      return *this;
    }
    rgp_iterator operator++(int) {
      iterator retval = *this;
      ++(*this);
      return retval;
    }

    rgp_iterator& operator+(int i) {
      idx_it += i;
      quad_it += i;
      return (*this);
    }

    bool operator==(rgp_iterator other) const {
      return idx_it == other.idx_it and quad_it == other.quad_it;
    }
    bool operator!=(rgp_iterator other) const { return not(*this == other); }

    value_type operator*() {
      return std::make_pair(std::make_pair(*idx_it, *(idx_it + 1)), *quad_it);
    }
  };

  using iterator = rgp_iterator<typename decltype(partition_idx_)::iterator,
                                typename decltype(quads_)::iterator>;
  using const_iterator =
      rgp_iterator<typename decltype(partition_idx_)::const_iterator,
                   typename decltype(quads_)::const_iterator>;

  iterator begin() { return iterator(partition_idx_.begin(), quads_.begin()); }
  iterator end() { return iterator(partition_idx_.end() - 1, quads_.end()); }

  const_iterator cbegin() const {
    return const_iterator(partition_idx_.cbegin(), quads_.cbegin());
  }
  const_iterator cend() const {
    return const_iterator(partition_idx_.cend() - 1, quads_.cend());
  }

  const_iterator begin() const {
    return const_iterator(partition_idx_.cbegin(), quads_.cbegin());
  }
  const_iterator end() const {
    return const_iterator(partition_idx_.cend() - 1, quads_.cend());
  }
};

template <typename RadialQuad, typename AngularQuad>
class PrunedSphericalQuadrature
    : public Quadrature<PrunedSphericalQuadrature<RadialQuad, AngularQuad>>,
      public SphericalQuadratureBase<typename AngularQuad::point_container,
                                     typename AngularQuad::weight_container> {
  using self_type = PrunedSphericalQuadrature<RadialQuad, AngularQuad>;
  using traits = quadrature_traits<self_type>;

  using sph_base =
      SphericalQuadratureBase<typename AngularQuad::point_container,
                              typename AngularQuad::weight_container>;

 public:
  using quad_base_type = Quadrature<self_type>;

  using point_type = typename traits::point_type;
  using weight_type = typename traits::weight_type;
  using point_container = typename traits::point_container;
  using weight_container = typename traits::weight_container;

 protected:
  const point_container& sph_points_adaptor() const override {
    return quad_base_type::points();
  }
  point_container& sph_points_adaptor() override {
    ;
    return quad_base_type::points();
  }
  const weight_container& sph_weights_adaptor() const override {
    ;
    return quad_base_type::weights();
  }
  weight_container& sph_weights_adaptor() override {
    ;
    return quad_base_type::weights();
  }

 public:
  PrunedSphericalQuadrature(const Quadrature<RadialQuad>& rq,
                            const RadialGridPartition<AngularQuad>& rgp,
                            const point_type cen = point_type({0., 0., 0.}))
      : quad_base_type(rq, rgp, cen), sph_base(cen) {}

  PrunedSphericalQuadrature(const PrunedSphericalQuadrature&) = default;
  PrunedSphericalQuadrature(PrunedSphericalQuadrature&&) noexcept = default;

  const auto& points() const { return quad_base_type::points(); }
  auto& points() { return quad_base_type::points(); }
  const auto& weights() const { return quad_base_type::weights(); }
  auto& weights() { return quad_base_type::weights(); }

  size_t npts() const { return quad_base_type::npts(); }

  std::shared_ptr<sph_base> clone() const override {
    return std::make_shared<self_type>(*this);
  }
};

template <typename RadialQuad, typename AngularQuad>
struct quadrature_traits<PrunedSphericalQuadrature<RadialQuad, AngularQuad>> {
  using radial_subquad_type = SubQuadrature<RadialQuad>;
  using prod_type = ProductQuadrature<
      detail::spherical_combine_op<typename RadialQuad::point_type>,
      radial_subquad_type, AngularQuad>;

  using radial_subquad_traits = quadrature_traits<radial_subquad_type>;
  using prod_traits = quadrature_traits<prod_type>;

  using radial_traits = quadrature_traits<RadialQuad>;
  using radial_point_container = typename radial_traits::point_container;

  using point_type = typename prod_traits::point_type;
  using weight_type = typename prod_traits::weight_type;
  using point_container = typename prod_traits::point_container;
  using weight_container = typename prod_traits::weight_container;

  inline static void shift_grid(point_container& points, point_type vector) {
    size_t npts = points.size();
    for(size_t i = 0; i < npts; ++i) {
      auto& p = points[i];
      p[0] += vector[0];
      p[1] += vector[1];
      p[2] += vector[2];
    }
  }

  inline static std::tuple<point_container, weight_container> generate(
      const Quadrature<RadialQuad>& r,
      const RadialGridPartition<AngularQuad>& rgp, const point_type& cen) {
    point_container points;
    weight_container weights;

    // Loop over angular quadratures + associated ranges
    for(const auto&& [r_range, a] : std::as_const(rgp)) {
      // Generate radial subquadrature
      radial_subquad_type r_subquad(r_range, r);

      // Generate product subquadrature
      auto [sub_points, sub_weights] = prod_traits::generate(r_subquad, a);

      // Include Spherical Jacobian
      const auto npr = r_subquad.npts();
      const auto npa = a.npts();

      const auto& rp = r_subquad.points();

      // Include Spherical Jacobian
      for(size_t j = 0; j < npa; ++j)
        for(size_t i = 0; i < npr; ++i) {
          const auto ij = i + j * npr;
          sub_weights[ij] *= 4 * M_PI * rp[i] * rp[i];
        }

      // Insert points/weights into final storage
      points.insert(points.end(), sub_points.begin(), sub_points.end());
      weights.insert(weights.end(), sub_weights.begin(), sub_weights.end());
    }

    // Recenter
    shift_grid(points, cen);

    return std::make_tuple(points, weights);
  }
};

}  // namespace IntegratorXX
