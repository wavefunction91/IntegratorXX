#pragma once

#include <integratorxx/quadratures/primitive/uniform.hpp>
#include <integratorxx/quadratures/radial/radial_transform.hpp>
#include <vector>
#include <cmath>

namespace IntegratorXX {


class MuraKnowlesRadialTraits : public RadialTraits {

  size_t npts_; ///< Number of grid points
  double R_; ///< Radial scaling factor

public:

  MuraKnowlesRadialTraits(size_t npts, double R = 1.0) : npts_(npts), R_(R) { }

  size_t npts() const noexcept { return npts_; }

  std::unique_ptr<RadialTraits> clone() const {
    return std::make_unique<MuraKnowlesRadialTraits>(*this);
  }

  bool compare(const RadialTraits& other) const noexcept {
    auto ptr = dynamic_cast<const MuraKnowlesRadialTraits*>(&other);
    return ptr ? *this == *ptr : false;
  }

  bool operator==(const MuraKnowlesRadialTraits& other) const noexcept {
    return npts_ == other.npts_ && R_ == other.R_;
  }

  template <typename PointType>
  inline auto radial_transform(PointType x) const noexcept {
    return -R_ * std::log(1.0 - x*x*x);
  }

  template <typename PointType>
  inline auto radial_jacobian(PointType x) const noexcept {
    const auto x2 = x*x;
    return R_ * 3.0 * x2 / (1.0 - x2 * x);
  }

}; 

template <typename PointType, typename WeightType>
using MuraKnowles = RadialTransformQuadrature<
  UniformTrapezoid<PointType,WeightType>,
  MuraKnowlesRadialTraits
>;


namespace detail {

template <typename QuadType>
static constexpr bool is_mk_v = std::is_same_v<
  QuadType, 
  MuraKnowles<typename QuadType::point_type, typename QuadType::weight_type>
>;

}

}
