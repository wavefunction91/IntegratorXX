#pragma once

#include <limits>

namespace IntegratorXX {

template <typename PointType, typename WeightType>
constexpr inline auto transform_minus_one_to_one(const PointType lo,
                                                 const PointType up,
                                                 PointType pt, WeightType wgt) {
  assert(lo < up);

  constexpr auto inf = std::numeric_limits<PointType>::infinity();
  const bool up_is_inf = up == inf;
  const bool lo_is_minf = lo == -inf;

  // Both upper and lower bounds are finite: Map (-1,1) -> (lo,up)
  if(not up_is_inf and not lo_is_minf) {
    // Factor jacobian into the weights
    // dx' = ((b-a)/2) * dx
    wgt *= (up - lo) / 2.;

    // Perform the coordinate shift
    // x' = ((b-a) / 2)*x + (a+b)/2
    pt = (up - lo) * pt / 2. + (up + lo) / 2.;

    // Lower bound is finite, upper is infinite: Map (-1,1) -> (lo,inf)
  } else if(not lo_is_minf) {
    // Factor jacobian into the weights
    // dx' = \frac{2}{(1-x)^2} dx
    wgt *= 2.0 / ((1 - pt) * (1 - pt));

    // Perform the coordinate shift
    // x' = a + \frac{(1+x)}{(1-x)}
    pt = lo + (1 + pt) / (1 - pt);

    // Upper bound is finite, low is infinite: Map (-1,1) -> (inf,up)
  } else if(not up_is_inf) {
    // Factor jacobian into the weights
    // dx' = - \frac{2}{(1+x)^2} dx
    wgt *= -2.0 / ((1 + pt) * (1 + pt));

    // Perform the coordinate shift
    // x' = b - \frac{(1-x)}{(1+x)}
    pt = up - (1 - pt) / (1 + pt);

    // Both infinite: Map (-1,1) -> (-inf,inf)
  } else {
    // Factor jacobian into the weights
    // dx' = \frac{1+x^2}{(x^2 - 1)^2} dx
    wgt *= (1 + pt * pt) / ((pt * pt - 1) * (pt * pt - 1));

    // Perform the coordinate shift
    // x' = x / (1+x) / (1-x)
    pt = pt / (1 + pt) / (1 - pt);
  }

  return std::tuple(pt, wgt);
}

}  // namespace IntegratorXX
