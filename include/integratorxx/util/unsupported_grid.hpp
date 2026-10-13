#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>

namespace IntegratorXX {
namespace detail {

/**
 *  @brief Report a request for an angular grid size that is not tabulated.
 *
 *  The angular quadratures are only defined for the specific point counts
 *  that integrate spherical harmonics exactly up to a given algebraic
 *  order. Requesting any other size is a programming error rather than a
 *  recoverable condition, but it is diagnosed at runtime because the size
 *  is typically read from user input.
 *
 *  @param[in] family Name of the angular quadrature family
 *  @param[in] npts   Unsupported number of points that was requested
 */
[[noreturn]] inline void throw_unsupported_grid_size(const char* family,
                                                     size_t npts) {
  throw std::runtime_error(
      std::string("IntegratorXX: ") + family + " does not tabulate a " +
      std::to_string(npts) +
      "-point grid. Supported sizes are given by "
      "quadrature_traits<" +
      family + ">::npts_by_algebraic_order(order).");
}

}  // namespace detail
}  // namespace IntegratorXX
