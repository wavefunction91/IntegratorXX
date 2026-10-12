#include "catch2/catch_all.hpp"
#include <integratorxx/util/legendre.hpp>
#include "test_functions.hpp"

/* Parity check for the hand-rolled associated Legendre polynomial in
 * test_functions.hpp.
 *
 * AssociatedLegendre::evaluate_fallback exists because libc++ does not
 * implement the C++17 mathematical special functions ([sf.cmath]). Where the
 * standard library DOES provide std::assoc_legendre, we check the fallback
 * against it, so that the code path the libc++ platforms depend on is
 * continuously verified rather than exercised only where it is the sole
 * option.
 */

#ifdef IXX_TEST_HAS_STD_ASSOC_LEGENDRE

TEST_CASE( "Associated Legendre matches std::assoc_legendre", "[assoc-legendre]" ) {

  constexpr int    max_l   = 12;
  constexpr int    npts    = 201;
  // Both sides carry rounding error from their own recurrences; measured
  // against an exact rational reference the fallback alone reaches ~1.2e-13
  // relative near l=8. 1e-11 leaves room for the standard library's error and
  // for compiler-to-compiler variation while staying orders of magnitude
  // tighter than any defect would be.
  constexpr double rel_tol = 1e-11;
  constexpr double abs_tol = 1e-13;

  for( int l = 0; l <= max_l; ++l )
  for( int m = 0; m <= l;     ++m ) {

    for( int i = 0; i < npts; ++i ) {
      const double x   = -1.0 + 2.0 * i / (npts - 1.0);
      const double ref = std::assoc_legendre(l, m, x);
      const double val = AssociatedLegendre::evaluate_fallback(l, m, x);

      const std::string msg = "P_l^m (L,M,X) = (" + std::to_string(l) + "," +
        std::to_string(m) + "," + std::to_string(x) + ")";

      // Relative agreement, with an absolute floor so values near a zero of the
      // polynomial are not held to an unattainable relative tolerance
      INFO( msg );
      REQUIRE_THAT( val, Catch::Matchers::WithinRel(ref, rel_tol) ||
                         Catch::Matchers::WithinAbs(ref, abs_tol) );
    }

  }

}

#endif

TEST_CASE( "Associated Legendre special values", "[assoc-legendre]" ) {

  // P_l^0 == P_l, checked against the recurrence in the library proper
  for( int l = 0; l <= 12; ++l )
  for( int i = 0; i < 51;  ++i ) {
    const double x = -1.0 + 2.0 * i / 50.0;
    double p_n;
    std::tie(p_n, std::ignore, std::ignore) = IntegratorXX::eval_Pn(x, l);
    INFO( "P_l^0 vs P_l, (L,X) = (" << l << "," << x << ")" );
    REQUIRE_THAT( AssociatedLegendre::evaluate_fallback(l, 0, x),
                  Catch::Matchers::WithinRel(p_n, 1e-11) ||
                  Catch::Matchers::WithinAbs(p_n, 1e-13) );
  }

  // P_l^m(+/-1) == 0 for m > 0
  for( int l = 1; l <= 12; ++l )
  for( int m = 1; m <= l;  ++m ) {
    REQUIRE( AssociatedLegendre::evaluate_fallback(l, m,  1.0) == 0.0 );
    REQUIRE( AssociatedLegendre::evaluate_fallback(l, m, -1.0) == 0.0 );
  }

  // Closed forms for the low orders
  const double x = 0.375;
  const double s = std::sqrt(1.0 - x*x);
  REQUIRE_THAT( AssociatedLegendre::evaluate_fallback(1, 1, x),
                Catch::Matchers::WithinRel(s, 1e-14) );
  REQUIRE_THAT( AssociatedLegendre::evaluate_fallback(2, 1, x),
                Catch::Matchers::WithinRel(3.0 * x * s, 1e-14) );
  REQUIRE_THAT( AssociatedLegendre::evaluate_fallback(2, 2, x),
                Catch::Matchers::WithinRel(3.0 * (1.0 - x*x), 1e-14) );

}
