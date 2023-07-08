#include <integratorxx/quadratures/gausslegendre.hpp>

#include "catch2/catch_all.hpp"

TEST_CASE("Constructor", "[quad]") {
  using quad_type = IntegratorXX::GaussLegendre<double, double>;
  using base_type = IntegratorXX::Quadrature<quad_type>;

  SECTION("Basic") { quad_type q(10, -1, 1); }

  SECTION("Copy To Base") {
    quad_type q(10, -1, 1);
    base_type b(q);

    CHECK(q.points() == b.points());
    CHECK(q.weights() == b.weights());
  }

  SECTION("Move to Base") {
    quad_type q(10, -1, 1);
    std::vector<double> pts = q.points();
    std::vector<double> wgt = q.weights();

    base_type b(std::move(q));

    CHECK(q.points().size() == 0);
    CHECK(q.weights().size() == 0);
    CHECK(b.points() == pts);
    CHECK(b.weights() == wgt);
  }
}
