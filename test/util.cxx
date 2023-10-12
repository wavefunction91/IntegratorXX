#include "catch2/catch_all.hpp"
#include <integratorxx/util/factorial.hpp>
#include <integratorxx/util/gamma.hpp>
#include <integratorxx/util/pow.hpp>
#include <integratorxx/util/pow_2.hpp>
#include <iostream>

uint64_t dumb_factorial(uint64_t n) {
  if( n == 0 or n == 1 ) return 1;
  else return n * dumb_factorial(n-1);
}

uint64_t dumb_double_factorial(uint64_t n) {
  if( n == 0 or n == 1 ) return 1;
  else return n * dumb_double_factorial(n-2);
}


TEST_CASE("Factorial", "[util]") {
  using namespace IntegratorXX;

  REQUIRE(std::is_same_v<
    decltype(factorial(std::declval<int>())), int>);
  REQUIRE(std::is_same_v<
    decltype(factorial(std::declval<unsigned>())), unsigned>);

  for(uint64_t i = 0; i <= 20; ++i) {
    REQUIRE(factorial(i) == dumb_factorial(i));
  }

  REQUIRE(double_factorial(-1) == 1);
  REQUIRE(double_factorial(-3) == -1);
  for(uint64_t i = 0; i <= 32; ++i) {
    REQUIRE(double_factorial(i) == dumb_double_factorial(i));
  }
}

TEST_CASE("Pow", "[util]") {

  using namespace IntegratorXX;
  using namespace Catch::Matchers;
  for(int i = 0; i < 10; ++i) {
    REQUIRE_THAT(integer_pow(2.0, i), WithinAbs(std::pow(2.0,i), 1e-16));
    REQUIRE_THAT(half_integer_pow(2.0, i), WithinAbs(std::pow(2.0,i/2.0), 1e-16));
  }

}

TEST_CASE("Pow 2", "[util]") {
  using namespace IntegratorXX;

  REQUIRE(std::is_same_v<
    decltype(integer_pow_2(std::declval<int>())), int>);
  REQUIRE(std::is_same_v<
    decltype(integer_pow_2(std::declval<unsigned>())), unsigned>);

  REQUIRE(std::is_same_v<
    decltype(half_integer_pow_2<float>(std::declval<int>())), float>);
  REQUIRE(std::is_same_v<
    decltype(half_integer_pow_2<float>(std::declval<unsigned>())), float>);
  REQUIRE(std::is_same_v<
    decltype(half_integer_pow_2<double>(std::declval<int>())), double>);
  REQUIRE(std::is_same_v<
    decltype(half_integer_pow_2<double>(std::declval<unsigned>())), double>);

  for(uint64_t i = 0; i < 64; ++i) {
    REQUIRE(integer_pow_2(i) == (1ul << i));
  }

  // 126 to avoid INT64_T overflow
  for(int64_t i = 0; i < 126; i += 2) {
    const auto v = half_integer_pow_2<double>(i);
    REQUIRE( v == (1ul << i/2));
    REQUIRE(half_integer_pow_2<double>(-i) == 1.0/v);
    REQUIRE(half_integer_pow_2<double>(i+1) == M_SQRT2*v);
    REQUIRE(half_integer_pow_2<double>(-(i+1)) == 1.0/(M_SQRT2*v));
  }

}

TEST_CASE("Gamma", "[util]") {
  using namespace IntegratorXX;

  REQUIRE(std::is_same_v<
    decltype(integer_tgamma<float>(std::declval<int>())), float>);
  REQUIRE(std::is_same_v<
    decltype(integer_tgamma<float>(std::declval<unsigned>())), float>);
  REQUIRE(std::is_same_v<
    decltype(integer_tgamma<double>(std::declval<int>())), double>);
  REQUIRE(std::is_same_v<
    decltype(integer_tgamma<double>(std::declval<unsigned>())), double>);
  
  REQUIRE(std::isnan(integer_tgamma<float>(-1)));
  REQUIRE(std::isnan(integer_tgamma<float>(0)));
  for(uint64_t i = 1; i <= 21; ++i) {
    REQUIRE(integer_tgamma<double>(i) == factorial(i-1));
  }

  REQUIRE(std::isnan(half_integer_tgamma<float>(0)));
  for(int64_t i = 1; i <= 34; ++i) {
    const auto v = half_integer_tgamma<long double>(i);
    const auto r = std::tgamma(i / 2.0l);
    REQUIRE_THAT(v, Catch::Matchers::WithinAbs(r, 1e-20));  
  }

}

