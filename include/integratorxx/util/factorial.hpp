#pragma once

namespace IntegratorXX {
constexpr size_t factorial(size_t n) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
}
