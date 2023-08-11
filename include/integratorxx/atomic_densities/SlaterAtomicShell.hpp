#pragma once
#include <cassert>
#include <cmath>
#include <vector>

class SlaterTypeAtomicShell {
  using int_container = std::vector<int>;
  using real_container = std::vector<double>;
  using coefficient_container = std::vector<real_container>;

  /// Angular momentum
  unsigned int _angular_momentum;
  /// Exponents
  real_container _exponents;
  /// Principal quantum numbers
  int_container _quantum_numbers;
  /// Contraction _orbital_coefficients
  coefficient_container _orbital_coefficients;
  /// Alpha orbital occupations
  int_container _alpha_occupations;
  /// Beta orbital occupations
  int_container _beta_occupations;

 public:
  /// Dummy constructor
  SlaterTypeAtomicShell() = default;
  /// Constructor
  SlaterTypeAtomicShell(unsigned int angular_momentum,
                        const real_container &exponents,
                        const int_container &quantum_numbers,
                        const coefficient_container &coefficients,
                        const int_container &alpha_occupations,
                        const int_container &beta_occupations)
      : _angular_momentum(angular_momentum),
        _exponents(exponents),
        _quantum_numbers(quantum_numbers),
        _orbital_coefficients(coefficients),
        _alpha_occupations(alpha_occupations),
        _beta_occupations(beta_occupations) {
    // Size checks
    assert(_exponents.size() == _quantum_numbers.size());
    assert(_exponents.size() == _orbital_coefficients.size());
    assert(_alpha_occupations.size() == _orbital_coefficients[0].size());
    assert(_beta_occupations.size() == _orbital_coefficients[0].size());
  };

  /// Evaluate number of basis functions
  size_t number_of_basis_functions() const { return _exponents.size(); }

  /// Evaluate number of orbitals
  size_t number_of_orbitals() const { return _orbital_coefficients[0].size(); }

  /// Evaluates the basis functions
  void evaluate_basis_functions(double r, double *array) {
    for(size_t ix = 0; ix < _exponents.size(); ix++) {
      double normalization =
          std::pow(2.0 * _exponents[ix], _quantum_numbers[ix] + 0.5) *
          std::tgamma(2 * _quantum_numbers[ix] + 1);
      array[ix] = normalization * std::pow(r, _quantum_numbers[ix] - 1) *
                  std::exp(-_exponents[ix] * r);
    }
  }

  /// Evaluates the basis functions; does memory allocation
  std::vector<double> evaluate_basis_functions(double r) {
    std::vector<double> bf(_exponents.size());
    evaluate_basis_functions(r, bf.data());
    return bf;
  }

  /// Evaluates the orbitals' values
  void evaluate_orbitals(const double *bf, double *orbs) {
    for(size_t iorb = 0; iorb < _alpha_occupations.size(); iorb++) {
      orbs[iorb] = 0.0;
      for(size_t ix = 0; ix < _exponents.size(); ix++)
        orbs[iorb] += bf[ix] * _orbital_coefficients[ix][iorb];
    }
  }

  /// Same but with allocations
  std::vector<double> evaluate_orbitals(double r) {
    std::vector<double> bf(evaluate_basis_functions(r));
    std::vector<double> orbs(_orbital_coefficients.size());
    evaluate_orbitals(bf.data(), orbs.data());
    return orbs;
  }

  /// Helper to evaluate electron densities
  double evaluate_density(const double *orbs, const int_container &occs) {
    double density = 0.0;
    for(size_t iorb = 0; iorb < occs.size(); iorb++) {
      double orbital_density = std::pow(orbs[iorb], 2);
      density += occs[iorb] * orbital_density;
    }
    return density;
  }

  /// Evaluates alpha electron density from computed orbitals
  double evaluate_alpha_density(const double *orbs) {
    return evaluate_density(orbs, _alpha_occupations);
  }

  /// Evaluates alpha electron density; computes orbitals
  double evaluate_alpha_density(double r) {
    std::vector<double> orbs(evaluate_orbitals(r));
    return evaluate_alpha_density(orbs.data());
  }

  /// Evaluates beta electron density from computed orbitals
  double evaluate_beta_density(const double *orbs) {
    return evaluate_density(orbs, _beta_occupations);
  }

  /// Evaluates beta electron density; computes orbitals
  double evaluate_beta_density(double r) {
    std::vector<double> orbs(evaluate_orbitals(r));
    return evaluate_beta_density(orbs.data());
  }

  /// Evaluate density gradient
  double evaluate_alpha_density_gradient(double r);
  /// Evaluate density gradient
  double evaluate_beta_density_gradient(double r);
  /// Evaluate kinetic energy density tau
  double evaluate_alpha_kinetic_energy_density(double r);
  /// Evaluate kinetic energy density tau
  double evaluate_beta_kinetic_energy_density(double r);
  /// Evaluate density Laplacian
  double evaluate_alpha_density_laplacian(double r);
  /// Evaluate density Laplacian
  double evaluate_beta_density_laplacian(double r);

  /// Return angular momentum
  auto angular_momentum() const { return _angular_momentum; }

  /// Return pointer to quantum number array
  auto *exponents_data() const { return _exponents.data(); }
  /// Fetch quantum number of i:th basis function
  auto exponent(size_t i) const { return _exponents[i]; }

  /// Return pointer to quantum number array
  auto *quantum_numbers_data() const { return _quantum_numbers.data(); }
  /// Fetch quantum number of i:th basis function
  auto quantum_number(size_t i) const { return _quantum_numbers[i]; }

  /// Return pointer to orbital coefficients for ix:th exponent
  auto *orbital_coefficients_data(int ix) const {
    return _orbital_coefficients[ix].data();
  }
  /// Fetch orbital coefficient of ix:th basis function in iorb:th orbital
  auto quantum_number(size_t ix, size_t iorb) const {
    return _orbital_coefficients[ix][iorb];
  }
};

class SlaterTypeAtom {
  /// Atomic number
  unsigned int _Z;
  /// Shells
  std::vector<SlaterTypeAtomicShell> _shells;

 public:
  /// Constructor
  SlaterTypeAtom(unsigned int Z,
                 const std::vector<SlaterTypeAtomicShell> &shells)
      : _Z(Z), _shells(shells){};
  /// Deconstructor
  ~SlaterTypeAtom(){};

  /// Evaluate density
  double evaluate_alpha_density(double r) const {
    double density = 0.0;
    for(auto shell : _shells) density += shell.evaluate_alpha_density(r);
    return density;
  }
  /// Evaluate density
  double evaluate_beta_density(double r) const {
    double density = 0.0;
    for(auto shell : _shells) density += shell.evaluate_beta_density(r);
    return density;
  }

  /// Get shells
  auto *shells() const { return _shells.data(); }
  /// Get i:th shell
  auto shell(size_t i) const { return _shells[i]; }
};
