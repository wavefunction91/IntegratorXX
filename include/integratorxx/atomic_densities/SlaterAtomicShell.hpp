#pragma once
#include <cassert>
#include <cmath>
#include <integratorxx/util/factorial.hpp>
#include <vector>

class SlaterTypeAtomicShell {
  using int_container = std::vector<int>;
  using real_container = std::vector<double>;
  using coefficient_container = std::vector<real_container>;

  /// Angular momentum
  unsigned int angular_momentum_;
  /// Exponents
  real_container exponents_;
  /// Principal quantum numbers
  int_container quantum_numbers_;
  /// Contraction coefficients
  coefficient_container orbital_coefficients_;
  /// Alpha orbital occupations
  int_container alpha_occupations_;
  /// Beta orbital occupations
  int_container beta_occupations_;
  /// Normalization coefficients
  real_container normalization_;

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
      : angular_momentum_(angular_momentum),
        exponents_(exponents),
        quantum_numbers_(quantum_numbers),
        orbital_coefficients_(coefficients),
        alpha_occupations_(alpha_occupations),
        beta_occupations_(beta_occupations) {
    // Size checks
    assert(exponents_.size() == quantum_numbers_.size());
    assert(exponents_.size() == orbital_coefficients_.size());
    assert(alpha_occupations_.size() == orbital_coefficients_[0].size());
    assert(beta_occupations_.size() == orbital_coefficients_[0].size());
    // Basis function normalization factors
    normalization_.resize(exponents_.size());
    for(size_t ix = 0; ix < exponents_.size(); ix++) {
      normalization_[ix] =
          std::pow(2.0 * exponents_[ix], quantum_numbers_[ix] + 0.5) *
          factorial(2 * quantum_numbers_[ix]);
    }
  };

  /// Evaluate number of basis functions
  size_t number_of_basis_functions() const { return exponents_.size(); }

  /// Evaluate number of orbitals
  size_t number_of_orbitals() const { return orbital_coefficients_[0].size(); }

  /// Evaluates the basis functions
  void evaluate_basis_functions(double r, double *array) {
    for(size_t ix = 0; ix < exponents_.size(); ix++) {
      array[ix] = normalization_[ix] * std::pow(r, quantum_numbers_[ix] - 1) *
                  std::exp(-exponents_[ix] * r);
    }
  }

  /// Evaluates the basis functions; does memory allocation
  std::vector<double> evaluate_basis_functions(double r) {
    std::vector<double> bf(exponents_.size());
    evaluate_basis_functions(r, bf.data());
    return bf;
  }

  /// Evaluates the orbitals' values
  void evaluate_orbitals(const double *bf, double *orbs) {
    for(size_t iorb = 0; iorb < alpha_occupations_.size(); iorb++) {
      orbs[iorb] = 0.0;
      for(size_t ix = 0; ix < exponents_.size(); ix++)
        orbs[iorb] += bf[ix] * orbital_coefficients_[ix][iorb];
    }
  }

  /// Same but with allocations
  std::vector<double> evaluate_orbitals(double r) {
    std::vector<double> bf(evaluate_basis_functions(r));
    std::vector<double> orbs(orbital_coefficients_.size());
    evaluate_orbitals(bf.data(), orbs.data());
    return orbs;
  }

  /// Helper to evaluate electron densities
  double evaluate_density(const double *orbs, const int_container &occs) {
    double density = 0.0;
    for(size_t iorb = 0; iorb < occs.size(); iorb++) {
      double orbital_density = std::abs(orbs[iorb]);
      orbital_density *= orbital_density;
      density += occs[iorb] * orbital_density;
    }
    return density;
  }

  /// Evaluates alpha electron density from computed orbitals
  double evaluate_alpha_density(const double *orbs) {
    return evaluate_density(orbs, alpha_occupations_);
  }

  /// Evaluates alpha electron density; computes orbitals
  double evaluate_alpha_density(double r) {
    std::vector<double> orbs(evaluate_orbitals(r));
    return evaluate_alpha_density(orbs.data());
  }

  /// Evaluates beta electron density from computed orbitals
  double evaluate_beta_density(const double *orbs) {
    return evaluate_density(orbs, beta_occupations_);
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
  auto angular_momentum() const { return angular_momentum_; }

  /// Return pointer to quantum number array
  auto *exponents_data() const { return exponents_.data(); }
  /// Fetch quantum number of i:th basis function
  auto exponent(size_t i) const { return exponents_[i]; }

  /// Return pointer to quantum number array
  auto *quantum_numbers_data() const { return quantum_numbers_.data(); }
  /// Fetch quantum number of i:th basis function
  auto quantum_number(size_t i) const { return quantum_numbers_[i]; }

  /// Return pointer to orbital coefficients for ix:th exponent
  auto *orbital_coefficients_data(int ix) const {
    return orbital_coefficients_[ix].data();
  }
  /// Fetch orbital coefficient of ix:th basis function in iorb:th orbital
  auto quantum_number(size_t ix, size_t iorb) const {
    return orbital_coefficients_[ix][iorb];
  }
};

class SlaterTypeAtom {
  /// Atomic number
  unsigned int Z_;
  /// Shells
  std::vector<SlaterTypeAtomicShell> shells_;

 public:
  /// Constructor
  SlaterTypeAtom(unsigned int Z,
                 const std::vector<SlaterTypeAtomicShell> &shells)
      : Z_(Z), shells_(shells){};
  /// Deconstructor
  ~SlaterTypeAtom(){};

  /// Evaluate density
  double evaluate_alpha_density(double r) const {
    double density = 0.0;
    for(auto shell : shells_) density += shell.evaluate_alpha_density(r);
    return density;
  }
  /// Evaluate density
  double evaluate_beta_density(double r) const {
    double density = 0.0;
    for(auto shell : shells_) density += shell.evaluate_beta_density(r);
    return density;
  }

  /// Get shells
  auto *shells() const { return shells_.data(); }
  /// Get i:th shell
  auto shell(size_t i) const { return shells_[i]; }
};
