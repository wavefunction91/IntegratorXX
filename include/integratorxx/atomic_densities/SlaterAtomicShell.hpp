#pragma once
#include <cassert>
#include <cmath>
#include <integratorxx/util/factorial.hpp>
#include <vector>

class SlaterTypeAtomicShell {
  using int_container = std::vector<int>;
  using real_container = std::vector<double>;

  /// Angular momentum
  unsigned int angular_momentum_;
  /// Exponents
  real_container exponents_;
  /// Principal quantum numbers
  int_container quantum_numbers_;
  /// Contraction coefficients
  real_container orbital_coefficients_;
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
                        const real_container &coefficients,
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
    assert(alpha_occupations_.size() == beta_occupations_.size());
    assert(alpha_occupations_.size() * exponents_.size() ==
           orbital_coefficients_.size());
    // Basis function normalization factors; include angular factor here for
    // simplicity
    normalization_.resize(exponents_.size());
    for(size_t ix = 0; ix < exponents_.size(); ix++) {
      normalization_[ix] =
          std::pow(2.0 * exponents_[ix], quantum_numbers_[ix] + 0.5) /
          std::sqrt(4.0 * M_PI *
                    IntegratorXX::factorial(2 * quantum_numbers_[ix]));
    }
  };

  /// Evaluate number of basis functions
  size_t number_of_basis_functions() const { return exponents_.size(); }

  /// Evaluate number of orbitals
  size_t number_of_orbitals() const { return alpha_occupations_.size(); }

  /// Evaluates the basis functions
  void evaluate_basis_functions(double r, double *array) {
    for(size_t ix = 0; ix < exponents_.size(); ix++) {
      array[ix] = normalization_[ix] * std::pow(r, quantum_numbers_[ix] - 1) *
                  std::exp(-exponents_[ix] * r);
    }
  }

  /// Evaluates the basis function gradients
  void evaluate_basis_function_gradients(double r, double *array) {
    for(size_t ix = 0; ix < exponents_.size(); ix++) {
      array[ix] = -exponents_[ix] * std::pow(r, quantum_numbers_[ix] - 1);
      if(quantum_numbers_[ix] > 0) {
        array[ix] +=
            (quantum_numbers_[ix] - 1) * std::pow(r, quantum_numbers_[ix] - 2);
      }
      array[ix] *= normalization_[ix] * std::exp(-exponents_[ix] * r);
    }
  }

  /// Evaluates the basis function second derivatives
  void evaluate_basis_function_laplacians(double r, double *array) {
    for(size_t ix = 0; ix < exponents_.size(); ix++) {
      array[ix] = exponents_[ix] * exponents_[ix] *
                  std::pow(r, quantum_numbers_[ix] - 1);
      if(quantum_numbers_[ix] > 1) {
        array[ix] -= 2.0 * exponents_[ix] * (quantum_numbers_[ix] - 1) *
                     std::pow(r, quantum_numbers_[ix] - 2);
        if(quantum_numbers_[ix] > 2) {
          array[ix] += (quantum_numbers_[ix] * quantum_numbers_[ix] -
                        3 * quantum_numbers_[ix] + 2) *
                       std::pow(r, quantum_numbers_[ix] - 3);
        }
      }
      array[ix] *= normalization_[ix] * std::exp(-exponents_[ix] * r);
    }
  }

  /// Evaluates the orbitals' values
  void evaluate_orbitals(const double *bf, double *orbs) {
    for(size_t iorb = 0; iorb < alpha_occupations_.size(); iorb++) {
      orbs[iorb] = 0.0;
      for(size_t ix = 0; ix < exponents_.size(); ix++)
        orbs[iorb] +=
            bf[ix] * orbital_coefficients_[iorb * exponents_.size() + ix];
    }
  }

  /// Helper to evaluate electron densities
  double evaluate_density(const double *orbs, const int_container &occs) {
    double density = 0.0;
    for(size_t iorb = 0; iorb < occs.size(); iorb++) {
      density += occs[iorb] * orbs[iorb] * orbs[iorb];
    }
    return density;
  }

  /// Helper to evaluate electron density gradient
  double evaluate_density_gradient(const double *orbs, const double *dorbs,
                                   const int_container &occs) {
    double gradient = 0.0;
    for(size_t iorb = 0; iorb < occs.size(); iorb++) {
      gradient += 2.0 * occs[iorb] * orbs[iorb] * dorbs[iorb];
    }
    return gradient;
  }

  /// Helper to evaluate electron density gradient
  double evaluate_tau(const double *dorbs, const int_container &occs) {
    double tau = 0.0;
    for(size_t iorb = 0; iorb < occs.size(); iorb++) {
      tau += occs[iorb] * dorbs[iorb] * dorbs[iorb];
    }
    tau *= 0.5;
    return tau;
  }

  /// Helper to evaluate electron density laplacian
  double evaluate_density_laplacian(const double *orbs, const double *dorbs,
                                    const double *lorbs,
                                    const int_container &occs) {
    double lapl = 0.0;
    for(size_t iorb = 0; iorb < occs.size(); iorb++) {
      lapl += 2.0 * occs[iorb] *
              (dorbs[iorb] * dorbs[iorb] + orbs[iorb] * lorbs[iorb]);
    }
    return lapl;
  }

  /// Evaluates alpha electron density from computed orbitals
  double evaluate_alpha_density(const double *orbs) {
    return evaluate_density(orbs, alpha_occupations_);
  }

  /// Evaluates beta electron density from computed orbitals
  double evaluate_beta_density(const double *orbs) {
    return evaluate_density(orbs, beta_occupations_);
  }

  /// Evaluates alpha electron density from computed orbitals
  double evaluate_alpha_density_gradient(const double *orbs,
                                         const double *dorbs) {
    return evaluate_density_gradient(orbs, dorbs, alpha_occupations_);
  }

  /// Evaluates beta electron density from computed orbitals
  double evaluate_beta_density_gradient(const double *orbs,
                                        const double *dorbs) {
    return evaluate_density_gradient(orbs, dorbs, beta_occupations_);
  }

  /// Evaluates alpha electron density from computed orbitals
  double evaluate_alpha_tau(const double *dorbs) {
    return evaluate_tau(dorbs, alpha_occupations_);
  }

  /// Evaluates beta electron density from computed orbitals
  double evaluate_beta_tau(const double *dorbs) {
    return evaluate_tau(dorbs, beta_occupations_);
  }

  /// Evaluates alpha electron density from computed orbitals
  double evaluate_alpha_density_laplacian(const double *orbs,
                                          const double *dorbs,
                                          const double *lorbs) {
    return evaluate_density_laplacian(orbs, dorbs, lorbs, alpha_occupations_);
  }

  /// Evaluates beta electron density from computed orbitals
  double evaluate_beta_density_laplacian(const double *orbs,
                                         const double *dorbs,
                                         const double *lorbs) {
    return evaluate_density_laplacian(orbs, dorbs, lorbs, beta_occupations_);
  }

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
  auto *orbital_coefficients_data() const {
    return orbital_coefficients_.data();
  }
  /// Fetch orbital coefficient of ix:th basis function in iorb:th orbital
  auto quantum_number(size_t ix, size_t iorb) const {
    return orbital_coefficients_[iorb * exponents_.size() + ix];
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

  /// Determine maximum number of basis functions
  auto maximum_number_of_basis_functions() {
    size_t max_bf = 0;
    for(auto shell : shells_)
      max_bf = std::max(max_bf, shell.number_of_basis_functions());
    return max_bf;
  }

  /// Determine maximum number of orbitals
  auto maximum_number_of_orbitals() {
    size_t max_orb = 0;
    for(auto shell : shells_)
      max_orb = std::max(max_orb, shell.number_of_orbitals());
    return max_orb;
  }

  /// Get shells
  auto shells() const { return shells_; }
  /// Get i:th shell
  auto shell(size_t i) const { return shells_[i]; }
};

/// Helper class for evaluating orbitals
class SlaterEvaluator {
  /// Atom
  SlaterTypeAtom atom_;
  /// Array for basis function data
  std::vector<double> bf_;
  /// Array for basis function gradient data
  std::vector<double> df_;
  /// Array for basis function laplacian data
  std::vector<double> lf_;

  /// Array for orbital data
  std::vector<double> orbs_;
  /// Array for orbital gradient data
  std::vector<double> dorbs_;
  /// Array for orbital laplacian data
  std::vector<double> lorbs_;

 public:
  SlaterEvaluator(const SlaterTypeAtom &atom) : atom_(atom) {
    bf_.resize(atom_.maximum_number_of_basis_functions());
    df_.resize(atom_.maximum_number_of_basis_functions());
    lf_.resize(atom_.maximum_number_of_basis_functions());
    orbs_.resize(atom_.maximum_number_of_orbitals());
    dorbs_.resize(atom_.maximum_number_of_orbitals());
    lorbs_.resize(atom_.maximum_number_of_orbitals());
  }

  /// Evaluate density
  double evaluate_alpha_density(double r) {
    double density = 0.0;
    for(auto shell : atom_.shells()) {
      shell.evaluate_basis_functions(r, bf_.data());
      shell.evaluate_orbitals(bf_.data(), orbs_.data());
      density += shell.evaluate_alpha_density(orbs_.data());
    }
    return density;
  }
  /// Evaluate density
  double evaluate_beta_density(double r) {
    double density = 0.0;
    for(auto shell : atom_.shells()) {
      shell.evaluate_basis_functions(r, bf_.data());
      shell.evaluate_orbitals(bf_.data(), orbs_.data());
      density += shell.evaluate_beta_density(orbs_.data());
    }
    return density;
  }
  /// Evaluate density gradient
  double evaluate_alpha_density_gradient(double r) {
    double gradient = 0.0;
    for(auto shell : atom_.shells()) {
      shell.evaluate_basis_functions(r, bf_.data());
      shell.evaluate_basis_function_gradients(r, df_.data());
      shell.evaluate_orbitals(bf_.data(), orbs_.data());
      shell.evaluate_orbitals(df_.data(), dorbs_.data());
      gradient +=
          shell.evaluate_alpha_density_gradient(orbs_.data(), dorbs_.data());
    }
    return gradient;
  }
  /// Evaluate density gradient
  double evaluate_beta_density_gradient(double r) {
    double gradient = 0.0;
    for(auto shell : atom_.shells()) {
      shell.evaluate_basis_functions(r, bf_.data());
      shell.evaluate_basis_function_gradients(r, df_.data());
      shell.evaluate_orbitals(bf_.data(), orbs_.data());
      shell.evaluate_orbitals(df_.data(), dorbs_.data());
      gradient +=
          shell.evaluate_beta_density_gradient(orbs_.data(), dorbs_.data());
    }
    return gradient;
  }
  /// Evaluate kinetic energy density
  double evaluate_alpha_tau(double r) {
    double tau = 0.0;
    for(auto shell : atom_.shells()) {
      shell.evaluate_basis_function_gradients(r, df_.data());
      shell.evaluate_orbitals(df_.data(), dorbs_.data());
      tau += shell.evaluate_alpha_tau(dorbs_.data());
    }
    return tau;
  }
  /// Evaluate kinetic energy density
  double evaluate_beta_tau(double r) {
    double tau = 0.0;
    for(auto shell : atom_.shells()) {
      shell.evaluate_basis_function_gradients(r, df_.data());
      shell.evaluate_orbitals(df_.data(), dorbs_.data());
      tau += shell.evaluate_beta_tau(dorbs_.data());
    }
    return tau;
  }
  /// Evaluate density laplacian
  double evaluate_alpha_density_laplacian(double r) {
    double lapl = 0.0;
    for(auto shell : atom_.shells()) {
      shell.evaluate_basis_functions(r, bf_.data());
      shell.evaluate_basis_function_gradients(r, df_.data());
      shell.evaluate_basis_function_laplacians(r, lf_.data());
      shell.evaluate_orbitals(bf_.data(), orbs_.data());
      shell.evaluate_orbitals(df_.data(), dorbs_.data());
      shell.evaluate_orbitals(lf_.data(), lorbs_.data());
      lapl += shell.evaluate_alpha_density_laplacian(
          orbs_.data(), dorbs_.data(), lorbs_.data());
    }
    return lapl;
  }
  /// Evaluate density laplacian
  double evaluate_beta_density_laplacian(double r) {
    double lapl = 0.0;
    for(auto shell : atom_.shells()) {
      shell.evaluate_basis_functions(r, bf_.data());
      shell.evaluate_basis_function_gradients(r, df_.data());
      shell.evaluate_basis_function_laplacians(r, lf_.data());
      shell.evaluate_orbitals(bf_.data(), orbs_.data());
      shell.evaluate_orbitals(df_.data(), dorbs_.data());
      shell.evaluate_orbitals(lf_.data(), lorbs_.data());
      lapl += shell.evaluate_beta_density_laplacian(orbs_.data(), dorbs_.data(),
                                                    lorbs_.data());
    }
    return lapl;
  }
};
