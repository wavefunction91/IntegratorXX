#pragma once
#include <vector>
#include <cmath>
#include <cassert>

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
  SlaterTypeAtomicShell(unsigned int angular_momentum, const real_container & exponents, const int_container & quantum_numbers, const coefficient_container & coefficients, const int_container & alpha_occupations, const int_container & beta_occupations) : _angular_momentum(angular_momentum), _exponents(exponents), _quantum_numbers(quantum_numbers), _orbital_coefficients(coefficients), _alpha_occupations(alpha_occupations), _beta_occupations(beta_occupations) {
    // Size checks
    assert(_exponents.size() == _quantum_numbers.size());
    assert(_exponents.size() == _orbital_coefficients.size());
    assert(_alpha_occupations.size() == _orbital_coefficients[0].size());
    assert(_beta_occupations.size() == _orbital_coefficients[0].size());
  };

  /// Evaluate number of basis functions
  size_t number_of_basis_functions() const {
    return _exponents.size();
  }

  /// Evaluate number of orbitals
  size_t number_of_orbitals() const {
    return _orbital_coefficients[0].size();
  }

  /// Evaluates the basis functions
  std::vector<double> evaluate_basis_functions(double r) {
    std::vector<double> bf(_exponents.size());
    evaluate_basis_functions(r, bf.data());
    return bf;
  }

  /// Same with external storage
  void evaluate_basis_functions(double r, double *array) {
    for(size_t ix=0;ix<bf.size();ix++) {
      double normalization = std::pow(2.0*_exponents[ix],_quantum_numbers[ix]+0.5) * std::tgamma(2*_quantum_numbers[ix]+1);
      array[ix] = normalization * std::pow(r, _quantum_numbers[ix]-1) * std::exp(-_exponents[ix]*r);
    }
  }

  /// Evaluates the orbitals' values
  std::vector<double> evaluate_orbitals(double r) {
    std::vector<double> bf(evaluate_basis_functions(r));
    std::vector<double> orbs(_orbital_coefficients.size());
    for(size_t iorb=0;iorb<_alpha_occupations.size();iorb++) {
      orbs[iorb] = 0.0;
      for(size_t ix=0;ix<bf.size();ix++)
        orbs[iorb] += bf[ix]*_orbital_coefficients[ix][iorb];
    }
    return orbs;
  }

  /// Evaluate density
  std::pair<double,double> evaluate_density(double r) {
    std::pair<double,double> pair;
    pair.first = pair.second = 0.0;

    // Evaluate orbitals
    std::vector<double> orbs(evaluate_orbitals(r));
    for(size_t iorb=0;iorb<orbs.size();iorb++) {
      double orbital_density = std::pow(orbs[iorb],2);
      pair.first += _alpha_occupations[iorb] * orbital_density;
      pair.second += _beta_occupations[iorb] * orbital_density;
    }

    return pair;
  }

  /// Evaluate density gradient
  std::pair<double,double> evaluate_density_gradient(double r);
  /// Evaluate kinetic energy density tau
  std::pair<double,double> evaluate_kinetic_energy_density(double r);
  /// Evaluate density Laplacian
  std::pair<double,double> evaluate_density_laplacian(double r);

  /// Return angular momentum
  auto angular_momentum() const { return _angular_momentum; }

  /// Return pointer to quantum number array
  auto* exponents_data() const { return _exponents.data(); }
  /// Fetch quantum number of i:th basis function
  auto exponent(size_t i) const { return _exponents[i]; }

  /// Return pointer to quantum number array
  auto* quantum_numbers_data() const { return _quantum_numbers.data(); }
  /// Fetch quantum number of i:th basis function
  auto quantum_number(size_t i) const { return _quantum_numbers[i]; }

  /// Return pointer to orbital coefficients for ix:th exponent
  auto* orbital_coefficients_data(int ix) const { return _orbital_coefficients[ix].data(); }
  /// Fetch orbital coefficient of ix:th basis function in iorb:th orbital
  auto quantum_number(size_t ix, size_t iorb) const { return _orbital_coefficients[ix][iorb]; }
};

class SlaterTypeAtom {
  /// Atomic number
  unsigned int _Z;
  /// Shells
  std::vector<SlaterTypeAtomicShell> _shells;
 public:
  /// Constructor
 SlaterTypeAtom(unsigned int Z, const std::vector<SlaterTypeAtomicShell> & shells) : _Z(Z), _shells(shells) {};
  /// Deconstructor
  ~SlaterTypeAtom() {};

  /// Evaluate density
  std::pair<double,double> evaluate_density(double r) {
    std::pair<double,double> p{0.0, 0.0};
    for(auto shell: _shells) {
      std::pair<double,double> sp = shell.evaluate_density(r);
      p.first += sp.first;
      p.second += sp.second;
    }
    return p;
  }

  /// Get shells
  auto * shells() const { return _shells.data(); }
  /// Get i:th shell
  auto shell(size_t i) const { return _shells[i]; }

  /// Evaluate orbitals
  std::vector<double> evaluate_orbitals(double r);
  /// Evaluate density gradient
  std::pair<double,double> evaluate_density_gradient(double r);
  /// Evaluate kinetic energy density tau
  std::pair<double,double> evaluate_kinetic_energy_density(double r);
  /// Evaluate density Laplacian
  std::pair<double,double> evaluate_density_laplacian(double r);
};
