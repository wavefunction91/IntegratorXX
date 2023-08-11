#pragma once
#include <vector>
#include <cmath>

class SlaterTypeAtomicShell {
  /// Angular momentum
  int angular_momentum;
  /// Exponents
  std::vector<double> exponents;
  /// Principal quantum numbers
  std::vector<int> quantum_numbers;
  /// Contraction orbital_coefficients
  std::vector<std::vector<double>> orbital_coefficients;
  /// Alpha orbital occupations
  std::vector<int> alpha_occupations;
  /// Beta orbital occupations
  std::vector<int> beta_occupations;
 public:
  /// Dummy constructor
  SlaterTypeAtomicShell() {};
  /// Constructor
  SlaterTypeAtomicShell(int angular_momentum_, const std::vector<double> & exponents_, const std::vector<int> & quantum_numbers_, const std::vector<std::vector<double>> & coefficients_, const std::vector<int> & alpha_occupations_, const std::vector<int> & beta_occupations_) : angular_momentum(angular_momentum_), exponents(exponents_), quantum_numbers(quantum_numbers_), orbital_coefficients(coefficients_), alpha_occupations(alpha_occupations_), beta_occupations(beta_occupations_) {};
  /// Deconstructor
  ~SlaterTypeAtomicShell() {};

  /// Evaluates the basis functions
  std::vector<double> evaluate_basis_functions(double r) {
    std::vector<double> bf(exponents.size());
    for(size_t ix=0;ix<bf.size();ix++) {
      double normalization = std::pow(2.0*exponents[ix],quantum_numbers[ix]+0.5) * std::tgamma(2*quantum_numbers[ix]+1);
      bf[ix] = normalization * std::pow(r, quantum_numbers[ix]-1) * std::exp(-exponents[ix]*r);
    }
    return bf;
  }

  /// Evaluates the orbitals' values
  std::vector<double> evaluate_orbitals(double r) {
    std::vector<double> bf(evaluate_basis_functions(r));
    std::vector<double> orbs(orbital_coefficients.size());
    for(size_t iorb=0;iorb<alpha_occupations.size();iorb++) {
      orbs[iorb] = 0.0;
      for(size_t ix=0;ix<bf.size();ix++)
        orbs[iorb] += bf[ix]*orbital_coefficients[ix][iorb];
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
      pair.first += alpha_occupations[iorb] * orbital_density;
      pair.second += beta_occupations[iorb] * orbital_density;
    }
    
    return pair;
  }
    
  /// Evaluate density gradient
  std::pair<double,double> evaluate_density_gradient(double r);
  /// Evaluate kinetic energy density tau
  std::pair<double,double> evaluate_kinetic_energy_density(double r);
  /// Evaluate density Laplacian
  std::pair<double,double> evaluate_density_laplacian(double r);
};

class SlaterTypeAtom {
  /// Atomic number
  int Z;
  /// Shells
  std::vector<SlaterTypeAtomicShell> shells;
 public:
  /// Constructor
 SlaterTypeAtom(int Z_, const std::vector<SlaterTypeAtomicShell> & shells_) : Z(Z_), shells(shells_) {};
  /// Deconstructor
  ~SlaterTypeAtom() {};

  /// Evaluate density
  std::pair<double,double> evaluate_density(double r) {
    std::pair<double,double> p{0.0, 0.0};
    for(auto shell: shells) {
      std::pair<double,double> sp = shell.evaluate_density(r);
      p.first += sp.first;
      p.second += sp.second;
    }
    return p;
  }
  /// Evaluate orbitals
  std::vector<double> evaluate_orbitals(double r);
  /// Evaluate density gradient
  std::pair<double,double> evaluate_density_gradient(double r);
  /// Evaluate kinetic energy density tau
  std::pair<double,double> evaluate_kinetic_energy_density(double r);
  /// Evaluate density Laplacian
  std::pair<double,double> evaluate_density_laplacian(double r);
};

