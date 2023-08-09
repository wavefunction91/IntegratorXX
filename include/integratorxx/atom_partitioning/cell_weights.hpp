#pragma once
#include <limits>
#include <cmath>

namespace IntegratorXX {

template <typename Derived>
struct HyperbolicCellFunctionWeights {
  const double* RAB = nullptr;
  size_t LDRAB = 0;
  double weight_tol = std::numeric_limits<double>::epsilon();

  template <typename T>
  static constexpr bool early_exit(T r, T dist_nearest) { 
    return Derived::early_exit(r, dist_nearest); 
  }

  template <typename T>
  void compute_unnormalized_weight_function(size_t natoms, const T* r, T* w) {
    std::fill_n(w, natoms, 1.0);
    for(size_t i = 0; i < natoms; ++i)
    for(size_t j = 0; j < i;      ++j) 
    if( w[i] > weight_tol or w[j] > weight_tol ) {
      const auto mu = (r[i] - r[j]) / RAB[j + i*LDRAB];
      Derived::apply_s_function(mu, w[i], w[j]);
    }
  }
};

template <unsigned NBecke>
struct BeckeWeightsTraits : 
  public HyperbolicCellFunctionWeights<BeckeWeightsTraits<NBecke>> {

  template <typename T>
  static constexpr bool early_exit(T r, T dist_nearest) { return false; }

  template <typename T>
  static constexpr auto h_function(T x) { 
    return 1.5 * x - 0.5 * x*x*x;
  }

  template <typename T>
  static constexpr auto g_function(T x) {
    auto r = h_function(x);
    for(int i = 1; i < NBecke; ++i) r = h_function(r);
    return r;
  }

  template <typename T>
  static constexpr auto apply_s_function(T mu, T& wi, T& wj) { 
    const auto s = 0.5 * (1.0 - g_function(mu)); 
    wi *= s;
    wj *= 1.0 - s;
  }

};

struct SSFWeightsTraits :
  public HyperbolicCellFunctionWeights<SSFWeightsTraits> {

  // Magic number beloe Eq(14) of reference
  inline static constexpr double ssf_factor = 0.64;

  // Eq (15) of reference
  template <typename T>
  static constexpr bool early_exit(T r, T dist_nearest) {
    const auto dist_cutoff = 0.5 * (1.0 - ssf_factor) * dist_nearest;
    return r < dist_cutoff;
  }

  // Eq (14) of reference 
  template <typename T>
  static constexpr auto g_function(T mu) {
    const double s_mu  = mu / ssf_factor;
    const double s_mu2 = s_mu  * s_mu;
    const double s_mu3 = s_mu  * s_mu2;
    const double s_mu5 = s_mu3 * s_mu2;
    const double s_mu7 = s_mu5 * s_mu2;

    return (35.*(s_mu - s_mu3) + 21.*s_mu5 - 5.*s_mu7) / 16.;
  }

  template <typename T>
  static constexpr auto apply_s_function(T mu, T& wi, T& wj) { 
    if(mu <= -ssf_factor)     wj = 0.;
    else if(mu >= ssf_factor) wi = 0.;
    else {
      const auto s = 0.5 * (1.0 - g_function(mu)); 
      wi *= s;
      wj *= 1.0 - s;
    }
  }

};

template <typename WeightTraits, typename RealType>
void cell_weights(size_t npts, size_t natoms, size_t iParent, 
  const RealType* coords, size_t LDC, const RealType* points, size_t LDP, 
  RealType* weights, WeightTraits traits, RealType* WORK = nullptr, 
  size_t LDW = 0) {

  // Handle work space
  std::vector<RealType> LOCAL_WORK;
  RealType* weightFunction;
  RealType* atomDist;
  if(WORK == nullptr or LDW < 2*natoms) {
    LOCAL_WORK.resize(2*natoms);
    WORK = LOCAL_WORK.data();
  }
  weightFunction = WORK;
  atomDist = weightFunction + natoms;

  // Loop over grid points
  for(size_t ipt = 0; ipt < npts; ++ipt) {

    const auto* point = points + ipt*LDP;

    // Compute distances of each center to each point 
    for(size_t iA = 0; iA < natoms; ++iA) {
      const auto* coord_A = coords + iA * LDC;
      const auto da_x = point[0] - coord_A[0]; 
      const auto da_y = point[1] - coord_A[1]; 
      const auto da_z = point[2] - coord_A[2]; 
      atomDist[iA] = std::hypot(da_x, da_y, da_z);
    }

    auto dist_nearest = *std::min_element(atomDist, atomDist + natoms);
    if(traits.early_exit(atomDist[iParent], dist_nearest)) continue;

    // Evaluate unnormalized partition functions
    traits.compute_unnormalized_weight_function(natoms, atomDist, weightFunction);

    // Normalization
    RealType sum = 0.0;
    for(size_t iA = 0; iA < natoms; ++iA) sum += weightFunction[iA];

    // Update Weight
    weights[ipt] *= weightFunction[iParent] / sum;
    
  } // ipts loop

}

}
