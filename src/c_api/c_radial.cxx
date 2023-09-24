#include <integratorxx/quadratures/primitive.hpp>

#include "c_internal.h"
#include <cstddef>

template <typename QuadType, typename... Args>
auto generate_quadrature(Args&&... args) {
  using alloc_type = std::allocator<QuadType>;
  using alloc_traits = std::allocator_traits<alloc_type>;
  
  alloc_type alloc;
  auto ptr = alloc_traits::allocate(alloc, 1);
  alloc_traits::construct(alloc, ptr, std::forward<Args>(args)...);

  return ptr;
}

template <typename QuadType, typename... Args>
int intxx_generate_impl(intxx_quad_type* p, Args&&... args) {
  if(p == NULL) return INTXX_NULL_QUADPTR;
  if(p->info == NULL) return INTXX_NULL_INFOPTR;

  using quad_type = QuadType;

  int npts = p->npoints;
  if( npts <= 0 ) {
    // Return error code
  }

  if(p->_state != NULL) {
    // Check if we need to regenerate
  }

  // Allocate and construct the quadrature instance
  auto ptr = generate_quadrature<quad_type>(npts, std::forward<Args>(args)... );

  // Store the pointer
  p->_state = (void*)ptr;
    
  return INTXX_SUCCESS;

} 

template <typename QuadType>
void destroy_quadrature(void* ptr) {
  using alloc_type = std::allocator<QuadType>;
  using alloc_traits = std::allocator_traits<alloc_type>;

  alloc_type alloc;

  // Cross your fingers...
  auto quad_ptr = reinterpret_cast<QuadType*>(ptr);

  // Destroy and deallocate
  alloc_traits::destroy(alloc, quad_ptr);
  alloc_traits::deallocate(alloc, quad_ptr, 1);
}

template <typename QuadType>
int intxx_destroy_impl(intxx_quad_type* p) {
  if(p == NULL) return INTXX_NULL_QUADPTR;
  if(p->info == NULL) return INTXX_NULL_INFOPTR;
  
  using namespace IntegratorXX;
  using quad_type = QuadType;

  if(p->_state != NULL) { 
    // Destroy quadrature
    destroy_quadrature<quad_type>(p->_state);
    
    // Null out state
    p->_state = NULL;
  }
  return INTXX_SUCCESS;
}

extern "C" {

int intxx_get_uniform_info(intxx_quad_info_type* p);
int intxx_get_gaussleg_info(intxx_quad_info_type* p);
int intxx_get_gausslob_info(intxx_quad_info_type* p);
int intxx_get_gausscheb1_info(intxx_quad_info_type* p);
int intxx_get_gausscheb2_info(intxx_quad_info_type* p);
int intxx_get_gausscheb2m_info(intxx_quad_info_type* p);
int intxx_get_gausscheb3_info(intxx_quad_info_type* p);

int intxx_generate_uniform(intxx_quad_type* p);
int intxx_destroy_uniform(intxx_quad_type* p);
int intxx_generate_gaussleg(intxx_quad_type* p);
int intxx_destroy_gaussleg(intxx_quad_type* p);
int intxx_generate_gausslob(intxx_quad_type* p);
int intxx_destroy_gausslob(intxx_quad_type* p);

void intxx_default_quad_info(intxx_quad_info_type* p) {
  p->number = 0; // Invalid
  p->kind   = 0; // Invalid
  p->name   = "Default"; // Default state
  p->generate = NULL;
  p->destroy  = NULL;
}

int intxx_get_prmq_info(intxx_quad_info_type* p, int quad) {
  if(p == NULL) return INTXX_NULL_INFOPTR;

  quad = quad & INTXX_PRMQ_MASK;
  p->number = quad;
  p->kind   = INTXX_PRM_QUAD;
  switch(quad) {
    case INTXX_PRMQ_UNIFORM:
      p->name = "UNIFORM";
      return intxx_get_uniform_info(p);
    case INTXX_PRMQ_GAUSSLEG:
      p->name = "GAUSS_LEGENDRE";
      return intxx_get_gaussleg_info(p);
    case INTXX_PRMQ_GAUSSCHEB_1:
      p->name = "GAUSS_CHEBYSHEV_1";
      return intxx_get_gausscheb1_info(p);
    case INTXX_PRMQ_GAUSSCHEB_2:
      p->name = "GAUSS_CHEBYSHEV_2";
      return intxx_get_gausscheb2_info(p);
    case INTXX_PRMQ_GAUSSCHEB_2MOD:
      p->name = "GAUSS_CHEBYSHEV_2MOD";
      return intxx_get_gausscheb2m_info(p);
    case INTXX_PRMQ_GAUSSCHEB_3:
      p->name = "GAUSS_CHEBYSHEV_3";
      return intxx_get_gausscheb3_info(p);
    case INTXX_PRMQ_GAUSSLOB:
      p->name = "GAUSS_LOBATTO";
      return intxx_get_gausslob_info(p);
    default:
      return INTXX_INVALID_QUAD;
  }
}

int intxx_noparam_info(intxx_quad_info_type* p, 
  int (*g)(intxx_quad_type*),
  int (*d)(intxx_quad_type*)) {
  if(p == NULL) return INTXX_NULL_INFOPTR;

  p->ext_params.n = 0; /// No External Parameters
  p->generate = g;
  p->destroy  = d;

  return INTXX_SUCCESS;
}

int intxx_get_uniform_info(intxx_quad_info_type* p) {
  return intxx_noparam_info(p, &intxx_generate_uniform,
    &intxx_destroy_uniform);
}
int intxx_get_gausslob_info(intxx_quad_info_type* p) {
  return intxx_noparam_info(p, &intxx_generate_gausslob,
    &intxx_destroy_gausslob);
}
int intxx_get_gaussleg_info(intxx_quad_info_type* p) {
  return intxx_noparam_info(p, &intxx_generate_gaussleg,
    &intxx_destroy_gaussleg);
}
int intxx_get_gausscheb1_info(intxx_quad_info_type* p) {
  return intxx_noparam_info(p, NULL, NULL);
}
int intxx_get_gausscheb2_info(intxx_quad_info_type* p) {
  return intxx_noparam_info(p, NULL, NULL);
}
int intxx_get_gausscheb2m_info(intxx_quad_info_type* p) {
  return intxx_noparam_info(p, NULL, NULL);
}
int intxx_get_gausscheb3_info(intxx_quad_info_type* p) {
  return intxx_noparam_info(p, NULL, NULL);
}

using uniform_quad_type = IntegratorXX::UniformTrapezoid<double,double>;
int intxx_generate_uniform(intxx_quad_type* p) {
  return intxx_generate_impl<uniform_quad_type>(p, 0.0, 1.0);
}
int intxx_destroy_uniform(intxx_quad_type* p) {
  return intxx_destroy_impl<uniform_quad_type>(p);
}


#define INTXX_GD_BASIC_IMPL(cname, qname)      \
int intxx_generate_##cname(intxx_quad_type* p) { \
  return intxx_generate_impl<qname>(p);        \
}                                              \
int intxx_destroy_##cname(intxx_quad_type* p) {  \
  return intxx_destroy_impl<qname>(p);         \
}

using gaussleg_quad_type = IntegratorXX::GaussLegendre<double,double>;
using gausslob_quad_type = IntegratorXX::GaussLobatto<double,double>;
INTXX_GD_BASIC_IMPL(gaussleg, gaussleg_quad_type);
INTXX_GD_BASIC_IMPL(gausslob, gausslob_quad_type);

#undef INTXX_GD_BASIC_IMPL

}
