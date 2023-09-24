#include <integratorxx/quadratures/primitive.hpp>
#include <integratorxx/quadratures/radial.hpp>

#include "c_internal.h"
#include <cstddef>

template <typename DataType, typename... Args>
auto generate_quadrature(Args&&... args) {
  using alloc_type = std::allocator<DataType>;
  using alloc_traits = std::allocator_traits<alloc_type>;
  
  alloc_type alloc;
  auto ptr = alloc_traits::allocate(alloc, 1);
  alloc_traits::construct(alloc, ptr, std::forward<Args>(args)...);

  return ptr;
}

template <typename QuadType, typename... Args>
int intxx_generate_noparam_impl(intxx_quad_type* p, Args&&... args) {
  if(p == NULL) return INTXX_NULL_QUADPTR;
  if(p->info == NULL) return INTXX_NULL_INFOPTR;

  using quad_type = QuadType;

  int npts = p->npoints;
  if( npts <= 0 ) {
    // Return error code
  }

  if(p->_state_quad != NULL) {
    // Check if we need to regenerate
  }

  // Allocate and construct the quadrature instance
  auto ptr = generate_quadrature<quad_type>(npts, std::forward<Args>(args)... );

  // Store the pointer
  p->_state_quad = (void*)ptr;
    
  return INTXX_SUCCESS;

} 


template <typename QuadType, typename ParamType, typename... Args>
int intxx_generate_param_impl(intxx_quad_type* p, Args&&... args) {
  if(p == NULL) return INTXX_NULL_QUADPTR;

  if(p->_state_parm == NULL) {
    // Return error code
  }

  // Cross your fingers...
  const auto _params = reinterpret_cast<ParamType*>(p->_state_parm);
  return intxx_generate_noparam_impl<QuadType>(p, *_params, std::forward<Args>(args)... );

} 

template <typename DataType>
void destroy_quadrature(void* ptr) {
  using alloc_type = std::allocator<DataType>;
  using alloc_traits = std::allocator_traits<alloc_type>;

  alloc_type alloc;

  // Cross your fingers...
  auto quad_ptr = reinterpret_cast<DataType*>(ptr);

  // Destroy and deallocate
  alloc_traits::destroy(alloc, quad_ptr);
  alloc_traits::deallocate(alloc, quad_ptr, 1);
}

template <typename QuadType>
int intxx_destroy_impl(intxx_quad_type* p) {
  if(p == NULL) return INTXX_NULL_QUADPTR;
  
  using quad_type = QuadType;

  if(p->_state_quad != NULL) { 
    // Destroy quadrature
    destroy_quadrature<quad_type>(p->_state_quad);
    
    // Null out state
    p->_state_quad = NULL;
  }
  return INTXX_SUCCESS;
}

template <typename TraitsType, typename... Args>
int intxx_generate_params_impl(intxx_quad_type* p, Args&&... args) {
  if(p == NULL) return INTXX_NULL_QUADPTR;
  if(p->info == NULL) return INTXX_NULL_INFOPTR;

  if(p->_state_parm != NULL) {
    return INTXX_SUCCESS;
  }

  // Allocate and construct the quadrature instance
  auto ptr = generate_quadrature<TraitsType>(std::forward<Args>(args)... );

  // Store the pointer
  p->_state_parm = (void*)ptr;
    
  return INTXX_SUCCESS;
  
}

template <typename TraitsType>
int intxx_destroy_params_impl(intxx_quad_type* p) {
  if(p == NULL) return INTXX_NULL_QUADPTR;

  if(p->_state_parm != NULL) { 
    // Destroy traits object
    destroy_quadrature<TraitsType>(p->_state_parm);
    
    // Null out state
    p->_state_parm = NULL;
  }
  return INTXX_SUCCESS;
}


extern "C" {

#define FWD_PRM(cname) \
int intxx_get_##cname##_info(intxx_quad_info_type* p); \
int intxx_generate_##cname(intxx_quad_type* p); \
int intxx_destroy_##cname(intxx_quad_type* p); 

FWD_PRM(uniform)
FWD_PRM(gaussleg)
FWD_PRM(gausslob)
FWD_PRM(gausscheb1)
FWD_PRM(gausscheb2)
FWD_PRM(gausscheb2m)
FWD_PRM(gausscheb3)

#undef FWD_PRM

#define FWD_RAD(cname) \
int intxx_get_##cname##_quad_info(intxx_quad_info_type* p); \
int intxx_generate_##cname##_quad(intxx_quad_type* p); \
int intxx_destroy_##cname##_quad(intxx_quad_type* p); \
int intxx_generate_##cname##_quad_params(intxx_quad_type* p); \
int intxx_destroy_##cname##_quad_params(intxx_quad_type* p);

FWD_RAD(becke)
FWD_RAD(mhl)
FWD_RAD(ta)
FWD_RAD(mk)

#undef FWD_RAD





void intxx_default_quad_info(intxx_quad_info_type* p) {
  p->number = 0; // Invalid
  p->kind   = 0; // Invalid
  p->name   = "Default"; // Default state
  p->generate = NULL; // Invalid
  p->destroy  = NULL; // Invalid
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


int intxx_get_radq_info(intxx_quad_info_type* p, int quad) {
  if(p == NULL) return INTXX_NULL_INFOPTR;

  quad = quad & INTXX_RADQ_MASK;
  p->number = quad;
  p->kind   = INTXX_RAD_QUAD;
  switch(quad) {
    case INTXX_RADQ_BECKE:
      p->name = "BECKE";
      return intxx_get_becke_quad_info(p);
    case INTXX_RADQ_MHL:
      p->name = "MURRAY_HANDY_LAMING";
      return intxx_get_mhl_quad_info(p);
    case INTXX_RADQ_TA:
      p->name = "TREUTLER_AHLRICHS";
      return intxx_get_ta_quad_info(p);
    case INTXX_RADQ_MK:
      p->name = "MURA_KNOWLES";
      return intxx_get_mk_quad_info(p);
    default:
      return INTXX_INVALID_QUAD;
  }
  return INTXX_SUCCESS;
}

/******************************************************/
/****** Info generation for quads w/o parameters ******/
/******************************************************/

int intxx_noparam_info(intxx_quad_info_type* p, 
  int (*g)(intxx_quad_type*),
  int (*d)(intxx_quad_type*)) {
  if(p == NULL) return INTXX_NULL_INFOPTR;

  p->ext_params.n = 0; /// No External Parameters
  p->generate = g;
  p->destroy  = d;

  return INTXX_SUCCESS;
}

#define INTXX_NOPARAM_GET_INFO_IMPL(cname) \
int intxx_get_##cname##_info(intxx_quad_info_type* p) { \
  return intxx_noparam_info(p, &intxx_generate_##cname, \
    &intxx_destroy_##cname);                            \
}

INTXX_NOPARAM_GET_INFO_IMPL(uniform);
INTXX_NOPARAM_GET_INFO_IMPL(gaussleg);
INTXX_NOPARAM_GET_INFO_IMPL(gausslob);
INTXX_NOPARAM_GET_INFO_IMPL(gausscheb1);
INTXX_NOPARAM_GET_INFO_IMPL(gausscheb2);
INTXX_NOPARAM_GET_INFO_IMPL(gausscheb2m);
INTXX_NOPARAM_GET_INFO_IMPL(gausscheb3);

#undef INTXX_NOPARAM_GET_INFO_IMPL

/********************************************************/
/****** Info generation for quads w rad parameters ******/
/********************************************************/


int intxx_radscal_info(intxx_quad_info_type* p, 
  int (*gparm)(intxx_quad_type*),
  int (*dparm)(intxx_quad_type*),
  int (*gquad)(intxx_quad_type*),
  int (*dquad)(intxx_quad_type*)) {
  if(p == NULL) return INTXX_NULL_INFOPTR;

  static const char* names[] = {"RAD_SCAL"};
  static const char* desc[]  = {"Radial Scaling Factor"};
  p->ext_params = {
    1, names, desc, 
    NULL, /* set_name */
    NULL, /* get_name */
    gparm, /* generate */
    dparm  /* destroy */
  };

  p->generate = gquad;
  p->destroy  = dquad;

  return INTXX_SUCCESS;
}

#define INTXX_RADSCAL_GET_INFO_IMPL(cname) \
int intxx_get_##cname##_quad_info(intxx_quad_info_type* p) { \
  return intxx_radscal_info(p, \
    intxx_generate_##cname##_quad_params, \
    intxx_destroy_##cname##_quad_params,  \
    intxx_generate_##cname##_quad, \
    intxx_destroy_##cname##_quad  \
    ); \
}

INTXX_RADSCAL_GET_INFO_IMPL(becke)
INTXX_RADSCAL_GET_INFO_IMPL(mhl)
INTXX_RADSCAL_GET_INFO_IMPL(ta)
INTXX_RADSCAL_GET_INFO_IMPL(mk)

#undef INTXX_RADSCAL_GET_INFO_IMPL


/***********************************************/
/****** Allocate/Deallocate Radial Traits ******/
/***********************************************/

#define INTXX_RAD_TRAITS_IMPL(cname, traits_type)              \
int intxx_generate_##cname##_quad_params(intxx_quad_type* p) { \
  return intxx_generate_params_impl<traits_type>(p);           \
}                                                              \
int intxx_destroy_##cname##_quad_params(intxx_quad_type* p) {  \
  return intxx_destroy_params_impl<traits_type>(p);            \
} 

using becke_traits_type = IntegratorXX::BeckeRadialTraits;
using mhl_traits_type   = IntegratorXX::MurrayHandyLamingRadialTraits<2>;
using ta_traits_type    = IntegratorXX::TreutlerAhlrichsRadialTraits;
using mk_traits_type    = IntegratorXX::MuraKnowlesRadialTraits;

INTXX_RAD_TRAITS_IMPL(becke, becke_traits_type);
INTXX_RAD_TRAITS_IMPL(mhl,   mhl_traits_type  );
INTXX_RAD_TRAITS_IMPL(ta,    ta_traits_type   );
INTXX_RAD_TRAITS_IMPL(mk,    mk_traits_type   );

#undef INTXX_RAD_TRAITS_IMPL




#define INTXX_GD_BASIC_IMPL(cname, qname, ...)   \
int intxx_generate_##cname(intxx_quad_type* p) { \
  return intxx_generate_noparam_impl<qname>(p, ##__VA_ARGS__);  \
}                                                \
int intxx_destroy_##cname(intxx_quad_type* p) {  \
  return intxx_destroy_impl<qname>(p);   \
}

using uniform_quad_type     = IntegratorXX::UniformTrapezoid<double,double>;
using gaussleg_quad_type    = IntegratorXX::GaussLegendre<double,double>;
using gausslob_quad_type    = IntegratorXX::GaussLobatto<double,double>;
using gausscheb1_quad_type  = IntegratorXX::GaussChebyshev1<double,double>;
using gausscheb2_quad_type  = IntegratorXX::GaussChebyshev2<double,double>;
using gausscheb3_quad_type  = IntegratorXX::GaussChebyshev3<double,double>;
using gausscheb2m_quad_type = IntegratorXX::GaussChebyshev2Modified<double,double>;

INTXX_GD_BASIC_IMPL(uniform,      uniform_quad_type, 0.0, 1.0);
INTXX_GD_BASIC_IMPL(gaussleg,      gaussleg_quad_type     );
INTXX_GD_BASIC_IMPL(gausslob,      gausslob_quad_type     );
INTXX_GD_BASIC_IMPL(gausscheb1,    gausscheb1_quad_type   );
INTXX_GD_BASIC_IMPL(gausscheb2,    gausscheb2_quad_type   );
INTXX_GD_BASIC_IMPL(gausscheb3,    gausscheb3_quad_type   );
INTXX_GD_BASIC_IMPL(gausscheb2m, gausscheb2m_quad_type);

#undef INTXX_GD_BASIC_IMPL

#define INTXX_GD_RAD_IMPL(cname, qname, tname) \
int intxx_generate_##cname##_quad(intxx_quad_type* p) { \
  return intxx_generate_param_impl<qname, tname>(p); \
}\
int intxx_destroy_##cname##_quad(intxx_quad_type* p) {  \
  return intxx_destroy_impl<qname>(p);   \
}

using becke_quad_type     = IntegratorXX::Becke<double, double>;
using mhl_quad_type     = IntegratorXX::MurrayHandyLaming<double, double>;
using ta_quad_type     = IntegratorXX::TreutlerAhlrichs<double, double>;
using mk_quad_type     = IntegratorXX::MuraKnowles<double, double>;

INTXX_GD_RAD_IMPL(becke, becke_quad_type, becke_traits_type);
INTXX_GD_RAD_IMPL(mhl, mhl_quad_type, mhl_traits_type);
INTXX_GD_RAD_IMPL(ta, ta_quad_type, ta_traits_type);
INTXX_GD_RAD_IMPL(mk, mk_quad_type, mk_traits_type);

}
