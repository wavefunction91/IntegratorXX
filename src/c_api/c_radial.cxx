#include <integratorxx/quadratures/radial.hpp>

#include "c_internal.h"
#include "c_api_util.hpp"


extern "C" {


// C++ types for internal quadrature states
using becke_traits_type = IntegratorXX::BeckeRadialTraits;
using mhl_traits_type   = IntegratorXX::MurrayHandyLamingRadialTraits<2>;
using ta_traits_type    = IntegratorXX::TreutlerAhlrichsRadialTraits;
using mk_traits_type    = IntegratorXX::MuraKnowlesRadialTraits;

using becke_quad_type  = IntegratorXX::Becke<double, double>;
using mhl_quad_type    = IntegratorXX::MurrayHandyLaming<double, double>;
using ta_quad_type     = IntegratorXX::TreutlerAhlrichs<double, double>;
using mk_quad_type     = IntegratorXX::MuraKnowles<double, double>;

/**********************************************************/
/****** Forward declare C-API for Radial Quadratures ******/
/**********************************************************/

#define FWD_RAD(cname) \
int intxx_get_##cname##_quad_info(intxx_quad_info_type* p); \
int intxx_generate_##cname##_quad(intxx_quad_type* p); \
int intxx_destroy_##cname##_quad(intxx_quad_type* p); \
int intxx_generate_##cname##_quad_params(intxx_quad_type* p); \
int intxx_destroy_##cname##_quad_params(intxx_quad_type* p); \
int intxx_get_name_##cname##_quad(intxx_quad_type* p, const char* name, double *v); \
int intxx_set_name_##cname##_quad(intxx_quad_type* p, const char* name, double v);

FWD_RAD(becke)
FWD_RAD(mhl)
FWD_RAD(ta)
FWD_RAD(mk)

#undef FWD_RAD

/**********************************************************/
/****** Runtime Generator for Radial Quadrature Info ******/
/**********************************************************/

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

/****************************************************/
/****** Info generation for Radial Quadratures ******/
/****************************************************/

#define INTXX_RADSCAL_GET_INFO_IMPL(cname) \
int intxx_get_##cname##_quad_info(intxx_quad_info_type* p) { \
  return intxx_radscal_info(p, \
    intxx_set_name_##cname##_quad, \
    intxx_get_name_##cname##_quad, \
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

INTXX_RAD_TRAITS_IMPL(becke, becke_traits_type);
INTXX_RAD_TRAITS_IMPL(mhl,   mhl_traits_type  );
INTXX_RAD_TRAITS_IMPL(ta,    ta_traits_type   );
INTXX_RAD_TRAITS_IMPL(mk,    mk_traits_type   );

#undef INTXX_RAD_TRAITS_IMPL

/**********************************/
/****** Radial Quad GET_NAME ******/
/**********************************/

#define INTXX_RAD_GETSET_IMPL(cname, traits_type) \
int intxx_get_name_##cname##_quad(intxx_quad_type* p, const char* name, double *v) { \
  if(strcmp(name, "RAD_SCAL")){ return INTXX_INVALID_OPT; } \
  return intxx_get_rad_scal_impl<traits_type>(p, v);\
}\
int intxx_set_name_##cname##_quad(intxx_quad_type* p, const char* name, double v) { \
  if(strcmp(name, "RAD_SCAL")){ return INTXX_INVALID_OPT; } \
  return intxx_set_rad_scal_impl<traits_type>(p, v);\
}

INTXX_RAD_GETSET_IMPL(becke, becke_traits_type);
INTXX_RAD_GETSET_IMPL(mhl,   mhl_traits_type  );
INTXX_RAD_GETSET_IMPL(ta,    ta_traits_type   );
INTXX_RAD_GETSET_IMPL(mk,    mk_traits_type   );

#undef INTXX_RAD_GETSET_IMPL

/**********************************************************/
/****** Allocate/Deallocate Radial Quadrature States ******/
/**********************************************************/

#define INTXX_GD_RAD_IMPL(cname, qname, tname) \
int intxx_generate_##cname##_quad(intxx_quad_type* p) { \
  return intxx_generate_param_impl<qname, tname>(p); \
}\
int intxx_destroy_##cname##_quad(intxx_quad_type* p) {  \
  return intxx_destroy_impl<qname>(p);   \
}

INTXX_GD_RAD_IMPL(becke, becke_quad_type, becke_traits_type);
INTXX_GD_RAD_IMPL(mhl, mhl_quad_type, mhl_traits_type);
INTXX_GD_RAD_IMPL(ta, ta_quad_type, ta_traits_type);
INTXX_GD_RAD_IMPL(mk, mk_quad_type, mk_traits_type);

#undef INTXX_GD_RAD_IMPL

}
