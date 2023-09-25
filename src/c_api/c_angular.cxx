#include <integratorxx/quadratures/s2.hpp>

#include "c_internal.h"
#include "c_api_util.hpp"


extern "C" {

using leb_quad_type    = IntegratorXX::LebedevLaikov<double>;
using delley_quad_type = IntegratorXX::Delley<double>;
using wom_quad_type    = IntegratorXX::Womersley<double>;
using ab_quad_type     = IntegratorXX::AhrensBeylkin<double>;


/*******************************************************/
/****** Forward declare C-API for ANG Quadratures ******/
/*******************************************************/

#define FWD_ANG(cname) \
int intxx_get_##cname##_info(intxx_quad_info_type* p); \
int intxx_generate_##cname(intxx_quad_type* p); \
int intxx_destroy_##cname(intxx_quad_type* p); \

FWD_ANG(leb)
FWD_ANG(delley)
FWD_ANG(wom)
FWD_ANG(ab)

#undef FWD_ANG

/*******************************************************/
/****** Runtime Generator for ANG Quadrature Info ******/
/*******************************************************/

int intxx_get_angq_info(intxx_quad_info_type* p, int quad) {
  if(p == NULL) return INTXX_NULL_INFOPTR;

  quad = quad & INTXX_ANGQ_MASK;
  p->number = quad;
  p->kind   = INTXX_ANG_QUAD;
  p->dim    = 3;
  switch(quad) {
    case INTXX_ANGQ_LEB:
      p->name = "LEBEDEV_LAIKOV";
      return intxx_get_leb_info(p);
    case INTXX_ANGQ_DEL:
      p->name = "DELLEY";
      return intxx_get_delley_info(p);
    case INTXX_ANGQ_AB:
      p->name = "AHRENS_BEYLKIN";
      return intxx_get_ab_info(p);
    case INTXX_ANGQ_WOM:
      p->name = "WOMERSLEY";
      return intxx_get_wom_info(p);
    default:
      return INTXX_INVALID_QUAD;
  }
  return INTXX_SUCCESS;
}

/**********************************************/
/****** NPTS Setters for ANG Quadratures ******/
/**********************************************/

#define INTXX_ANG_SET_NPTS(cname, nsp) \
int intxx_##cname##_set_npts(intxx_quad_type* p, int npts) {\
  using namespace IntegratorXX::detail::nsp;\
  if(algebraic_order_by_npts(npts) > 0) {\
    p->npoints = npts;\
    return INTXX_SUCCESS;\
  } else {\
    p->npoints = -1;\
    return INTXX_INVALID_ARG;\
  }\
}

INTXX_ANG_SET_NPTS(leb, lebedev)
INTXX_ANG_SET_NPTS(delley, delley)
INTXX_ANG_SET_NPTS(ab, ahrensbeylkin)
INTXX_ANG_SET_NPTS(wom, womersley)

/*************************************************/
/****** Info generation for ANG Quadratures ******/
/*************************************************/

#define INTXX_NOPARAM_GET_INFO_IMPL(cname) \
int intxx_get_##cname##_info(intxx_quad_info_type* p) { \
  return intxx_noparam_info(p, &intxx_generate_##cname, \
    &intxx_destroy_##cname, intxx_##cname##_set_npts);  \
}

INTXX_NOPARAM_GET_INFO_IMPL(leb);
INTXX_NOPARAM_GET_INFO_IMPL(delley);
INTXX_NOPARAM_GET_INFO_IMPL(ab);
INTXX_NOPARAM_GET_INFO_IMPL(wom);

#undef INTXX_NOPARAM_GET_INFO_IMPL

/*************************************************************/
/****** Allocate/Deallocate ANG Quadrature States ******/
/*************************************************************/

#define INTXX_GD_BASIC_IMPL(cname, qname, ...)   \
int intxx_generate_##cname(intxx_quad_type* p) { \
  return intxx_generate_noparam_impl<qname>(p, ##__VA_ARGS__);  \
}                                                \
int intxx_destroy_##cname(intxx_quad_type* p) {  \
  return intxx_destroy_impl<qname>(p);   \
}

INTXX_GD_BASIC_IMPL(leb, leb_quad_type);
INTXX_GD_BASIC_IMPL(delley, delley_quad_type);
INTXX_GD_BASIC_IMPL(ab, ab_quad_type);
INTXX_GD_BASIC_IMPL(wom, wom_quad_type);

#undef INTXX_GD_BASIC_IMPL
}
