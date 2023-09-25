#include "c_internal.h"
#include <cstddef>

extern "C" {

void intxx_default_quad_info(intxx_quad_info_type* p) {
  p->number = 0; // Invalid
  p->kind   = 0; // Invalid
  p->name   = "Default"; // Default state
  p->generate = NULL; // Invalid
  p->destroy  = NULL; // Invalid
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

/********************************************************/
/****** Info generation for quads w rad parameters ******/
/********************************************************/

int intxx_radscal_info(intxx_quad_info_type* p, 
  int (*sn)(intxx_quad_type*,const char*,double),
  int (*gn)(intxx_quad_type*,const char*,double*),
  int (*gparm)(intxx_quad_type*),
  int (*dparm)(intxx_quad_type*),
  int (*gquad)(intxx_quad_type*),
  int (*dquad)(intxx_quad_type*)) {
  if(p == NULL) return INTXX_NULL_INFOPTR;

  static const char* names[] = {"RAD_SCAL"};
  static const char* desc[]  = {"Radial Scaling Factor"};
  p->ext_params = {
    1, names, desc, 
    sn, /* set_name */
    gn, /* get_name */
    gparm, /* generate */
    dparm  /* destroy */
  };

  p->generate = gquad;
  p->destroy  = dquad;

  return INTXX_SUCCESS;
}

}
