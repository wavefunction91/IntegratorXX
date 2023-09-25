#pragma once
#include <integratorxx/c_api.h>

#ifdef __cplusplus
extern "C" {
#endif

void intxx_default_quad_info(intxx_quad_info_type* p);
int  intxx_get_radq_info(intxx_quad_info_type* p, int quad);
int  intxx_get_angq_info(intxx_quad_info_type* p, int quad);
int  intxx_get_prmq_info(intxx_quad_info_type* p, int quad);


int intxx_noparam_info(intxx_quad_info_type* p, 
  int (*g)(intxx_quad_type*),
  int (*d)(intxx_quad_type*),
  int (*s)(intxx_quad_type*, int)); 

int intxx_radscal_info(intxx_quad_info_type* p, 
  int (*sn)(intxx_quad_type*,const char*,double),
  int (*gn)(intxx_quad_type*,const char*,double*),
  int (*gparm)(intxx_quad_type*),
  int (*dparm)(intxx_quad_type*),
  int (*gquad)(intxx_quad_type*),
  int (*dquad)(intxx_quad_type*)); 

#ifdef __cplusplus
}
#endif

