#pragma once
#include <integratorxx/c_api.h>

#ifdef __cplusplus
extern "C" {
#endif

void intxx_default_quad_info(intxx_quad_info_type* p);
int  intxx_get_radq_info(intxx_quad_info_type* p, int quad);
int  intxx_get_angq_info(intxx_quad_info_type* p, int quad);
int  intxx_get_prmq_info(intxx_quad_info_type* p, int quad);

#ifdef __cplusplus
}
#endif

