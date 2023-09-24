#include "c_internal.h"
#include <stdlib.h>
#include <assert.h>



int intxx_quad_init(intxx_quad_type* p, int quad) {
  // Local type defs
  typedef intxx_quad_info_type info_type;

  // Sanity check
  if(p == NULL) return INTXX_NULL_QUADPTR;
  p->info    = NULL;
  p->_state  = NULL;
  p->npoints = -1;

  if(quad < 0) return INTXX_INVALID_QUAD;

  // Determine quadrature class
  int is_prmq = quad & INTXX_PRMQ_MASK;
  int is_radq = quad & INTXX_RADQ_MASK;
  int is_angq = quad & INTXX_ANGQ_MASK;
  int is_sphq = is_radq && is_angq;

  // Passed quadrature had to be something sane
  if(!is_prmq && !is_radq && !is_angq) 
    return INTXX_INVALID_QUAD;

  // Primitive quadratures cannot be mixed with angular ones
  if(is_prmq && is_radq) return INTXX_INVALID_QUAD;
  if(is_prmq && is_angq) return INTXX_INVALID_QUAD;

  // Get info
  info_type* finfo = (info_type*)malloc(sizeof(info_type));
  intxx_default_quad_info(finfo); // set (invalid) state
  p->info = finfo;

  int error;
  if(is_prmq) {
    // Get primitive quadrature info
    error = intxx_get_prmq_info(finfo, quad);
  } else if(is_sphq) {
    // TODO: Get spherical quadrature info
  } else if(is_radq) {
    // TODO: Get radial quadrature info
  } else {
    // Angular by exclusion
    // TODO: Get angular quadrature info
  }

  return error;
}

void intxx_quad_end(intxx_quad_type* p) {
  // Stateless - just bail
  if(p == NULL || p->info == NULL) return;

  // Free info
  free((void*)p->info);
  p->info = NULL;
}



int intxx_quad_set_npts(intxx_quad_type* p, int npts) {
  if(p == NULL) return INTXX_NULL_QUADPTR;
  if(p->info == NULL) return INTXX_NULL_INFOPTR;

  // NPTS must be > 0
  if(npts <= 0) return INTXX_INVALID_ARG;

  // TODO: Handle the case when NPTS is derived from params
  p->npoints = npts;
  return INTXX_SUCCESS;
}

int intxx_quad_get_npts(intxx_quad_type* p, int* npts) {
  if(p == NULL) return INTXX_NULL_QUADPTR;
  if(p->info == NULL) return INTXX_NULL_INFOPTR;

  // Return memory must be valid
  if(npts == NULL) return INTXX_INVALID_ARG;

  // TODO: Handle the case when NPTS is derived from params
  *npts = p->npoints;
  return *npts > 0 ? INTXX_SUCCESS : INTXX_INVALID_OUT;
}

int intxx_generate_quad(intxx_quad_type* p) {
  if(p == NULL) return INTXX_NULL_QUADPTR;
  if(p->info == NULL) return INTXX_NULL_INFOPTR;

  return p->info->generate(p);
}
  
int intxx_destroy_quad(intxx_quad_type* p) {
  if(p == NULL) return INTXX_NULL_QUADPTR;
  if(p->info == NULL) return INTXX_NULL_INFOPTR;

  return p->info->destroy(p);
}
