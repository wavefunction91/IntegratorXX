#pragma once
#include "c_internal.h"
#include <cstddef>
#include <cstring>

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

template <typename TraitsType>
int intxx_get_rad_scal_impl(intxx_quad_type* p, double* R) {
  if(p == NULL) return INTXX_NULL_QUADPTR;
  if(p->info == NULL) return INTXX_NULL_INFOPTR;

  if(R == NULL) {
    // Return error code
  }

  auto ext_param = p->info->ext_params;

  // Check that this quadrature has a radial scaling factor
  bool r_found = false;
  for(int i = 0; i < ext_param.n; ++i) {
    if(r_found) break;
    const auto name = ext_param.names[i];
    r_found = strcmp(name, "RAD_SCAL");
  }

  if(not r_found) {
    // Return error code
  }
  
  if(p->_state_parm == NULL) {
    // Return error code
  }

  // Read data
  *R = reinterpret_cast<TraitsType*>(p->_state_parm)->R();

  return INTXX_SUCCESS;
}

template <typename TraitsType>
int intxx_set_rad_scal_impl(intxx_quad_type* p, double R) {
  if(p == NULL) return INTXX_NULL_QUADPTR;
  if(p->info == NULL) return INTXX_NULL_INFOPTR;

  if(R <= 0.0) {
    // Return error code
  }

  auto ext_param = p->info->ext_params;

  // Check that this quadrature has a radial scaling factor
  bool r_found = false;
  for(int i = 0; i < ext_param.n; ++i) {
    if(r_found) break;
    const auto name = ext_param.names[i];
    r_found = strcmp(name, "RAD_SCAL");
  }

  if(not r_found) {
    // Return error code
  }
  
  if(p->_state_parm == NULL) {
    // Return error code
  }

  // Overwrite
  *reinterpret_cast<TraitsType*>(p->_state_parm) = TraitsType(R);

  // Regenerate if quadrature is populated
  if(p->_state_quad) {
    intxx_destroy_quad(p);
    intxx_generate_quad(p);
  }

  return INTXX_SUCCESS;
}


