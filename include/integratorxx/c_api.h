#pragma once

/*** Error codes ***/
#define INTXX_SUCCESS       0
#define INTXX_INVALID_QUAD -1
#define INTXX_NULL_QUADPTR -2
#define INTXX_NULL_INFOPTR -3
#define INTXX_INVALID_OPT  -4
#define INTXX_INVALID_ARG  -5
#define INTXX_INVALID_OUT  -5

/*** Quadrature Classes ***/
#define INTXX_PRM_QUAD 1
#define INTXX_RAD_QUAD 2
#define INTXX_ANG_QUAD 3
#define INDXX_SPH_QUAD 4

/*** Primitive Quadratures ***/
#define INTXX_PRMQ_MASK           0x0000FF
#define INTXX_PRMQ_UNIFORM        0x000001
#define INTXX_PRMQ_GAUSSLEG       0x000002
#define INTXX_PRMQ_GAUSSCHEB_1    0x000003
#define INTXX_PRMQ_GAUSSCHEB_2    0x000004
#define INTXX_PRMQ_GAUSSCHEB_2MOD 0x000005
#define INTXX_PRMQ_GAUSSCHEB_3    0x000006
#define INTXX_PRMQ_GAUSSLOB       0x000007

/*** Radial Quadratures ***/
#define INTXX_RADQ_MASK  0x00FF00
#define INTXX_RADQ_BECKE 0x000100 // Becke
#define INTXX_RADQ_MHL   0x000200 // Murray-Handy-Laming
#define INTXX_RADQ_TA    0x000300 // Treuter-Ahlrichs
#define INTXX_RADQ_MK    0x000400 // Mura-Knowles

/*** Angular (S2) Quadratures ***/
#define INTXX_ANGQ_MASK 0xFF0000
#define INTXX_ANGQ_LEB  0x010000 // Lebedev-Laikov
#define INTXX_ANGQ_DEL  0x020000 // Delley
#define INTXX_ANGQ_AB   0x030000 // Ahrens-Beylkin
#define INTXX_ANGQ_WOM  0x040000 // Womersley

#ifdef __cplusplus
extern "C" {
#endif

// Forward defs of types
struct intxx_quad_type;

typedef struct {
  int n; ///< Number of parameters
  
  const char** names;        ///< Names of the params
  const char** descriptions; ///< Long descriptions of params
  const double *values;      ///< Default values of params

  //void (*set)(struct intxx_quad_type* p, const double** v);
    ///< Set function
} intxx_quad_params_type;

typedef struct {
  int number; ///< Quadrature identifier
  int kind;   ///< Type of quadrature (PRM, RAD, ANG, SPH)

  const char* name; ///< Name of the functional, e.g. "Becke"
  // TODO References

  intxx_quad_params_type ext_params; ///< External params

  //void (*init)(struct intxx_quad_type* p);
} intxx_quad_info_type;

struct intxx_quad_type {
  int npoints;
  const intxx_quad_info_type* info;
};

typedef struct intxx_quad_type intxx_quad_type;



/**
 *  @brief Initialize a specified quadrature
 *
 *  See error code returns for how to interpret failures.
 *
 *  @param[out] p Quadrature instance of the specified type
 *  @param[in]  quad Type of the specified quadrature
 *
 *  @returns INTXX_SUCCESS:      no errors were encountered
 *           INTXX_NULL_QUADPTR: `p == NULL`
 *           INTXX_INVALID_QUAD: `quad` is invalid
 */
int  intxx_quad_init(intxx_quad_type* p, int quad);

/// Frees a quadrature instance. No throw gurantee
void intxx_quad_end (intxx_quad_type* p);

/**
 *  @brief Set the number of quadrature points
 *
 *  Sets the number of nodes for the passed quadrature 
 *  instance. Only sensible for quadratures in which
 *  the size is a free parameter (most).
 *
 *  See error code returns for how to interpret failures.
 *
 *  @param[out] p The quadrature for which to set npts.
 *  @param[in]  npts The number of points
 *
 *  @returns INTXX_SUCCESS:      no errors were encountered
 *           INTXX_NULL_QUADPTR: `p == NULL`
 *           INTXX_NULL_INFOPTR: `p->info == NULL` (uninit)
 *           INTXX_INVALID_ARG:  invalid npts (e.g. npts < 0)
 */
int intxx_quad_set_npts(intxx_quad_type* p, int npts);

/**
 *  @brief Retrieve the number of points for a quadrature
 *
 *  Retreieves the number of grid points for a specified 
 *  quadrature. In necessicary, this will be initialized 
 *  from the default parameters, otherwise it will be
 *  read from p->npoints.
 *
 *  See error code returns for how to interpret failures.
 *
 *  @param[in]   p    The quadrature for which to get npts.
 *  @param[out]  npts The number of points
 *
 *  @returns INTXX_SUCCESS:      no errors were encountered
 *           INTXX_NULL_QUADPTR: `p == NULL`
 *           INTXX_NULL_INFOPTR: `p->info == NULL` (uninit)
 *           INTXX_INVALID_ARG:  `npts == NULL` 
 *           INTXX_INVALID_OUT:  `npts` is not valid (npts < 0)
 */ 
int intxx_quad_get_npts(intxx_quad_type* p, int* npts);

#ifdef __cplusplus
}
#endif
