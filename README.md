# IntegratorXX

## Synopsis

IntegratorXX is a modern C++ library for the generation of atomic and molecular
grids for quantum chemistry calculations. It provides a uniform interface for
the generation of primitive, radial and solid angle quadratures, as well as
there combination into spherical grids.

## Design Goals 

* Provide stable, reusable, and reproducible implementations of the various atomic and
molecular grids commonly encountered in quantum chemistry calculations
* Develop a modern, modular, extensible C++ API to allow for the implementation
and validation new atomic and molecular quadrature schemes.
* Provide complete C and Python interfaces to allow reusing the implementation in projects written in these languages

## Dependencies

* CMake (3.17+)
* Modern C++ compiler (C++17 compliant)

## Major Contributors

* David Williams-Young - LBNL (dbwy at lbl dot gov)
* Susi Lehtola - University of Helsinki

## Implemented Quadratures

A note on quadrature classifications:
* Primitive quadratures are those generated on a finite bound (e.g. Gauss quadrature rules). The general software design pattern of IntegratorXX is to build up higher-order quadrature rules (e.g. radial transformation, etc) from these primitive quadratures.
* Radial quadratures are convolutions of primitive quadrature rules with a radial transformation scheme (mapping the natural domain of the primitive quadrature to positive semi-indefinite). The jacobian of the transformation *is* included while the radial component of the spherical jacobian is *not*
* Angular quadratures integrate over $S^2$ (solid angle). It is typically the case that these are manually constructed to integrate spherical harmonics up to a specific order, and are thus only compatible witch *specific* grid orders.


Below are the quadratures currently implemented in IntegratorXX. Please refer to the
source for appropriate references.

| Quadrature Name                 | Quadrature Class | Domain  | C++ Class           |
|---------------------------------|------------------|---------|---------------------|
| Gauss-Chebyshev (First Kind)    | Primitive        | (-1,1)  | `GaussChebyshev1`   |
| Gauss-Chebyshev (Second Kind)   | Primitive        | [-1,1]  | `GaussChebyshev2`   |
| Gauss-Chebyshev (Third Kind)    | Primitive        | (0,1)   | `GaussChebyshev3`   |
| Gauss-Legendre                  | Primitive        | [-1,1]  | `GaussLegendre`     |
| Gauss-Lobatto                   | Primitive        | [-1,1]  | `GaussLobatto`      |
| Trapezoid Rule                  | Primitive        | [0,1]   | `UniformTrapezoid`  |
| Becke                           | Radial           | (0,inf) | `Becke`             |
| Murray-Handy-Laming (MHL)       | Radial           | (0,inf) | `MurrayHandyLaming` |
| Mura-Knowles (MK)               | Radial           | (0,inf) | `MuraKnowles`       |
| Treutler-Ahlrichs (TA, M3 + M4) | Radial           | (0,inf) | `TreutlerAhlrichs`  |
| Ahrens-Beylkin                  | Angular          | $S^2$   | `AhrensBeylkin`     |
| Delley                          | Angular          | $S^2$   | `Delley`            |
| Lebedev-Laikov                  | Angular          | $S^2$   | `LebedevLaikov`     |
| Womersley                       | Angular          | $S^2$   | `Womersley`         |


For the generation of spherical quadratures, IntegratorXX additionally supports the following radial pruning schemes:
| Name      | Description                          | C++ Specifier             |
|-----------|--------------------------------------|---------------------------|
| Unpruned  | Do not perform pruning               | `PruningScheme::Unpruned` |
| Robust    | The Psi4 "robust" pruning scheme     | `PruningScheme::Robust`   |
| Treutler  | The Treutler-Ahlrichs pruning scheme | `PruntinScheme::Treutler` |

### A Note on Angular Quadratures

All of the currently implemented angular quadrature schemes are only compatible with *specific* grid
orders corresponding to *specific* algebraic orders of spherical harmonics they integrate exactly.
The construction of the angular grids takes the number of points as argument, and will fail if the
grid order is incompatible. As these *magic numbers* are different for each of the quadratures,
we provide a set of look-up functions which can safely produce compatible grid orders:

```
using angular_type = LebedevLaikov<double>; // FP64 LL grid, similar for other implementations
using traits = IntegratorXX::quadrature_traits<angular_type>;

auto npts  = traits::npts_by_algebraic_order(order); // Return the grid order associated with a particular algebraic order
auto order = traits::algebraic_order_by_npts(npts);  // Return the algebratic order associated with a grid order
auto next_order = traits::next_algebraic_order(order); // Return the next largest (inclusive) algebratic order compatible with `angular_type` 
```


## Example Usage

Many example usages for 1-d quadratures (i.e. primitive and radial) can be found in `test/1d_quadratures.cxx` and `test/spherical_generator.cxx`. Below is a simple invocation example for the generation of an atomic sphere via the runtime generator:
```
using namespace IntegratorXX;               // Import namespace
auto rad_spec = RadialQuad::MuraKnowles;    // MK Radial scheme
auto ang_spec = AngularQuad::AhrensBeylkin; // AH Angular scheme
size_t nrad   = 99;
size_t nang   = 372;
double rscal  = 2.0;

// Generate Grid Specification
UnprunedSphericalGridSpecification unp{
  rad_spec, nrad, rscal, ang_spec, nang
};
auto pruning_spec = create_pruned_spec(PruningScheme::Robust, unp 

// Generate Quadrature
auto sph_quad = SphericalGridFactory::generate_grid(pruning_spec);

size_t npts = sph_quad->npts();
const auto& points  = sph_quad->points();  // std::vector<std::array<double,3>>
const auto& weights = sph_quad->weights(); // std::vector<double>

```

## Contributing and Bug Reports

We welcome any and all contributions and encourage bug reports. Please use the
[Issue](https://github.com/wavefunction91/IntegratorXX/issues) and 
[Pull Request](https://github.com/wavefunction91/IntegratorXX/pulls) features as appropriate.

# License

IntegratorXX is made freely available under the terms of the 3-Clause BSD license. See
LICENSE.txt for details.
