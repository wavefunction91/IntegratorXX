#pragma once

#include <integratorxx/quadrature.hpp>

namespace IntegratorXX {

template <typename CombineOp, typename Quad1, typename Quad2>
class ProductQuadrature { 

private:

  using traits = quadrature_traits<ProductQuadrature<CombineOp,Quad1,Quad2>>;

protected:

  std::unique_ptr<Quadrature<Quad1>> quad1_;
  std::unique_ptr<Quadrature<Quad2>> quad2_;

public:

  using point_type       = typename traits::point_type;
  using weight_type      = typename traits::weight_type;
  using point_container  = typename traits::point_container;
  using weight_container = typename traits::weight_container;

  ProductQuadrature( const Quadrature<Quad1>& q1, const Quadrature<Quad2>& q2 ) :
    quad1_( std::make_unique<Quadrature<Quad1>>(q1) ),
    quad2_( std::make_unique<Quadrature<Quad2>>(q2) ){ }
  
  ProductQuadrature( Quadrature<Quad1>&& q1, Quadrature<Quad2>&& q2 ) :
    quad1_( std::make_unique<Quadrature<Quad1>>(std::move(q1)) ),
    quad2_( std::make_unique<Quadrature<Quad2>>(std::move(q2)) ){ }


  ProductQuadrature( const ProductQuadrature& other ) :
    ProductQuadrature( *other.quad1_, *other.quad2_ ) { }
  ProductQuadrature( ProductQuadrature&& other ) noexcept :
    ProductQuadrature( std::move(*other.quad1_), std::move(*other.quad2_) ) { }

};



template <typename CombineOp, typename Quad1, typename Quad2>
struct quadrature_traits< ProductQuadrature<CombineOp,Quad1,Quad2> > {

  using quad1_traits = quadrature_traits<Quad1>;
  using quad2_traits = quadrature_traits<Quad2>;

  using quad1_point_type  = typename quad1_traits::point_type;
  using quad1_weight_type = typename quad1_traits::weight_type;
  using quad2_point_type  = typename quad2_traits::point_type;
  using quad2_weight_type = typename quad2_traits::weight_type;

  using point_type = std::invoke_result_t<
    CombineOp, const quad1_point_type&, const quad2_point_type&
  >;
  using weight_type = std::invoke_result_t<
    CombineOp, const quad1_weight_type&, const quad2_weight_type&
  >;

  using point_container  = std::vector< point_type >;
  using weight_container = std::vector< weight_type >;

  inline static std::tuple<point_container,weight_container>
    generate( const Quadrature<Quad1>& q1, const Quadrature<Quad2>& q2 ) {

    const auto npts1 = q1.npts();
    const auto npts2 = q2.npts();

    const auto npts = npts1 * npts2;
    point_container  points( npts );
    weight_container weights( npts );
    
    const auto& points1  = q1.points();
    const auto& weights1 = q1.weights();
    const auto& points2  = q2.points();
    const auto& weights2 = q2.weights();

    CombineOp op;

    for( size_t j = 0; j < npts2; ++j ) 
    for( size_t i = 0; i < npts1; ++i ) {

      const auto ij = i + j*npts1;
      points[ij]  = op(points1[i],  points2[j] );
      weights[ij] = op(weights1[i], weights2[j]);

    }

    return std::tuple( points, weights );
  }



};

}
