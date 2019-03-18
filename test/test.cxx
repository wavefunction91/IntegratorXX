#include <functional>
#include <type_traits>
#include <integratorxx/quadrature.hpp>
#include <integratorxx/batch.hpp>
#include <integratorxx/integrator.hpp>

#include <algorithm>
#include <numeric>
#include <type_traits>
#include <iostream>
#include <Eigen/Core>


template <typename T>
class ctor_print_wrapper : public T {

public:

  template <typename... Args>
  ctor_print_wrapper(Args&&... args):
    T(std::forward<Args>(args)...){

    std::cout << "DIRECT CTOR" << std::endl;

  };

  ctor_print_wrapper( const ctor_print_wrapper &other ) :
    T(other){

    std::cout << "COPY CTOR" << std::endl;
  }

  ctor_print_wrapper( ctor_print_wrapper &&other ) :
    T(std::move(other)){

    std::cout << "MOVE CTOR" << std::endl;

  }

};

template <typename T>
using vector_print = ctor_print_wrapper< std::vector<T> >;


int main() {

  using namespace IntegratorXX;

  GaussLegendre<double> quad(10,-1,1);
  GaussChebyshev1<double> quad2(10,-1,1);
  GaussChebyshev2<double> quad3(10,-1,1);
  EulerMaclaurin<double> quad4(100);
  Lebedev<double> quad5(302);



  const auto& [pts, wgt] = std::tuple( quad4.points(), quad4.weights() );


  //auto f = []( double x ){ return x*x; };
  auto f = []( double x ){ return x*x*std::exp(-x*x); };

  double res = 0.;
  for( auto i = 0; i < pts.size(); ++i )
    res += wgt[i] * f(pts[i]);

  std::cout << res-std::sqrt(M_PI)/4. << std::endl;


/*
  QuadratureBatch batch(quad,3);

  std::cout << "Raw Loop" << std::endl;
  for( auto i = 0; i < quad.points().size(); ++i )
    std::cout << quad.points()[i] << ", " << quad.weights()[i] << std::endl;

  std::cout << std::endl << "BATCHES" << std::endl;
  size_t cnt = 0;
  for( auto&& [points, weights] : batch ) {

    std::cout << "BATCH " << cnt << ": " << std::endl;
    for( auto i = 0; i < points.size(); ++i ) 
      std::cout << "  " << points[i] << ", " << weights[i] << std::endl;
    ++cnt;
  }
  std::cout << cnt << ", " << batch.n_batches() << std::endl;

  auto it = batch.begin() + 3;
  auto&& [p,w] = *it;

  std::cout << std::endl;
  if( it != batch.end() )
  for( auto i = 0; i < p.size(); ++i )
    std::cout << p[i] << ", " << w[i] << std::endl;
  else std::cout << "ENDL" << std::endl;
  
*/



  auto batch2 = QuadratureBatch2D<std::array<double,2>> ( quad, quad, 2, 3 );

  auto f2 = []( double x, double y ) {
    return x*x*y*y;
  };


  size_t cnt = 0;
  res = 0;
  for( auto&& [points, weights] : batch2 ) {

    for( auto i = 0; i < points.size(); ++i ) 
      res += weights[i] * f2( points[i][0], points[i][1] );

    cnt++;
  }

  std::cout << res - 4./9.<< std::endl;

  std::cout << cnt << ", " << batch2.n_batches() << std::endl;


  const auto &leb_w = quad5.weights();
  std::cout << std::accumulate( leb_w.begin(), leb_w.end(), 0. ) << std::endl;

  auto solid_gaussian = []( double x, double y, double z ) {
    return std::exp( -x*x - y*y - z*z );
  };

  auto solid_gaussian_shift = []( double x, double y, double z ) {
    return std::exp( -(x-4.)*(x-4.) - (y-6.)*(y-6.) - (z-1.)*(z-1.) );
  };

  res = 0;
  for(auto i = 0; i < quad5.nPts(); ++i)
  for(auto j = 0; j < quad4.nPts(); ++j) {
    auto pt = quad5.points()[i];
    for(auto &x : pt) x *= quad4.points()[j];

    auto rsq = (pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2]); 
    res += quad5.weights()[i] * quad4.weights()[j] * rsq * solid_gaussian( pt[0], pt[1], pt[2] );
  }

  res = res *4 * M_PI;

  std::cout << res - M_PI * std::sqrt(M_PI) << std::endl;
  //std::cout << res  << std::endl;

  auto sph_batch = SphericalBatch( quad4, quad5, {0.,0.,0.}, 1., 5, quad5.nPts() );
  auto sph_batch_shift = SphericalBatch( quad4, quad5, {4.,6.,1.}, 1., 5, quad5.nPts() );

  auto sph_micro_batch = SphericalMicroBatch( quad4, quad5 );


  res = 0.;
  cnt = 0;
  for( auto&& [points, weights] : sph_batch ) {

    for( auto i = 0; i < points.size(); ++i ) {
      res += weights[i] * solid_gaussian( points[i][0], points[i][1], points[i][2] );
    }
    ++cnt;
  }
  //res = res *4 * M_PI;

  std::cout << cnt << ", " << res - M_PI * std::sqrt(M_PI) << std::endl;



  res = 0;
  for( auto&& [lo, up, points, weights] : sph_micro_batch ) {

    for( auto i = 0; i < points.size(); ++i ) {
      res += weights[i] * solid_gaussian( points[i][0], points[i][1], points[i][2] );
    }
    ++cnt;
  }

  std::cout << cnt << ", " << res - M_PI * std::sqrt(M_PI) << std::endl;
 
  res = 0;
  for( auto&& [points, weights] : sph_batch_shift ) {

    for( auto i = 0; i < points.size(); ++i ) {
      res += weights[i] * solid_gaussian_shift( points[i][0], points[i][1], points[i][2] );
    }
    ++cnt;
  }

  std::cout << cnt << ", " << res - M_PI * std::sqrt(M_PI) << std::endl;



  Integrator1D int1d( quad4, 4 );
  auto batch_f = [&]( double& init, const auto& points, const auto& weights ) {
    for( auto i = 0; i < points.size(); ++i)
      init += weights[i] * f(points[i]);
  };
  auto new_result = int1d.integrate( batch_f, double(0.) ); 
  std::cout << new_result-std::sqrt(M_PI)/4. << std::endl;


  using Matrix = Eigen::MatrixXd;
  //using Matrix = ctor_print_wrapper<Eigen::MatrixXd>;

  Matrix mat_result(3,3);
  mat_result.fill(0.);

  auto batch_mat_f = [&]( Matrix& init, const auto& points,const auto& weights ) {

    
    for( auto i = 0; i < points.size(); ++i) {

      Eigen::VectorXd tmp(3);
      auto eval = f(points[i]);
      tmp(0) = eval;
      tmp(1) = eval;
      tmp(2) = eval;

      init.noalias() += weights[i] * tmp * tmp.transpose();
    }

  };


  mat_result = int1d.integrate<Matrix>( batch_mat_f, Matrix::Zero(3,3) ); 
  for( int i = 0; i < 9; ++i)
    mat_result.data()[i] -= 3.*std::sqrt(M_PI)/(8*std::pow(2.,2.5));
  std::cout << std::endl << mat_result << std::endl;

//std::function<double(double&)> test_f = [&]( double& x ){ return x;};
//std::cout << std::boolalpha << is_increment_function_sig< decltype(test_f), double, void >::value << std::endl;


  return 0;
}
