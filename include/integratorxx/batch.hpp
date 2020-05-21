#ifndef __INCLUDED_INTEGRATORXX_BATCH_HPP__
#define __INCLUDED_INTEGRATORXX_BATCH_HPP__

#include <iterator>
#include <numeric>
#include <iomanip>

#include "quadrature.hpp"

namespace IntegratorXX {


  /**
   *  \brief Generic implementation of batching in 1D Quadratures
   *
   *  Provides an iterator which generates quadrature batches on the
   *  fly.
   */ 
  template <typename QuadratureType>
  class QuadratureBatch {

    const QuadratureType& quad; ///< Base quadrature
    const size_t          batch_sz; ///< Batch size

    /**
     *  \brief 1D Quadrature Batch Iterator
     */ 
    class iterator {

      using point_type  = typename QuadratureType::point_type;
      using weight_type = typename QuadratureType::weight_type;
      using point_container  = typename QuadratureType::point_container;
      using weight_container = typename QuadratureType::weight_container;

      const QuadratureType* quad; ///< Base quadrature

      size_t cur_indx;
      size_t batch_sz;

    public:

      using value_type  = 
        std::tuple< point_container, weight_container >;
      using pointer           = const value_type*;
      using reference         = value_type&&;
      using difference_type   = size_t;
      using iterator_catagory = std::forward_iterator_tag;
      
      iterator() = delete;

      iterator( 
        const QuadratureType &q, 
        const size_t indx,
        const size_t bsz = 1 
      ) : quad(&q), cur_indx(indx), batch_sz(bsz) { }

      iterator& operator++(){ cur_indx += std::min(batch_sz,quad->nPts()-cur_indx); return *this; }
      iterator  operator++(int){ iterator retval = *this; ++(*this); return retval; }


      iterator& operator+(int i){
        cur_indx = std::min( cur_indx + i*batch_sz, quad->nPts() );
        return *this;
      }

      bool operator==(iterator other){ return cur_indx == other.cur_indx and quad == other.quad; }
      bool operator!=(iterator other){ return not (*this == other); }

      auto range() {
      
        size_t ncpy = 
          std::min(batch_sz, quad->nPts() - cur_indx);

        auto pts_st = quad->points().begin() + cur_indx;
        auto pts_en = pts_st + ncpy;

        auto wgt_st = quad->weights().begin() + cur_indx;
        auto wgt_en = wgt_st + ncpy;

        return std::tuple( pts_st, pts_en, wgt_st, wgt_en );
                 
      }

      auto operator*(){ 

        size_t ncpy = 
          std::min(batch_sz, quad->nPts() - cur_indx);

        point_container  pts( ncpy );
        weight_container wgt( ncpy );

        std::copy( quad->points().begin() + cur_indx,
                   quad->points().begin() + cur_indx + ncpy,
                   pts.begin() );

        std::copy( quad->weights().begin() + cur_indx,
                   quad->weights().begin() + cur_indx + ncpy,
                   wgt.begin() );

        return std::tuple( std::move(pts),std::move(wgt) ); 

      }

    };

  public:

    using point_type  = typename QuadratureType::point_type;
    using weight_type = typename QuadratureType::weight_type;
    using point_container = 
      typename QuadratureType::point_container;
    using weight_container = 
      typename QuadratureType::weight_container;

    QuadratureBatch( 
      const QuadratureType &q, 
      const size_t bsz = 1  
    ): batch_sz(bsz), quad(q){ };


    iterator begin() const { return iterator(quad,0, batch_sz); }
    iterator end() const {   return iterator(quad,quad.nPts(),batch_sz); }

    inline size_t n_batches() const {
      if( batch_sz > quad.nPts() ) return 1;

      size_t n = quad.nPts() / batch_sz;
      if( quad.nPts() % batch_sz ) n++;
      return n;
    }

    inline size_t max_batch_size() const {
      return batch_sz;
    }

  };






  namespace detail {

    struct standard_combine_op {
      template <typename T, typename U, typename V>
      static T combine( const U& u, const V& v ) {
        return {u,v};
      }
    };


    struct spherical_from_radial_cart_combine_op {
      template <typename T, typename U, typename V>
      static 
        std::enable_if_t<
          std::is_same_v< cartesian_pt_t<U>, T> and
          std::is_same_v< cartesian_pt_t<U>, V>,
          T
        > combine( const U& r, const V& v ) {
        return { r*v[0], r*v[1], r*v[2] };
      }
    };

  }




  /**
   *  \brief Generic implementation of batching in 2D Quadratures
   *
   *  Provides an iterator which generates 2D quadrature batches on the
   *  fly.
   */ 
  template <typename CombinedType, typename QuadratureType1, typename QuadratureType2,
            typename CombineOp = detail::standard_combine_op >
  class QuadratureBatch2D_t {

  protected:

    using batch_iterator_1 =
      decltype( std::declval<QuadratureBatch<QuadratureType1>>().begin() );
    using batch_iterator_2 =
      decltype( std::declval<QuadratureBatch<QuadratureType2>>().begin() );




    const QuadratureBatch<QuadratureType1> q_batch_1;
    const QuadratureBatch<QuadratureType2> q_batch_2;

    const size_t batch_sz_1;
    const size_t batch_sz_2;


    class iterator {

    protected: 

      using batch_iterator_1 = 
        typename QuadratureBatch2D_t::batch_iterator_1;
      using batch_iterator_2 = 
        typename QuadratureBatch2D_t::batch_iterator_2;

      using weight_type = 
        decltype( 
          std::declval<typename QuadratureType1::weight_type>() *
          std::declval<typename QuadratureType2::weight_type>()
        );

      using weight_container = typename QuadratureType1::weight_container;
      using point_container  = std::vector<CombinedType>;

      const QuadratureBatch<QuadratureType1>& q_batch_1;
      const QuadratureBatch<QuadratureType2>& q_batch_2;

      size_t cur_indx;
      batch_iterator_1 it_1;
      batch_iterator_2 it_2;

    public:

      using value_type = 
        std::tuple< CombinedType, weight_container >;
      using pointer           = const value_type*;
      using reference         = value_type&&;
      using difference_type   = size_t;
      using iterator_catagory = std::forward_iterator_tag;

      iterator() = delete;

      iterator( 
        const decltype(q_batch_1)& q1,
        const decltype(q_batch_2)& q2,
        const size_t i1,
        const size_t i2
      ) :
        q_batch_1(q1), q_batch_2(q2),
        // Col major ordering of cartesian product
        cur_indx( i1 + i2*q1.n_batches() ),
        it_1( q1.begin() + i1 ),
        it_2( q2.begin() + i2 ){ }


      iterator& operator++(){ 
        it_1++;
        if( it_1 == q_batch_1.end() ) {
          it_1 = q_batch_1.begin();
          it_2++;
        }

        return *this;
      }
      iterator  operator++(int){ iterator retval = *this; ++(*this); return retval; }


      bool operator==(iterator other){ return it_1 == other.it_1 and it_2 == other.it_2; }
      bool operator!=(iterator other){ return not (*this == other); }

      auto operator*(){ 

        auto [pts1_st, pts1_en, wgt1_st, wgt1_en] = it_1.range();
        auto [pts2_st, pts2_en, wgt2_st, wgt2_en] = it_2.range();

        auto ncpy_1 = std::distance( pts1_st, pts1_en );
        auto ncpy_2 = std::distance( pts2_st, pts2_en );

        point_container  pts( ncpy_1 * ncpy_2 );
        weight_container wgt( ncpy_1 * ncpy_2 );


        auto pts_loop_init = std::tuple(pts2_st, pts.begin());
        for(auto [i2, indx] = pts_loop_init; i2 != pts2_en; ++i2        )
        for(auto i1 = pts1_st;               i1 != pts1_en; ++i1, ++indx)
          *indx = CombineOp::template combine<CombinedType>( *i1, *i2 );

        auto wgt_loop_init = std::tuple(wgt2_st, wgt.begin());
        for(auto [i2, indx] = wgt_loop_init; i2 != wgt2_en; ++i2        )
        for(auto i1 = wgt1_st;               i1 != wgt1_en; ++i1, ++indx)
          *indx = (*i1) * (*i2);

        return std::tuple( std::move(pts),std::move(wgt) ); 

      }
    };
    

  public:

    QuadratureBatch2D_t(
      const QuadratureType1& q1,
      const QuadratureType2& q2,
      const size_t bsz1,
      const size_t bsz2
    ) : q_batch_1(q1,bsz1), q_batch_2(q2,bsz2),
        batch_sz_1(bsz1), batch_sz_2(bsz2) { }


    iterator begin(){ return iterator(q_batch_1, q_batch_2,0,0); }
    iterator end(){ return iterator(q_batch_1, q_batch_2, 0, q_batch_2.n_batches() ); }

    inline size_t n_batches() const {
      return q_batch_1.n_batches() * q_batch_2.n_batches();
    }

    inline size_t max_batch_size() const {
      return q_batch_1.max_batch_size() *
             q_batch_2.max_batch_size();
    }

  };


  /**
   *  \brief Generic implementation of batching in 2D Quadratures
   *
   *  Provides an iterator which generates 2D quadrature batches on the
   *  fly.
   */ 
  template <typename CombinedType, typename RadialQuadrature, typename AngularQuadrature,
            typename CombineOp = detail::spherical_from_radial_cart_combine_op >
  class SphericalBatch_t : 
    public QuadratureBatch2D_t<CombinedType, RadialQuadrature, AngularQuadrature, CombineOp> {

    using batch_base = QuadratureBatch2D_t<CombinedType, RadialQuadrature, AngularQuadrature, CombineOp>;

    const cartesian_pt_t<typename RadialQuadrature::point_type> center;
    const typename RadialQuadrature::point_type                 scale;

    class iterator : public batch_base::iterator {


      const cartesian_pt_t<typename RadialQuadrature::point_type> center;
      const typename RadialQuadrature::point_type                 scale;

    public:

      using point_container = typename batch_base::iterator::point_container;
      using weight_container = typename batch_base::iterator::weight_container;

      template <typename... Args>
      iterator( const decltype(center) cen, const decltype(scale) scl, Args&&... args ) :
        center(cen), scale(scl), batch_base::iterator(std::forward<Args>(args)...){ }

      auto radius_bounds() {

        auto [pts1_st, pts1_en, wgt1_st, wgt1_en] = this->it_1.range();

        return std::tuple( *pts1_st, *pts1_en );

      }

      auto operator*(){ 

        auto [pts1_st, pts1_en, wgt1_st, wgt1_en] = this->it_1.range();
        auto [pts2_st, pts2_en, wgt2_st, wgt2_en] = this->it_2.range();

        auto ncpy_1 = std::distance( pts1_st, pts1_en );
        auto ncpy_2 = std::distance( pts2_st, pts2_en );

        point_container  pts( ncpy_1 * ncpy_2 );
        weight_container wgt( ncpy_1 * ncpy_2 );


        auto pts_loop_init = std::tuple(pts2_st, pts.begin());
        for(auto [i2, indx] = pts_loop_init; i2 != pts2_en; ++i2        )
        for(auto i1 = pts1_st;               i1 != pts1_en; ++i1, ++indx)
          *indx = CombineOp::template combine<CombinedType>( *i1, *i2 );

        auto wgt_loop_init = std::tuple(wgt2_st, wgt.begin());
        for(auto [i2, indx] = wgt_loop_init; i2 != wgt2_en; ++i2        )
        for(auto i1 = wgt1_st;               i1 != wgt1_en; ++i1, ++indx)
          *indx = (*i1) * (*i2);

        // Modify weights and points
        auto r_it = pts1_st;
        auto mod_loop_init = std::tuple( pts.begin(), wgt.begin() );
        for( auto [pi, wi] = mod_loop_init; pi != pts.end(); ++pi, ++wi ) {

          *wi *= 4*M_PI * (*r_it) * (*r_it) * scale * scale * scale; ++r_it;
          if( r_it == pts1_en ) r_it = pts1_st;

          (*pi)[0] = scale*(*pi)[0] + center[0];
          (*pi)[1] = scale*(*pi)[1] + center[1];
          (*pi)[2] = scale*(*pi)[2] + center[2];

        }

        return std::tuple( std::move(pts),std::move(wgt) ); 

      }
    };
    

  public:


    template <typename... Args>
    SphericalBatch_t( const decltype(center) cen, const decltype(scale) scl, Args&&... args ) :
      center(cen), scale(scl), batch_base(std::forward<Args>(args)...){ }

    iterator begin(){ return iterator(center,scale,this->q_batch_1, this->q_batch_2,0,0); }
    iterator end(){ return iterator(center,scale,this->q_batch_1, this->q_batch_2, 0, this->q_batch_2.n_batches() ); }


  };


  template <typename T>
  std::ostream& operator<<( std::ostream& out, const cartesian_pt_t<T>& pt ) {

    out << "{" << pt[0] << ", " << pt[1] << ", " << pt[2] << "}";

    return out;

  }

  
  template <typename CombinedType, typename RadialQuadrature, typename AngularQuadrature,
            typename CombineOp = detail::spherical_from_radial_cart_combine_op >
  class SphericalMicroBatch_t {


    using radial_point_type = typename RadialQuadrature::point_type;
    using point_type = cartesian_pt_t< radial_point_type >;

    SphericalBatch_t<CombinedType, RadialQuadrature, AngularQuadrature, CombineOp> 
      sphere_generator;


    decltype( *sphere_generator.begin() ) sphere;


    const point_type         center;
    const radial_point_type  box_dim;
    const int64_t            n_macro_subdivide;
    const int64_t            n_micro_subdivide;


    std::vector<
      std::tuple<
        point_type, // low bound
        point_type,  // up bound
        std::vector< point_type >, // points
        std::vector< double >      // weights
      >
    > micro_batches;

    void generate_micro_batches() {




      assert( (n_macro_subdivide - n_micro_subdivide) % 2 == 0 );
      const auto macro_step = box_dim / n_macro_subdivide;
      const auto micro_offset = ((n_macro_subdivide - n_micro_subdivide) / 2) * macro_step;

      const auto micro_step = micro_offset / n_micro_subdivide / 2;

      const auto& sphere_points  = std::get<0>(sphere);
      const auto& sphere_weights = std::get<1>(sphere);

      std::vector< std::tuple< point_type, double > > 
        sphere_zip(sphere_points.size());

      for(auto i = 0; i < sphere_points.size(); ++i)
        sphere_zip[i] = { sphere_points[i], sphere_weights[i] };



      const point_type box_lo_bound = {
        center[0] - 0.5*box_dim,
        center[1] - 0.5*box_dim,
        center[2] - 0.5*box_dim
      };

      const point_type box_up_bound = {
        center[0] + 0.5*box_dim,
        center[1] + 0.5*box_dim,
        center[2] + 0.5*box_dim
      };

      const point_type micro_box_lo_bound = {
        box_lo_bound[0] + micro_offset,
        box_lo_bound[1] + micro_offset,
        box_lo_bound[2] + micro_offset
      };

      const point_type micro_box_up_bound = {
        box_up_bound[0] - micro_offset,
        box_up_bound[1] - micro_offset,
        box_up_bound[2] - micro_offset
      };

/*
      std::cout << std::setprecision(5) << std::scientific;
      std::cout << "BOX DIM    = " << box_dim << std::endl;
      std::cout << "BOX BOUNDS:  " << box_lo_bound << " -> " << box_up_bound << std::endl;
      std::cout << "MACRO STEP " << macro_step << std::endl;
      std::cout << "MICRO STEP " << micro_step << std::endl;
      std::cout << "MICRO OFFSET " << micro_offset << std::endl;
      std::cout << "MICRO BOX BOUNDS:  " << micro_box_lo_bound << " -> " << micro_box_up_bound << std::endl;
*/


      constexpr auto point_in_box = 
        []( const auto& pt, const auto& box_lo, const auto& box_hi ) {

          return (pt[0] >= box_lo[0] and pt[0] < box_hi[0] ) and
                 (pt[1] >= box_lo[1] and pt[1] < box_hi[1] ) and
                 (pt[2] >= box_lo[2] and pt[2] < box_hi[2] );

        };


      auto point_in_micro_box = 
        [&]( const auto& pt ) {
          return point_in_box( pt, micro_box_lo_bound, micro_box_up_bound );
        };



      constexpr auto segregate_boxes = 
        []( const auto& lo_bnd, const auto& up_bnd, const auto nseg ) {


        std::vector< std::tuple< point_type, point_type > > boxes;
        const auto step = (up_bnd[0] - lo_bnd[0]) / nseg;
        

        for( auto i = 0; i < nseg; ++i ) {

          const auto x_lo = lo_bnd[0] + i * step;
          const auto x_hi = x_lo + step;

        for( auto j = 0; j < nseg; ++j ) {

          const auto y_lo = lo_bnd[1] + j * step;
          const auto y_hi = y_lo + step;

        for( auto k = 0; k < nseg; ++k ) {

          const auto z_lo = lo_bnd[2] + k * step;
          const auto z_hi = z_lo + step;

          const point_type lo = {x_lo, y_lo, z_lo};
          const point_type up = {x_hi, y_hi, z_hi};
        
          boxes.push_back( { lo, up } );

        }
        }
        }

        return boxes;
      };



      const size_t n_macro_dim = (n_macro_subdivide - n_micro_subdivide) / 2;
      const size_t n_micro_dim = n_micro_subdivide * 2;

      // Create big box place holders
      const auto big_boxes = segregate_boxes( box_lo_bound, box_up_bound, 3 );

      // Get the index of the center box
      const auto center_box = std::find_if( big_boxes.begin(), big_boxes.end(),
        [&]( const auto& bx ) { return point_in_box( point_type({0.,0.,0.}), std::get<0>(bx), std::get<1>(bx) ); } );

      // Loop over boxes and bisect
      for( auto bx_it = big_boxes.begin(); bx_it != big_boxes.end(); bx_it++ ) {
        const auto& bx = *bx_it;
        const size_t nseg = (bx_it == center_box) ? n_micro_dim : n_macro_dim;

        auto boxes = segregate_boxes( std::get<0>(bx), std::get<1>(bx), nseg );

        for( const auto &b : boxes )
          micro_batches.push_back(
            { 
              std::get<0>(b), std::get<1>(b),
              std::vector<point_type>(), std::vector<double>()
            }
          );

      }

      // Assign points to boxes
      for( auto iPt = 0; iPt < sphere_points.size(); ++iPt )
      for( auto& bx : micro_batches ) 
      if( point_in_box( sphere_points[iPt], std::get<0>(bx), std::get<1>(bx) ) ) {
        std::get<2>(bx).emplace_back( sphere_points[iPt] );
        std::get<3>(bx).emplace_back( sphere_weights[iPt] );
      }

      


    

/*
      // Excluded point in big box
      for( auto iPt = 0; iPt < sphere_points.size(); ++iPt )
      if( (std::find(indx_keep.begin(), indx_keep.end(), iPt) == indx_keep.end()) and
          (point_in_box( sphere_points[iPt], box_lo_bound, box_up_bound ) or
           (sphere_points[iPt][0] == box_up_bound[0]) or
           (sphere_points[iPt][1] == box_up_bound[1]) or
           (sphere_points[iPt][2] == box_up_bound[2])
          ) )
        std::cout << iPt << ", " << sphere_points[iPt] << std::endl;
*/

/*
      std::cout << "IS  = " << indx_keep.size() << std::endl;
      std::unique(indx_keep.begin(), indx_keep.end() );
      std::cout << "ISU = " << indx_keep.size() << std::endl;

      std::cout << "NB  = " << micro_batches.size() << std::endl;
      std::cout << "NP  = " <<
        std::accumulate( micro_batches.begin(), micro_batches.end(), 0,
                         [](auto x, auto y){ return x + std::get<2>(y).size(); } ) << std::endl;

*/

    }

  public:

    SphericalMicroBatch_t( 
      const RadialQuadrature&   r,
      const AngularQuadrature&  l,
      const point_type cen, 
      const radial_point_type scl, 
      const radial_point_type b_dim,
      const int64_t n_macro,
      const int64_t n_micro
    ) :
      center(cen), 
      box_dim(b_dim),
      n_macro_subdivide( n_macro ),
      n_micro_subdivide( n_micro ),
      sphere_generator(cen, scl, r,l, r.nPts(), l.nPts()){ 

      // Generate the sphere
      sphere = std::move( *sphere_generator.begin() );

      // Generate the micro batches
      generate_micro_batches();
    }



    auto begin() const { return micro_batches.begin(); }
    auto end() const { return micro_batches.end(); }

  };


#if 0
  template <typename CombinedType, typename RadialQuadrature, typename AngularQuadrature,
            typename CombineOp = detail::spherical_from_radial_cart_combine_op >
  class FastSphericalMicroBatch_t {

  public:

    using radial_point_type = typename RadialQuadrature::point_type;
    using point_type = cartesian_pt_t< radial_point_type >;

    struct spherical_microbatch {
      point_type                lo_bound;
      point_type                up_bound;
      std::vector< point_type > points;
      std::vector< double >     weights;
    };

  private:


    SphericalBatch_t<CombinedType, RadialQuadrature, AngularQuadrature, CombineOp> 
      sphere_generator;

    decltype( *sphere_generator.begin() ) sphere;

    std::vector<spherical_micro_batch> micro_batches;

    void generate_micro_batches() {

      const auto& sphere_points  = std::get<0>(sphere);
      const auto& sphere_weights = std::get<1>(sphere);
      const auto  npts           = sphere_points.size();

      

      constexpr auto point_in_box = 
        []( const auto& pt, const auto& box_lo, const auto& box_hi ) {

          return (pt[0] >= box_lo[0] and pt[0] < box_hi[0] ) and
                 (pt[1] >= box_lo[1] and pt[1] < box_hi[1] ) and
                 (pt[2] >= box_lo[2] and pt[2] < box_hi[2] );

        };

     
      constexpr auto segregate_box_27 =
        []( const auto& points, auto p_idx_begin, auto p_idx_end, const auto& box_lo, const auto& box_up ) {

        const auto x_extent = box_up[0] - box_lo[0];
        const auto y_extent = box_up[1] - box_lo[1];
        const auto z_extent = box_up[2] - box_lo[2];

        const std::vector x_part_points = {
          box_lo[0],
          box_lo[0] + (x_extent / 3.),
          box_lo[0] + (x_extent / 6.),
          box_up[0] 
        };

        const std::vector y_part_points = {
          box_lo[1],
          box_lo[1] + (y_extent / 3.),
          box_lo[1] + (y_extent / 6.),
          box_up[1] 
        };

        const std::vector z_part_points = {
          box_lo[2],
          box_lo[2] + (z_extent / 3.),
          box_lo[2] + (z_extent / 6.),
          box_up[2] 
        };

        std::vector<decltype(p_idx_begin)> partition_its = { p_idx_begin };
        partition_its.reserve(28);

        std::vector< std::pair<point_type, point_type> > box_bounds;
        box_bounds.reserve(27);

        for( auto i = 0; i < 3; ++i )
        for( auto j = 0; j < 3; ++j )
        for( auto k = 0; k < 3; ++k ) {
          point_type local_box_lo = { x_part_points[i], y_part_points[j], z_part_points[k] };
          point_type local_box_up = { x_part_points[i+1], y_part_points[j+1], z_part_points[k+1] };
          auto it = std::partition( partition_its.back(), p_idx_end,
            [&](const auto& p_idx){ return point_in_box( points[p_idx], local_box_lo, local_box_up ); }
          );
          partition_its.emplace_back(it);
          box_bounds.push_back({local_box_lo}, {local_box_up}); 
        }

        return std::tuple( partition_its, box_bounds );

      };
 
      constexpr auto segregate_box_8 =
        []( const auto& points, auto p_idx_begin, auto p_idx_end, const auto& box_lo, const auto& box_up ) {

        const auto x_extent = box_up[0] - box_lo[0];
        const auto y_extent = box_up[1] - box_lo[1];
        const auto z_extent = box_up[2] - box_lo[2];

        const std::vector x_part_points = {
          box_lo[0],
          box_lo[0] + (x_extent / 2.),
          box_up[0] 
        };

        const std::vector y_part_points = {
          box_lo[1],
          box_lo[1] + (y_extent / 2.),
          box_up[1] 
        };

        const std::vector z_part_points = {
          box_lo[2],
          box_lo[2] + (z_extent / 2.),
          box_up[2] 
        };

        std::vector<decltype(p_idx_begin)> partition_its = { p_idx_begin };
        partition_its.reserve(9);
        for( auto i = 0; i < 2; ++i )
        for( auto j = 0; j < 2; ++j )
        for( auto k = 0; k < 2; ++k ) {
          point_type local_box_lo = { x_part_points[i], y_part_points[j], z_part_points[k] };
          point_type local_box_up = { x_part_points[i+1], y_part_points[j+1], z_part_points[k+1] };
          auto it = std::partition( partition_its.back(), p_idx_end,
            [&](const auto& p_idx){ return point_in_box( points[p_idx], local_box_lo, local_box_up ); }
          );
          partition_its.emplace_back(it);
        }

        return partition_its;

      };

      auto [min_x_it, max_x_it] = 
        *std::minmax_element( sphere_points.begin(), sphere_points.end(),
        [](const auto& p1, const auto& p2 ){
          return p1[0] < p2[0];
        } );
      auto [min_y_it, max_y_it] = 
        *std::minmax_element( sphere_points.begin(), sphere_points.end(),
        [](const auto& p1, const auto& p2 ){
          return p1[1] < p2[1];
        } );
      auto [min_z_it, max_z_it] = 
        *std::minmax_element( sphere_points.begin(), sphere_points.end(),
        [](const auto& p1, const auto& p2 ){
          return p1[2] < p2[2];
        } );

      const auto max_x = *max_x_it;
      const auto max_y = *max_y_it;
      const auto max_x = *max_x_it;

      const auto min_y = *min_y_it;
      const auto min_z = *min_z_it;
      const auto min_z = *min_z_it;

      const point_type big_box_lo_bound = {
        min_x, min_y, min_z
      };
        
      const point_type big_box_up_bound = {
        max_x, max_y, max_z
      };

    
      std::vector<int64_t> point_idx( npts );
      std::iota( point_idx.begin(), point_idx.end(), 0 );

      // Get initial partition
      auto [partition_its, box_bounds] = segregate_box_27( points, point_idx.begin(), point_idx.end(), 
        big_box_lo_bound, big_box_up_bound );
      assert( partition_its.size() == 28 );
      assert( box_bounds.size() == 27 );


      std::vector<point_type> new_points( npts );
      std::vector<double>     new_weights( npts );

      for( auto i = 0; i < npts; ++i ) {
        new_points[i]  = sphere_points[point_idx[i]];
        new_weights[i] = sphere_weights[point_idx[i]];
      }
      sphere_points  = std::move( new_points );
      sphere_weights = std::move( new_weights );

      for( auto i = 0; i < box_bounds.size(); ++i ) {

        auto [box_lo, box_up] = box_bounds[i];
        auto np_local = std::distance( partition_its[i], partition_its[i+1] );
        
        spherical_microbatch

      }
/*
      const size_t max_npts = 512;
      while( std::any_of( partition_its.begin(), partition_its.end()-1,
        [&](const auto it){ return std::distance( it, it+1 ) > max_npts; }
      )) {

        const auto npart = partition_its.size();
        decltype(partition_its) new_parts;
        for(auto ipart = 0; ipart < npart-1; ++ipart ) {
          auto box_st  = partition_its[ipart];
          auto box_end = partition_its[ipart+1];
          if( std::distance( partition_its[ipart], partition_its[ipart+1] ) > 512 ) {
            new_parts = segragate_box_27( points, box_st, box_end ); 
          }
        }

      }          
*/
    }


  }; // FastSphericalMicroBatch_t
#endif








  template <typename CombinedType, typename QuadratureType1, typename QuadratureType2,
            typename CombineOp = detail::standard_combine_op >
  QuadratureBatch2D_t<CombinedType, QuadratureType1, QuadratureType2, CombineOp>
  QuadratureBatch2D(
    const QuadratureType1& q1,
    const QuadratureType2& q2,
    const size_t bsz1 = 1,
    const size_t bsz2 = 1
  ) {

    return QuadratureBatch2D_t<CombinedType, QuadratureType1, QuadratureType2, CombineOp>(q1,q2,bsz1,bsz2);

  }


  template <typename RadialQuadrature, typename point_type = typename RadialQuadrature::point_type >
  SphericalBatch_t< cartesian_pt_t<point_type>, RadialQuadrature, Lebedev<point_type> >
  SphericalBatch(
    const RadialQuadrature&     r,
    const Lebedev<point_type>&  l,
    const cartesian_pt_t<point_type> cen = {0.,0.,0.},
    const point_type scale = 1.,
    const size_t bsz1 = 1,
    const size_t bsz2 = 1
  ) {
    return
    SphericalBatch_t< cartesian_pt_t<point_type>, RadialQuadrature, Lebedev<point_type> >(cen,scale,r,l,bsz1,bsz2); 
  }

  template <typename RadialQuadrature, typename point_type = typename RadialQuadrature::point_type >
  SphericalMicroBatch_t< cartesian_pt_t<point_type>, RadialQuadrature, Lebedev<point_type> >
  SphericalMicroBatch(
    const RadialQuadrature&     r,
    const Lebedev<point_type>&  l,
    const cartesian_pt_t<point_type> cen = {0.,0.,0.},
    const point_type scale = 1.
  ) {
    return
    SphericalMicroBatch_t< cartesian_pt_t<point_type>, RadialQuadrature, Lebedev<point_type> >(
      r,l,cen,scale,66.,6,2); 
  }
};

#endif
