/*
 * Clustering functions for common distance metrics on matrices.
 */
#ifndef MLEARNING__CLUSTERING__METRICS__MATRIX_METRICS_HH
#define MLEARNING__CLUSTERING__METRICS__MATRIX_METRICS_HH

#include "collections/abstract/collection.hh"
#include "collections/array_list.hh"
#include "lang/array.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/matrices/functors/matrix_distance_functors.hh"
#include "math/matrices/matrix.hh"
#include "mlearning/clustering/metrics/metric.hh"
#include "mlearning/clustering/metrics/scalar_metrics.hh"

namespace mlearning {
namespace clustering {
namespace metrics {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::array_list;
using lang::array;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;
using math::matrices::functors::matrix_distance_functors;
using math::matrices::functors::matrix_L1_distance;
using math::matrices::functors::matrix_L2_distance;
using math::matrices::matrix;

/*
 * Clustering functions for the L1 distance between matrices.
 */
template <typename T = double>
class matrix_L1_metric : public metric<matrix<T>,T> {
public:
   /*
    * Distance computation.
    * Return the L1 distance from the item to the centroid.
    */
   T distance(const matrix<T>& item, const matrix<T>& centroid) const {
      static const matrix_L1_distance<T>& L1_dist =
         matrix_distance_functors<T>::L1_distance();
      return L1_dist(item, centroid);
   }

   /*
    * Centroid computation.
    * Return the matrix of weighted medians along each position. 
    * Note that the cluster must not be empty.
    */
   auto_ptr< matrix<T> > centroid(
      const collection< matrix<T> >& items,
      const collection<double>&      weights) const
   {
      /* initialize centroid */
      auto_ptr< matrix<T> > m_ptr = matrix_L1_metric<T>::initialize_centroid(
         items
      );
      /* compute each element of centroid matrix */
      matrix<T>& m = *m_ptr;
      unsigned long n_elements = m.size();
      for (unsigned long n = 0; n < n_elements; n++)
         m[n] = matrix_L1_metric<T>::weighted_median(items, weights, n);
      return m_ptr;
   }

protected:
   /*
    * Initialize centroid to matrix of appropriate size.
    */
   static auto_ptr< matrix<T> > initialize_centroid(
      const collection< matrix<T> >& items)
   {
      auto_ptr< iterator< matrix<T> > > i_items = items.iter_create();
      return auto_ptr< matrix<T> >(
         new matrix<T>(i_items->next().dimensions())
      );
   }

   /*
    * Compute the weighted median of elements in the specified position.
    * Note that the collection of items must not be empty.
    */
   static T weighted_median(
      const collection< matrix<T> >& items,
      const collection<double>&      weights,
      unsigned long                  pos)
   {
      /* extract values at current position */
      array_list<T> pos_values;
      for (auto_ptr< iterator< matrix<T> > > i = items.iter_create();
           i->has_next(); )
         pos_values.add((i->next())[pos]);
      /* compute weighted median of values */
      static const scalar_L1_metric<T>& mtrc = scalar_metrics<T>::L1_metric();
      auto_ptr<T> cntrd = mtrc.centroid(pos_values, weights);
      return *cntrd;
   }
};

/*
 * Clustering functions for the L2 distance between matrices.
 */
template <typename T = double>
class matrix_L2_metric : public metric<matrix<T>,T> {
public:
   /*
    * Distance computation.
    * Return the L2 distance from the item to the centroid.
    */
   T distance(const matrix<T>& item, const matrix<T>& centroid) const {
      static const matrix_L2_distance<T>& L2_dist =
         matrix_distance_functors<T>::L2_distance();
      return L2_dist(item, centroid);
   }

   /*
    * Centroid computation.
    * Return the weighted mean of the items in the cluster.
    * Note that the cluster must not be empty.
    */
   auto_ptr< matrix<T> > centroid(
      const collection< matrix<T> >& items,
      const collection<double>&      weights) const
   {
      /* create iterators over items in cluster */
      auto_ptr< iterator< matrix<T> > > i_items = items.iter_create();
      auto_ptr< iterator<double> > i_weights = weights.iter_create();
      /* initialize total weight and centroid */
      T total_weight(i_weights->next());
      auto_ptr< matrix<T> > m_ptr(
         new matrix<T>(total_weight * (i_items->next()))
      );
      /* compute centroid */
      matrix<T>& m = *m_ptr;
      while (i_items->has_next()) {
         T curr_weight(i_weights->next());
         total_weight += curr_weight;
         m += curr_weight * (i_items->next());
      }
      m /= total_weight;
      return m_ptr;
   }
};

/*
 * Globally accessible set of useful matrix clustering metrics.
 */
template <typename T = double>
class matrix_metrics {
public:
   static const matrix_L1_metric<T>& L1_metric();
   static const matrix_L2_metric<T>& L2_metric();
};

template <typename T>
const matrix_L1_metric<T>& matrix_metrics<T>::L1_metric() {
   static const matrix_L1_metric<T>* f = new matrix_L1_metric<T>();
   return *f;
}

template <typename T>
const matrix_L2_metric<T>& matrix_metrics<T>::L2_metric() {
   static const matrix_L2_metric<T>* f = new matrix_L2_metric<T>();
   return *f;
}

} /* namespace metrics */
} /* namespace clustering */
} /* namespace mlearning */

#endif
