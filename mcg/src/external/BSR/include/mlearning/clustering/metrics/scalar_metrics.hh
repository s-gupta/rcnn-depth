/*
 * Clustering functions for common distance metrics on scalars.
 */
#ifndef MLEARNING__CLUSTERING__METRICS__SCALAR_METRICS_HH
#define MLEARNING__CLUSTERING__METRICS__SCALAR_METRICS_HH

#include "collections/abstract/collection.hh"
#include "collections/array_list.hh"
#include "lang/array.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"
#include "mlearning/clustering/metrics/metric.hh"

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

/*
 * Clustering functions for the L1 distance between scalars.
 */
template <typename T = double>
class scalar_L1_metric : public metric<T,T> {
public:
   /*
    * Distance computation.
    * Return the L1 distance from the item to the centroid.
    */
   T distance(const T& item, const T& centroid) const {
      return (item < centroid) ? (centroid - item) : (item - centroid);
   }

   /*
    * Centroid computation.
    * Return the weighted median.
    * Note that the cluster must not be empty.
    */
   auto_ptr<T> centroid(
      const collection<T>&      items,
      const collection<double>& weights) const
   {
      /* sort items (and put weights in corresponding order) */
      array_list<T> items_arr(items);
      array<unsigned long> idx = items_arr.sort_idx();
      array_list<double> weights_arr;
      {
         array_list<double> temp(weights);
         temp.subarray(idx, weights_arr);
      }
      /* compute cumulative weight */
      unsigned long n_items = items.size();
      array<double> cum_weight(n_items);
      double total_weight = 0;
      for (unsigned long n = 0; n < n_items; n++) {
         total_weight += weights_arr[n];
         cum_weight[n] = total_weight;
      }
      /* find bounding positions for the half weight mark */
      double half_weight = total_weight/2;
      unsigned long half_pos_right = 0;
      while ((cum_weight[half_pos_right] < half_weight) && 
             (half_pos_right < (n_items-1)))
         half_pos_right++;
      unsigned long half_pos_left = 
         (half_pos_right > 0) ? (half_pos_right - 1) : 0;
      /* check which position is preferable (bounding more weight) */
      double weight_left  = cum_weight[half_pos_left];
      double weight_right = total_weight - weight_left;
      if (weight_left > weight_right) {
         return auto_ptr<T>(new T(items_arr[half_pos_left]));
      } else if (weight_left < weight_right) {
         return auto_ptr<T>(new T(items_arr[half_pos_right]));
      } else {
         return auto_ptr<T>(
            new T((items_arr[half_pos_right] + items_arr[half_pos_left])/T(2))
         );
      }
   }
};

/*
 * Clustering functions for the L2 distance between scalars.
 */
template <typename T = double>
class scalar_L2_metric : public metric<T,T> {
public:
   /*
    * Distance computation.
    * Return the L2 distance from the item to the centroid.
    */
   T distance(const T& item, const T& centroid) const {
      return (item < centroid) ? (centroid - item) : (item - centroid);
   }

   /*
    * Centroid computation.
    * Return the weighted mean.
    * Note that the cluster must not be empty.
    */
   auto_ptr<T> centroid(
      const collection<T>&      items,
      const collection<double>& weights) const
   {
      /* create iterators over items in cluster */
      auto_ptr< iterator<T> >      i_items   = items.iter_create();
      auto_ptr< iterator<double> > i_weights = weights.iter_create();
      /* initialize total weight and centroid */
      T total_weight(i_weights->next());
      auto_ptr<T> c_ptr(new T(total_weight * (i_items->next())));
      /* compute centroid */
      T& c = *c_ptr;
      while (i_items->has_next()) {
         T curr_weight(i_weights->next());
         total_weight += curr_weight;
         c += curr_weight * (i_items->next());
      }
      c /= total_weight;
      return c_ptr;
   }
};

/*
 * Globally accessible set of useful scalar clustering metrics.
 */
template <typename T = double>
class scalar_metrics {
public:
   static const scalar_L1_metric<T>& L1_metric();
   static const scalar_L2_metric<T>& L2_metric();
};

template <typename T>
const scalar_L1_metric<T>& scalar_metrics<T>::L1_metric() {
   static const scalar_L1_metric<T>* f = new scalar_L1_metric<T>();
   return *f;
}

template <typename T>
const scalar_L2_metric<T>& scalar_metrics<T>::L2_metric() {
   static const scalar_L2_metric<T>* f = new scalar_L2_metric<T>();
   return *f;
}

} /* namespace metrics */
} /* namespace clustering */
} /* namespace mlearning */

#endif
