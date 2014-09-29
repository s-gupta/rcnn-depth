/*
 * Sample clusterer.
 *
 * This clusterer produces centroids by clustering a randomly selected subset
 * of the items.  It then uses these centroids to compute cluster assignments
 * for all of the items.
 *
 * Items are selected for inclusion in the subset to cluster with probablility
 * proportional to their weight.
 */
#ifndef MLEARNING__CLUSTERING__CLUSTERERS__GENERAL__SAMPLE_CLUSTERER_HH
#define MLEARNING__CLUSTERING__CLUSTERERS__GENERAL__SAMPLE_CLUSTERER_HH

#include "collections/abstract/collection.hh"
#include "collections/array_list.hh"
#include "collections/pointers/auto_collection.hh"
#include "lang/array.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "math/math.hh"
#include "math/random/util/randperm.hh"
#include "mlearning/clustering/clusterers/abstract/centroid_clusterer.hh"

namespace mlearning {
namespace clustering {
namespace clusterers {
namespace general {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::array_list;
using collections::pointers::auto_collection;
using lang::array;
using lang::exceptions::ex_invalid_argument;
using mlearning::clustering::clusterers::abstract::centroid_clusterer;

/*
 * Sample clusterer.
 */
template <typename T>
class sample_clusterer : public centroid_clusterer<T> {
public:
   /*
    * Constructor.
    * Create a clusterer that operates on a fixed size sample.
    */
   explicit sample_clusterer(
      const centroid_clusterer<T>&, /* underlying clusterer to employ */
      unsigned long                 /* max # items (0 for unlimited) */
   );

   /*
    * Constructor.
    * Create a clusterer that operates on a fraction of the total items.
    * However, if this fraction is less than the specified minimum sample 
    * size, then the clusterer operates on the minimum sample size.
    */
   explicit sample_clusterer(
      const centroid_clusterer<T>&, /* underlying clusterer to employ */
      double,                       /* fraction of items */
      unsigned long = 1             /* minimum sample size */
   );

   /*
    * Copy constructor.
    * Create a cluster with the same parameters.
    */
   sample_clusterer(const sample_clusterer<T>&);

   /*
    * Destructor.
    */
   virtual ~sample_clusterer();

   /*
    * Clustering.
    * Return the cluster assignments and (optionally) the cluster centroids.
    * If weights are not specified, then all items are set to have weight one.
    */
   using centroid_clusterer<T>::cluster;
   
   array<unsigned long> cluster(
      const collection<T>&,                  /* items to cluster */
      auto_collection< T, array_list<T> >&   /* returned cluster centroids */
   ) const;

   array<unsigned long> cluster(
      const collection<T>&,                  /* items to cluster */
      const collection<double>&,             /* weight of each item */
      auto_collection< T, array_list<T> >&   /* returned cluster centroids */
   ) const;
   
  /*
    * Cluster assignment.
    * Return the id of the cluster centroid closest to the item.
    * Note that the centroid array must not be empty.
    */
   unsigned long assign_cluster(
      const T&,                              /* item */
      const array_list<T>&                   /* centroids */
   ) const;
   
   /*
    * Cluster assignment.
    * For each item, return the id of the closest cluster centroid.
    * Note that the centroid array must not be empty.
    */
   array<unsigned long> assign_cluster(
      const collection<T>&,                  /* items */
      const array_list<T>&                   /* centroids */
   ) const;

protected:
   /*
    * Sample clusterer parameters.
    */
   const centroid_clusterer<T>& _clusterer; /* underlying clusterer */
   unsigned long                _max_items; /* max # items (0 for unlimited) */
   double                       _max_frac;  /* max fraction of items */
   unsigned long                _min_samp;  /* min # items for frac sample */

   /*
    * Compute size of the sample to cluster given the total number of items.
    */
   unsigned long sample_size(unsigned long) const;
};

/***************************************************************************
 * Sample clusterer implementation.
 ***************************************************************************/

/*
 * Constructor.
 * Create a clusterer that operates on a fixed size sample.
 */
template <typename T>
sample_clusterer<T>::sample_clusterer(
   const centroid_clusterer<T>& c,
   unsigned long                max_items)
 : _clusterer(c),
   _max_items(max_items),
   _max_frac(1),
   _min_samp(0)
{ }

/*
 * Constructor.
 * Create a clusterer that operates on a fraction of the total items.
 */
template <typename T>
sample_clusterer<T>::sample_clusterer(
   const centroid_clusterer<T>& c,
   double                       max_frac,
   unsigned long                min_samp)
 : _clusterer(c),
   _max_items(0),
   _max_frac(max_frac),
   _min_samp(min_samp)
{
   /* check maximum fraction argument */
   if ((_max_frac <= 0) || (_max_frac > 1)) {
      throw ex_invalid_argument(
         "fraction of items to cluster must be in (0,1]"
      );
   }
}

/*
 * Copy constructor.
 * Create a cluster with the same parameters.
 */
template <typename T>
sample_clusterer<T>::sample_clusterer(const sample_clusterer<T>& c)
 : _clusterer(c._clusterer),
   _max_items(c._max_items),
   _max_frac(c._max_frac),
   _min_samp(c._min_samp)
{ }

/*
 * Destructor.
 */
template <typename T>
sample_clusterer<T>::~sample_clusterer() {
   /* do nothing */
}

/*
 * Compute size of the sample to cluster given the total number of items.
 */
template <typename T>
unsigned long sample_clusterer<T>::sample_size(unsigned long n_items) const {
   /* initialize sample size */
   unsigned long samp_size = n_items;
   /* enforce maximum items constraint */
   if (_max_items > 0) {
      if (_max_items < samp_size) { samp_size = _max_items; }
   }
   /* enforce maximum fractional constraint */
   unsigned long n_frac = static_cast<unsigned long>(
      math::ceil(_max_frac * static_cast<double>(n_items))
   );
   if (n_frac < _min_samp) { n_frac = _min_samp; }
   if (n_frac < samp_size) { samp_size = n_frac; }
   return samp_size;
}

/*
 * Clustering.
 * Return the cluster assignments and cluster centroids.
 */
template <typename T>
array<unsigned long> sample_clusterer<T>::cluster(
   const collection<T>&                 items,
   auto_collection< T, array_list<T> >& centroids) const 
{
   /* compute number of items to cluster */
   unsigned long n_items = items.size();
   unsigned long samp_size = this->sample_size(n_items); 
   /* check if clustering the full set */
   if (samp_size == n_items) {
      return _clusterer.cluster(items, centroids);
   } else {
      /* randomize item order */
      array<unsigned long> idx_map = math::random::util::randperm(n_items);
      array_list<T> items_array;
      {
         array_list<T> temp_items(items);
         temp_items.subarray(idx_map, items_array);
      }
      /* create subset of items to cluster */
      array_list<T> items_subset;
      for (unsigned long n = 0; n < samp_size; n++)
         items_subset.add(items_array.remove_head());
      /* cluster subset */
      array<unsigned long> assign_subset = _clusterer.cluster(
         items_subset, centroids
      );
      /* compute assignments for remaining items */
      array<unsigned long> assign_array = _clusterer.assign_cluster(
         items_array, *centroids
      );
      /* combine assignment arrays */
      array<unsigned long> assign(n_items);
      for (unsigned long n = 0; n < samp_size; n++)
         assign[idx_map[n]] = assign_subset[n];
      for (unsigned long n = samp_size; n < n_items; n++)
         assign[idx_map[n]] = assign_array[n - samp_size];
      return assign;
   }
}

/*
 * Weighted clustering.
 * Return the cluster assignments and cluster centroids.
 */
template <typename T>
array<unsigned long> sample_clusterer<T>::cluster(
   const collection<T>&                 items,
   const collection<double>&            weights,
   auto_collection< T, array_list<T> >& centroids) const
{
   /* compute number of items to cluster */
   unsigned long n_items = items.size();
   unsigned long samp_size = this->sample_size(n_items); 
   /* check if clustering the full set */
   if (samp_size == n_items) {
      return _clusterer.cluster(items, weights, centroids);
   } else {
      /* randomize item order using weighted random permutation */
      array<unsigned long> idx_map = math::random::util::randperm(weights);
      array_list<T>      items_array;
      array_list<double> weights_array;
      {
         array_list<T>      temp_items(items);
         array_list<double> temp_weights(weights);
         temp_items.subarray(idx_map, items_array);
         temp_weights.subarray(idx_map, weights_array);
      }
      /* create subset of items to cluster */
      array_list<T>      items_subset;
      array_list<double> weights_subset;
      for (unsigned long n = 0; n < samp_size; n++) {
         items_subset.add(items_array.remove_head());
         weights_subset.add(weights_array.remove_head());
      }
      /* cluster subset */
      array<unsigned long> assign_subset = _clusterer.cluster(
         items_subset, weights_subset, centroids
      );
      /* compute assignments for remaining items */
      array<unsigned long> assign_array = _clusterer.assign_cluster(
         items_array, *centroids
      );
      /* combine assignment arrays */
      array<unsigned long> assign(n_items);
      for (unsigned long n = 0; n < samp_size; n++)
         assign[idx_map[n]] = assign_subset[n];
      for (unsigned long n = samp_size; n < n_items; n++)
         assign[idx_map[n]] = assign_array[n - samp_size];
      return assign;
   }
}

/*
 * Cluster assignment.
 * Return the id of the cluster centroid closest to the item.
 * Note that the centroid array must not be empty.
 */
template <typename T>
unsigned long sample_clusterer<T>::assign_cluster(
   const T&             item,
   const array_list<T>& centroids) const
{
   return _clusterer.assign_cluster(item, centroids);
}

/*
 * Cluster assignment.
 * For each item, return the id of the closest cluster centroid.
 * Note that the centroid array must not be empty.
 */
template <typename T>
array<unsigned long> sample_clusterer<T>::assign_cluster(
   const collection<T>& items,
   const array_list<T>& centroids) const
{
   return _clusterer.assign_cluster(items, centroids);
}

} /* namespace general */
} /* namespace clusterers */
} /* namespace clustering */
} /* namespace mlearning */

#endif
