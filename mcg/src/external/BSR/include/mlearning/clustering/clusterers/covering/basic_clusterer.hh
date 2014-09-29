/*
 * Basic covering clusterer.
 *
 * This template class implements a parallel clustering algorithm that works 
 * with any valid distance metric.  Clusters are closed spheres of a fixed 
 * radius centered at a sampled item.  Items are sampled for use as cluster 
 * centers with probability proportional to their weight.  Clustering is 
 * complete once all items are covered or the maximum number of samples have 
 * been drawn.
 *
 * To simply return a weighted random sample of the items, set the cluster 
 * radius to zero.
 *
 * To guarantee that all items are covered (within the radius of a cluster 
 * center), set the maximum number of samples to unlimited.
 */
#ifndef MLEARNING__CLUSTERING__CLUSTERERS__COVERING__BASIC_CLUSTERER_HH
#define MLEARNING__CLUSTERING__CLUSTERERS__COVERING__BASIC_CLUSTERER_HH

#include "collections/abstract/collection.hh"
#include "collections/array_list.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "lang/array.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/random/util/randperm.hh"
#include "mlearning/clustering/clusterers/abstract/metric_clusterer.hh"
#include "mlearning/clustering/metrics/metric.hh"


namespace mlearning {
namespace clustering {
namespace clusterers {
namespace covering {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::array_list;
using collections::list;
using collections::pointers::auto_collection;
using lang::array;
using lang::pointers::auto_ptr;
using mlearning::clustering::clusterers::abstract::metric_clusterer;
using mlearning::clustering::metrics::metric;

/*
 * Basic covering clusterer.
 * T is the type of object being clustered.
 * V is the numeric type used to represent distance.
 */
template <typename T, typename V = double>
class basic_clusterer : public metric_clusterer<T,V> {
public:
   /*
    * Constructor.
    */
   explicit basic_clusterer(
      V,                   /* cluster radius */
      unsigned long,       /* maximum # of samples (0 for unlimited) */
      const metric<T,V>&   /* metric to use for clustering */
   );

   /*
    * Copy constructor.
    * Create a clusterer with the same parameters.
    */
   basic_clusterer(const basic_clusterer<T,V>&);

   /*
    * Destructor.
    */
   virtual ~basic_clusterer();

   /*
    * Clustering.
    *
    * Return the cluster assignments and (optionally) the cluster centroids.
    *
    * The algorithm iteratively assigns items to existing clusters (processing
    * items in parallel) and then samples an item not covered for use as a new
    * cluster center.
    *
    * The set of clusters is initially empty and items are sampled with 
    * probability proportional to their weight.  Because the sampling order 
    * is randomized, different calls of this method on identical arguments may
    * produce different results.
    *
    * If weights are not specified, then all items are set to have weight one.
    *
    * Any cluster which becomes empty is dropped from the final result.
    */
   using metric_clusterer<T,V>::cluster;
   
   virtual array<unsigned long> cluster(
      const collection<T>&,                  /* items to cluster */
      auto_collection< T, array_list<T> >&   /* returned cluster centroids */
   ) const;

   virtual array<unsigned long> cluster(
      const collection<T>&,                  /* items to cluster */
      const collection<double>&,             /* weight of each item */
      auto_collection< T, array_list<T> >&   /* returned cluster centroids */
   ) const;

protected:
   /*
    * Cluster parameters.
    */
   V             _radius;        /* cluster radius */
   unsigned long _max_samples;   /* max # of samples (0 for unlimited) */

   /*
    * Weighted clustering using the specified random permutation.
    * Return the cluster assignments and cluster centroids.
    */
   array<unsigned long> cluster(
      const collection<T>&,                  /* items to cluster */
      const collection<double>&,             /* weight of each item */
      const array<unsigned long>&,           /* random permutation */
      auto_collection< T, array_list<T> >&   /* returned cluster centroids */
   ) const;
};

/***************************************************************************
 * Basic covering clusterer implementation.
 ***************************************************************************/

/*
 * Constructor.
 */
template <typename T, typename V>
basic_clusterer<T,V>::basic_clusterer(
   V                  radius,
   unsigned long      max_samples,
   const metric<T,V>& m)
 : metric_clusterer<T,V>(m),
   _radius(radius),
   _max_samples(max_samples)
{ }

/*
 * Copy constructor.
 * Create a clusterer with the same parameters.
 */
template <typename T, typename V>
basic_clusterer<T,V>::basic_clusterer(const basic_clusterer<T,V>& c)
 : metric_clusterer<T,V>(c),
   _radius(c._radius),
   _max_samples(c._max_samples)
{ }

/*
 * Destructor.
 */
template <typename T, typename V>
basic_clusterer<T,V>::~basic_clusterer() {
   /* do nothing */
}

/*
 * Clustering.
 * Return the cluster assignments and cluster centroids.
 */
template <typename T, typename V>
array<unsigned long> basic_clusterer<T,V>::cluster(
   const collection<T>&                 items,
   auto_collection< T, array_list<T> >& centroids) const 
{
   auto_collection< double, list<double> > weights(new list<double>());
   unsigned long n_items = items.size();
   for (unsigned long n = 0; n < n_items; n++) {
      auto_ptr<double> w(new double(1));
      weights->add(*w);
      w.release();
   }
   array<unsigned long> idx_map = math::random::util::randperm(n_items);
   return this->cluster(items, *weights, idx_map, centroids);
}

/*
 * Weighted clustering.
 * Return the cluster assignments and cluster centroids.
 */
template <typename T, typename V>
array<unsigned long> basic_clusterer<T,V>::cluster(
   const collection<T>&                 items,
   const collection<double>&            weights,
   auto_collection< T, array_list<T> >& centroids) const 
{
   array<unsigned long> idx_map = math::random::util::randperm(weights);
   return this->cluster(items, weights, idx_map, centroids);
}
  
/*
 * Weighted clustering using the specified random permutation.
 * Return the cluster assignments and cluster centroids.
 */
template <typename T, typename V>
array<unsigned long> basic_clusterer<T,V>::cluster(
   const collection<T>&                 items,
   const collection<double>&            weights,
   const array<unsigned long>&          idx_map,
   auto_collection< T, array_list<T> >& centroids) const 
{
   /* randomize item order using the specified random permutation */
   array_list<T>      items_array;
   array_list<double> weights_array;
   {
      array_list<T>      temp_items(items);
      array_list<double> temp_weights(weights);
      temp_items.subarray(idx_map, items_array);
      temp_weights.subarray(idx_map, weights_array);
   }
   /* allocate collection of centroids */
   centroids.reset(new array_list<T>());
   /* return now if there are no items to cluster */
   unsigned long n_items = items_array.size();
   if (n_items == 0)
      return array<unsigned long>();
   /* compute centroid of first cluster */
   {
      list<T>      clstr_items;
      list<double> clstr_weights;
      clstr_items.add(items_array[0]);
      clstr_weights.add(weights_array[0]);
      auto_ptr<T> cntrd = this->_metric.centroid(clstr_items, clstr_weights);
      centroids->add(*cntrd);
      cntrd.release();
   }
   /* allocate arrays for assignments and minimum item to cluster distance */
   array<unsigned long> assign(n_items);
   array<V> distance_min(n_items);
   /* allocate data structures for item to most recent cluster distance */
   array<unsigned long> changed_ids(1);   /* id of most recent cluster */
   array< array<V> > distances(n_items);  /* distance to most recent cluster */
   for (unsigned long n = 0; n < n_items; n++)
      distances[n].resize(1);
   /* iteratively update assignments and centroids */
   unsigned long n_samples = 1;
   bool added_cluster = true;
   do {
      /* compute distances to most recently added cluster */
      typename metric_clusterer<T,V>::distance_updater dist_updater(
         this->_metric,
         0,
         n_items - 1,
         items_array,
         *centroids,
         changed_ids,
         distances
      );
      dist_updater.run();
      /* update item assignments */
      for (unsigned long n = 0; n < n_items; n++) {
         V dist_recent = (distances[n])[0];
         if (dist_recent < distance_min[n]) {
            /* reassign to most recent cluster */
            distance_min[n] = dist_recent;
            assign[n] = 0;
         } else if (assign[n] == 0) {
            /* initialize minimum distance */
            distance_min[n] = dist_recent;
         }
      }
      /* check limit on number of samples */
      /* note: n_samples initialized to 1 */
      /*       so _max_samples set to 0   */
      /*       means unlimited samples    */
      if (n_samples == _max_samples)
         break;
      /* add a new cluster for the first (if any) item not covered */
      added_cluster = false;
      for (unsigned long n = 0; ((n < n_items) && (!added_cluster)); n++) {
         if (distance_min[n] > _radius) {
            /* add centroid for item to beginning of centroid array */
            list<T>      clstr_items;
            list<double> clstr_weights;
            clstr_items.add(items_array[n]);
            clstr_weights.add(weights_array[n]);
            auto_ptr<T> cntrd = this->_metric.centroid(
               clstr_items, 
               clstr_weights
            );
            centroids->prepend(*cntrd);
            cntrd.release();
            /* shift assignments */
            for (unsigned long nn = 0; nn < n_items; nn++)
               assign[nn]++;
            /* indicate that cluster has been added */
            n_samples++;
            added_cluster = true;
         }
      }
   } while (added_cluster);
   /* compute which clusters are nonempty */
   unsigned long K = centroids->size();
   array<bool> clstr_nonempty(K, false);
   for (unsigned long n = 0; n < n_items; n++)
      clstr_nonempty[assign[n]] = true;
   /* drop empty clusters */
   array<unsigned long> remap(K, K);   /* init remap to invalid cluster # K */
   unsigned long n_clusters = 0;
   for (unsigned long n = 0; n < K; n++) {
      auto_ptr<T> cntrd(&(centroids->remove_head()));
      if (clstr_nonempty[n]) {
         remap[n] = n_clusters;
         centroids->add(*cntrd);
         cntrd.release();
         n_clusters++;
      }
   }
   /* retrieve original assignment order */
   array<unsigned long> assignments(n_items);
   for (unsigned long n = 0; n < n_items; n++)
      assignments[idx_map[n]] = remap[assign[n]];
   return assignments;
}

} /* namespace covering */
} /* namespace clusterers */
} /* namespace clustering */
} /* namespace mlearning */

#endif
