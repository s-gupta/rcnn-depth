/*
 * Basic K-means clusterer.
 *
 * This template class implements a parallel K-means clustering algorithm that 
 * works with any valid distance metric.
 */
#ifndef MLEARNING__CLUSTERING__CLUSTERERS__KMEANS__BASIC_CLUSTERER_HH
#define MLEARNING__CLUSTERING__CLUSTERERS__KMEANS__BASIC_CLUSTERER_HH

#include "collections/abstract/collection.hh"
#include "collections/array_list.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "lang/array.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/random/util/randperm.hh"
#include "mlearning/clustering/clusterers/abstract/metric_clusterer.hh"
#include "mlearning/clustering/metrics/metric.hh"

namespace mlearning {
namespace clustering {
namespace clusterers {
namespace kmeans {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::array_list;
using collections::list;
using collections::pointers::auto_collection;
using lang::array;
using lang::exceptions::ex_invalid_argument;
using lang::pointers::auto_ptr;
using mlearning::clustering::clusterers::abstract::metric_clusterer;
using mlearning::clustering::metrics::metric;

/*
 * Base class containing code common to all K-means clusterers.
 */
class basic_clusterer_base {
public:
   /*
    * Destructor.
    */
   virtual ~basic_clusterer_base() = 0;

protected:
   /*
    * Compute ids of changed centroids given change flags.
    * Optionally specify the number of centroids to consider.
    */
   static array<unsigned long> compute_changed_ids(
      const array<bool>& /* change flags */
   );
   
   static array<unsigned long> compute_changed_ids(
      const array<bool>& /* change flags */, unsigned long /* # centroids */
   );
};

/*
 * Basic K-means clusterer.
 * T is the type of object being clustered.
 * V is the numeric type used to represent distance.
 */
template <typename T, typename V = double>
class basic_clusterer : public metric_clusterer<T,V>,
                        public basic_clusterer_base {
public:
   /*
    * Constructor.
    */
   explicit basic_clusterer(
      unsigned long,       /* K (desired number of clusters) */
      unsigned long,       /* maximum # of iterations (0 for unlimited) */
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
    * The algorithm iteratively assigns items to clusters (processing items 
    * in parallel) and then recomputes cluster centroids (processing clusters
    * in parallel), until the item assignments stabilize.
    *
    * Cluster centroids are initialized by solving a coarser K-means clustering
    * problem on a random subset of items sampled with probability proportional
    * to their weight.  Since initialization is randomized, different calls of
    * this method on identical arguments may produce different results.
    *
    * If weights are not specified, then all items are set to have weight one.
    *
    * Empty clusters are dropped from the final result.
    *
    * O(n + K + nK/p) time per iteration is required for clustering where n is
    * the number of items, K is the number of clusters, and p is the number of
    * available processors.
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
    * K-means parameters.
    */
   unsigned long _K;                /* desired number of clusters */
   unsigned long _max_iterations;   /* max # of iterations (0 for unlimited) */

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
 * Basic K-means clusterer implementation.
 ***************************************************************************/

/*
 * Constructor.
 */
template <typename T, typename V>
basic_clusterer<T,V>::basic_clusterer(
   unsigned long      K, 
   unsigned long      max_iterations,
   const metric<T,V>& m)
 : metric_clusterer<T,V>(m), 
   _K(K),
   _max_iterations(max_iterations)
{ } 

/*
 * Copy constructor.
 * Create a new clusterer with the same parameters.
 */
template <typename T, typename V>
basic_clusterer<T,V>::basic_clusterer(const basic_clusterer<T,V>& c)
 : metric_clusterer<T,V>(c),
   _K(c._K),
   _max_iterations(c._max_iterations)
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
   /* compute maximum number of clusters to return */
   unsigned long n_items = items_array.size();
   unsigned long K = (_K < n_items) ? _K : n_items;
   if ((n_items > 0) && (K == 0)) {
      throw ex_invalid_argument(
         "K must be nonzero during K-means clustering of nonempty collection"
      );
   }
   /* allocate collection of centroids */
   centroids.reset(new array_list<T>());
   /* return now if there are no items to cluster */
   if (n_items == 0)
      return array<unsigned long>();
   /* initialize assignments (to undefined cluster # K) */
   /* initialize cluster membership for centroids       */
   array<unsigned long>  assign(n_items, K);
   array< list<T> >      cluster_items(K);
   array< list<double> > cluster_weights(K);
   for (unsigned long n = 0; n < K; n++) {
      cluster_items[n].add(items_array[n]);
      cluster_weights[n].add(weights_array[n]);
   }
   /* initialize centroids */
   for (unsigned long n = 0; n < K; n++) {
      auto_ptr<T> cntrd = this->_metric.centroid(
         cluster_items[n],
         cluster_weights[n]
      );
      centroids->add(*cntrd);
      cntrd.release();
   }
   /* allocate arrays for distances */
   array< array<V> > distances(n_items);
   typename metric_clusterer<T,V>::distance_resizer dist_resizer(
      0, n_items - 1, K, distances
   );
   dist_resizer.run();
   /* initially mark the undefined cluster # K as changed */
   array<bool> has_changed(K+1, false);
   has_changed[K] = true;
   array<unsigned long> changed_ids;
   /* compute sizes of coarse problems */
   array<unsigned long> coarse_sizes =
      metric_clusterer<T,V>::coarse_problem_sizes(n_items, K);
   unsigned long n_problems = coarse_sizes.size();
   /* solve coarse to fine clustering problems */
   unsigned long n_items_prev = 0;
   for (unsigned long prob = 0; prob < n_problems; prob++) {
      /* get problem size */
      unsigned long n_items_curr = coarse_sizes[prob];
      /* compute distances between new items and unchanged clusters */
      {
         array<bool> has_not_changed(K);
         for (unsigned long n = 0; n < K; n++)
            has_not_changed[n] = !(has_changed[n]);
         array<unsigned long> unchanged_ids = 
            basic_clusterer<T,V>::compute_changed_ids(has_not_changed);
         typename metric_clusterer<T,V>::distance_updater dist_updater(
            this->_metric,
            n_items_prev,
            n_items_curr - 1,
            items_array,
            *centroids,
            unchanged_ids,
            distances
         );
         dist_updater.run();
      }
      /* initialize any empty clusters with new items */
      for (unsigned long n = 0, next_item = n_items_prev;
           ((n < K) && (next_item < n_items_curr)); n++)
      {
         if (cluster_items[n].is_empty()) {
            /* add item to cluster */
            cluster_items[n].add(items_array[next_item]);
            cluster_weights[n].add(weights_array[next_item]);
            /* compute centroid */
            auto_ptr<T> cntrd = this->_metric.centroid(
               cluster_items[n],
               cluster_weights[n]
            );
            T& cntrd_old = centroids->replace(n, *cntrd);
            cntrd.release();
            delete &cntrd_old;
            /* indicate that cluster has changed */
            has_changed[n] = true;
            next_item++;
         }
      }
      /* recompute changed ids to include any filled empty clusters */
      changed_ids = basic_clusterer<T,V>::compute_changed_ids(has_changed, K);
      /* iteratively update assignments and centroids */
      for (unsigned long n_iter = 0; 
           ((n_iter < _max_iterations) || (_max_iterations == 0));
           n_iter++)
      {
         /* store old assignments */
         array<unsigned long> assign_old = assign.subarray(0, n_items_curr-1);
         /* update distances */
         typename metric_clusterer<T,V>::distance_updater dist_updater(
            this->_metric,
            0,
            n_items_curr - 1,
            items_array,
            *centroids,
            changed_ids,
            distances
         );
         dist_updater.run();
         /* update assignments */
         typename metric_clusterer<T,V>::assignment_updater assign_updater(
            0,
            n_items_curr - 1,
            K,
            has_changed,
            changed_ids,
            distances,
            assign
         );
         assign_updater.run();
         /* compute which clusters have changed */
         for (unsigned long n = 0; n < K; n++)
            has_changed[n] = false;
         for (unsigned long n = 0; n < n_items_curr; n++) {
            unsigned long assign_id     = assign[n];
            unsigned long assign_id_old = assign_old[n];
            if (assign_id != assign_id_old) {
               has_changed[assign_id]     = true;
               has_changed[assign_id_old] = true;
            }
         }
         /* compute ids of changed clusters */
         changed_ids = basic_clusterer<T,V>::compute_changed_ids(
            has_changed, K
         );
         /* finish if no clusters changed */
         unsigned long n_changed = changed_ids.size();
         if (n_changed == 0)
            break;
         /* update cluster membership */
         for (unsigned long n = 0; n < K; n++) {
            cluster_items[n].clear();
            cluster_weights[n].clear();
         }
         for (unsigned long n = 0; n < n_items_curr; n++) {
            unsigned long assign_id = assign[n];
            cluster_items[assign_id].add(items_array[n]);
            cluster_weights[assign_id].add(weights_array[n]);
         }
         /* update centroids */
         typename metric_clusterer<T,V>::centroid_updater cntrd_updater(
            this->_metric,
            0,
            n_changed - 1,
            changed_ids,
            cluster_items,
            cluster_weights,
            centroids
         );   
         cntrd_updater.run();
      }
      /* update previous problem size */
      n_items_prev = n_items_curr;
   }
   /* drop empty clusters */
   array<unsigned long> remap(K, K);   /* init remap to invalid cluster # K */
   unsigned long n_clusters = 0;
   for (unsigned long n = 0; n < K; n++) {
      auto_ptr<T> cntrd(&(centroids->remove_head()));
      if (!(cluster_items[n].is_empty())) {
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

} /* namespace kmeans */
} /* namespace clusterers */
} /* namespace clustering */
} /* namespace mlearning */

#endif
