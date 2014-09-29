/*
 * Abstract metric clusterer.
 *
 * A metric clusterer is a centroid clusterer for which a distance metric on 
 * the space of items is used to determine the centroids.
 */
#ifndef MLEARNING__CLUSTERING__CLUSTERERS__ABSTRACT__METRIC_CLUSTERER_HH
#define MLEARNING__CLUSTERING__CLUSTERERS__ABSTRACT__METRIC_CLUSTERER_HH

#include "collections/abstract/collection.hh"
#include "collections/array_list.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "concurrent/threads/child_thread.hh"
#include "concurrent/threads/runnable.hh"
#include "concurrent/threads/thread.hh"
#include "lang/array.hh"
#include "lang/pointers/auto_ptr.hh"
#include "mlearning/clustering/clusterers/abstract/centroid_clusterer.hh"
#include "mlearning/clustering/metrics/metric.hh"

namespace mlearning {
namespace clustering {
namespace clusterers {
namespace abstract {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::array_list;
using collections::list;
using collections::pointers::auto_collection;
using concurrent::threads::child_thread;
using concurrent::threads::runnable;
using concurrent::threads::thread;
using lang::array;
using lang::pointers::auto_ptr;
using mlearning::clustering::metrics::metric;

/*
 * Abstract metric clusterer.
 * T is the type of object being clustered.
 * V is the numeric type used to represent distance.
 */
template <typename T, typename V = double>
class metric_clusterer : public centroid_clusterer<T> {
public:
   /*
    * Constructor.
    */
   explicit metric_clusterer(
      const metric<T,V>&         /* metric to use for clustering */
   );

   /*
    * Copy constructor.
    * Create a clusterer that uses the same metric.
    */
   metric_clusterer(const metric_clusterer<T,V>&);

   /*
    * Destructor.
    */
   virtual ~metric_clusterer() = 0;

   /*
    * Clustering.
    */
   using centroid_clusterer<T>::cluster;

   /*
    * Weighted clustering.
    * Return the cluster assignments and cluster centroids.
    */
   virtual array<unsigned long> cluster(
      const collection<T>&,                  /* items to cluster */
      const collection<double>&,             /* weight of each item */
      auto_collection< T, array_list<T> >&   /* returned cluster centroids */
   ) const = 0;
  
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
   using centroid_clusterer<T>::assign_cluster;

protected:
   /*
    * Clusterer parameters.
    */
   const metric<T,V>& _metric;   /* metric to use for clustering */

   /*
    * Runnable object for resizing distance arrays in parallel.
    */
   class distance_resizer : public runnable {
   public:
      /*
       * Constructor.
       */
      explicit distance_resizer(
         unsigned long,                /* first item to process */
         unsigned long,                /* last item to process */
         unsigned long,                /* new number of centroids */
         array< array<V> >&            /* item to centroid distances */
      );

      /*
       * Destructor.
       */
      virtual ~distance_resizer();

      /*
       * Resize the distance arrays.
       */
      virtual void run();

   protected:
      const unsigned long              _start;
      const unsigned long              _end;
      const unsigned long              _n_centroids;
      array< array<V> >&               _distances;
   };
   
   /*
    * Runnable object for updating item to centroid distances in parallel.
    */
   class distance_updater : public runnable {
   public:
      /*
       * Constructor.
       */
      explicit distance_updater(
         const metric<T,V>&,           /* metric */
         unsigned long,                /* first item to process */
         unsigned long,                /* last item to process */
         const array_list<T>&,         /* items */
         const array_list<T>&,         /* centroids */
         const array<unsigned long>&,  /* indices of changed centroids */
         array< array<V> >&            /* item to centroid distances */
      );

      /*
       * Destructor.
       */
      virtual ~distance_updater();

      /*
       * Update distances to changed centroids.
       */
      virtual void run();

   protected:
      const metric<T,V>&               _metric;
      const unsigned long              _start;
      const unsigned long              _end;
      const array_list<T>&             _items;
      const array_list<T>&             _centroids;
      const array<unsigned long>&      _changed_ids;
      array< array<V> >&               _distances;
   };
   
   /*
    * Runnable object for updating item assignments in parallel.
    */
   class assignment_updater : public runnable {
   public:
      /*
       * Constructor.
       */
      explicit assignment_updater(
         unsigned long,                /* first item to process */
         unsigned long,                /* last item to process */
         unsigned long,                /* number of centroids */
         const array<bool>&,           /* changed flag for each centroid */
         const array<unsigned long>&,  /* indices of changed centroids */
         const array< array<V> >&,     /* item to centroid distances */
         array<unsigned long>&         /* assignments */
      );

      /*
       * Destructor.
       */
      virtual ~assignment_updater();

      /*
       * Update assignments.
       * Each item is assigned to the cluster with the closest centroid.
       */
      virtual void run();

   protected:
      const unsigned long              _start;
      const unsigned long              _end;
      const unsigned long              _n_centroids;
      const array<bool>&               _has_changed;
      const array<unsigned long>&      _changed_ids;
      const array< array<V> >&         _distances;
      array<unsigned long>&            _assignments;
   };

   /*
    * Runnable object for updating cluster centroids in parallel.
    */
   class centroid_updater : public runnable {
   public:
      /*
       * Constructor.
       */
      explicit centroid_updater(
         const metric<T,V>&,                 /* metric */
         unsigned long,                      /* first changed cluster to process */
         unsigned long,                      /* last changed cluster to process */
         const array<unsigned long>&,        /* indices of changed clusters */
         const array< list<T> >&,            /* items in each cluster */
         const array< list<double> >&,       /* weight of each item */
         auto_collection<T, array_list<T> >& /* centroids */
      );

      /*
       * Destructor.
       */
      virtual ~centroid_updater();
      
      /*
       * Update centroids of changed clusters.
       */
      virtual void run();

   protected:
      const metric<T,V>&                     _metric;
      const unsigned long                    _start;
      const unsigned long                    _end;
      const array<unsigned long>&            _changed_ids;
      const array< list<T> >&                _cluster_items;
      const array< list<double> >&           _cluster_weights;
      auto_collection<T, array_list<T> >&    _centroids;
   };
};

/***************************************************************************
 * Metric clusterer implementation.
 ***************************************************************************/

/*
 * Constructor.
 */
template <typename T, typename V>
metric_clusterer<T,V>::metric_clusterer(const metric<T,V>& m)
 : _metric(m)
{ }

/*
 * Copy constructor.
 * Create a clusterer that uses the same metric.
 */
template <typename T, typename V>
metric_clusterer<T,V>::metric_clusterer(const metric_clusterer<T,V>& c)
 : _metric(c._metric)
{ }

/*
 * Pure virtual destructor.
 */
template <typename T, typename V>
metric_clusterer<T,V>::~metric_clusterer() { }

/*
 * Cluster assignment.
 * Return the id of the cluster centroid closest to the item.
 * Note that the centroid array must not be empty.
 */
template <typename T, typename V>
unsigned long metric_clusterer<T,V>::assign_cluster(
   const T&             item, 
   const array_list<T>& centroids) const
{
   /* initialize id and distance */
   unsigned long id = 0;
   V min_dist = _metric.distance(item, centroids[0]);
   /* find closest centroid */
   unsigned long n_items = centroids.size();
   for (unsigned long n = 1; n < n_items; n++) {
      V dist = _metric.distance(item, centroids[n]);
      if (dist < min_dist) {
         min_dist = dist;
         id = n;
      }
   }
   return id;
}

/***************************************************************************
 * Distance resizer implementation.
 ***************************************************************************/

/*
 * Constructor.
 */
template <typename T, typename V>
metric_clusterer<T,V>::distance_resizer::distance_resizer(
   unsigned long      start,
   unsigned long      end,
   unsigned long      n_centroids,
   array< array<V> >& distances)
 : _start(start),
   _end(end),
   _n_centroids(n_centroids),
   _distances(distances)
{ }

/*
 * Destructor.
 */
template <typename T, typename V>
metric_clusterer<T,V>::distance_resizer::~distance_resizer() {
   /* do nothing */
}

/*
 * Resize the distance arrays.
 */
template <typename T, typename V>
void metric_clusterer<T,V>::distance_resizer::run() {
   if ((thread::processors() > 1) && (_start < _end)) {
      /* recursively resize distances in parallel */
      unsigned long mid = (_start + _end)/2;
      distance_resizer dist_resizer_left(
         _start, mid, _n_centroids, _distances
      );
      distance_resizer dist_resizer_right(
         mid+1, _end, _n_centroids, _distances
      );
      child_thread::start(dist_resizer_left, dist_resizer_right);
   } else {
      /* resize distances sequentially */
      for (unsigned long n = _start; n <= _end; n++)
         _distances[n].resize(_n_centroids);
   }
}

/***************************************************************************
 * Distance updater implementation.
 ***************************************************************************/

/*
 * Constructor.
 */
template <typename T, typename V>
metric_clusterer<T,V>::distance_updater::distance_updater(
   const metric<T,V>&          m,
   unsigned long               start,
   unsigned long               end,
   const array_list<T>&        items,
   const array_list<T>&        centroids,
   const array<unsigned long>& changed_ids,
   array< array<V> >&          distances)
 : _metric(m),
   _start(start),
   _end(end),
   _items(items),
   _centroids(centroids),
   _changed_ids(changed_ids),
   _distances(distances)
{ }

/*
 * Destructor.
 */
template <typename T, typename V>
metric_clusterer<T,V>::distance_updater::~distance_updater() {
   /* do nothing */
}

/*
 * Update distances to changed centroids.
 */
template <typename T, typename V>
void metric_clusterer<T,V>::distance_updater::run() {
   if ((thread::processors() > 1) && (_start < _end)) {
      /* recursively update distances in parallel */
      unsigned long mid = (_start + _end)/2;
      distance_updater dist_updater_left(
         _metric, _start, mid, _items, _centroids, _changed_ids, _distances
      );
      distance_updater dist_updater_right(
         _metric, mid+1, _end, _items, _centroids, _changed_ids, _distances
      );
      child_thread::start(dist_updater_left, dist_updater_right);
   } else {
      /* update distances sequentially */
      unsigned long n_changed = _changed_ids.size();
      for (unsigned long n = _start; n <= _end; n++) {
         /* get item and distance array */
         const T& item      = _items[n];
         array<V>& distance = _distances[n];
         /* update distance to changed clusters */
         for (unsigned long n_id = 0; n_id < n_changed; n_id++) {
            unsigned long id = _changed_ids[n_id];
            distance[id] = _metric.distance(item, _centroids[id]);
         }
      }
   }
}

/***************************************************************************
 * Assignment updater implementation.
 ***************************************************************************/

/*
 * Constructor.
 */
template <typename T, typename V>
metric_clusterer<T,V>::assignment_updater::assignment_updater(
   unsigned long               start,
   unsigned long               end,
   unsigned long               n_centroids,
   const array<bool>&          has_changed,
   const array<unsigned long>& changed_ids,
   const array< array<V> >&    distances,
   array<unsigned long>&       assignments)
 : _start(start),
   _end(end),
   _n_centroids(n_centroids),
   _has_changed(has_changed),
   _changed_ids(changed_ids),
   _distances(distances),
   _assignments(assignments)
{ }

/*
 * Destructor.
 */
template <typename T, typename V>
metric_clusterer<T,V>::assignment_updater::~assignment_updater() {
   /* do nothing */
}

/*
 * Update assignments.
 * Each item is assigned to the cluster with the closest centroid.
 */
template <typename T, typename V>
void metric_clusterer<T,V>::assignment_updater::run() {
   if ((thread::processors() > 1) && (_start < _end)) {
      /* recursively update assignments in parallel */
      unsigned long mid = (_start + _end)/2;
      assignment_updater assign_updater_left(
         _start, mid, _n_centroids, _has_changed,
         _changed_ids, _distances, _assignments
      );
      assignment_updater assign_updater_right(
         mid+1, _end, _n_centroids, _has_changed,
         _changed_ids, _distances, _assignments
      );
      child_thread::start(assign_updater_left, assign_updater_right);
   } else {
      /* update assignments sequentially */
      unsigned long n_changed = _changed_ids.size();
      for (unsigned long n = _start; n <= _end; n++) {
         /* get distance array and current assignment */
         const array<V>& distance = _distances[n];
         unsigned long assign_id  = _assignments[n];
         /* check if cluster to which item is assigned has changed */
         if (_has_changed[assign_id]) {
            /* search all distances to find minimum */
            assign_id = 0;
            V min_dist = distance[0];
            for (unsigned long id = 1; id < _n_centroids; id++) {
               if (distance[id] < min_dist) {
                  min_dist = distance[id];
                  assign_id = id;
               }
            }
         } else {
            /* search only distances to current and changed clusters */
            V min_dist = distance[assign_id];
            for (unsigned long n_id = 0; n_id < n_changed; n_id++) {
               unsigned long id = _changed_ids[n_id];
               if (distance[id] < min_dist) {
                  min_dist = distance[id];
                  assign_id = id;
               }
            }
         }
         _assignments[n] = assign_id;
      }
   }
}

/***************************************************************************
 * Centroid updater implementation.
 ***************************************************************************/

/*
 * Constructor.
 */
template <typename T, typename V>
metric_clusterer<T,V>::centroid_updater::centroid_updater(
   const metric<T,V>&                  m,
   unsigned long                       start,
   unsigned long                       end,
   const array<unsigned long>&         changed_ids,
   const array< list<T> >&             cluster_items,
   const array< list<double> >&        cluster_weights,
   auto_collection<T, array_list<T> >& centroids)
 : _metric(m),
   _start(start),
   _end(end),
   _changed_ids(changed_ids),
   _cluster_items(cluster_items),
   _cluster_weights(cluster_weights),
   _centroids(centroids)
{ }

/*
 * Destructor.
 */
template <typename T, typename V>
metric_clusterer<T,V>::centroid_updater::~centroid_updater() {
   /* do nothing */
}

/*
 * Update centroids of changed clusters.
 */
template <typename T, typename V>
void metric_clusterer<T,V>::centroid_updater::run() {
   if ((thread::processors() > 1) && (_start < _end)) {
      /* recursively update centroids in parallel */
      unsigned long mid = (_start + _end)/2;
      centroid_updater cntrd_updater_left(
         _metric, _start, mid, _changed_ids,
         _cluster_items, _cluster_weights, _centroids
      );
      centroid_updater cntrd_updater_right(
         _metric, mid+1, _end, _changed_ids,
         _cluster_items, _cluster_weights, _centroids
      );
      child_thread::start(cntrd_updater_left, cntrd_updater_right);
   } else {
      /* update centroids sequentially */
      for (unsigned long n = _start; n <= _end; n++) {
         unsigned long id = _changed_ids[n];
         if (!(_cluster_items[id].is_empty())) {
            auto_ptr<T> cntrd = _metric.centroid(
               _cluster_items[id],
               _cluster_weights[id]
            );
            T& cntrd_old = _centroids->replace(id, *cntrd);
            cntrd.release();
            delete &cntrd_old;
         }
      }
   }
}

} /* namespace abstract */
} /* namesapce clusterers */
} /* namespace clustering */
} /* namespace mlearning */

#endif
