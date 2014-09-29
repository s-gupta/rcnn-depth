/*
 * Abstract centroid clusterer.
 *
 * This abstract template class provides an interface common to all clustering
 * algorithms which define clusters by centroids produced by a weighted 
 * clustering process.
 */
#ifndef MLEARNING__CLUSTERING__CLUSTERERS__ABSTRACT__CENTROID_CLUSTERER_HH
#define MLEARNING__CLUSTERING__CLUSTERERS__ABSTRACT__CENTROID_CLUSTERER_HH

#include "collections/abstract/collection.hh"
#include "collections/array_list.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "concurrent/threads/child_thread.hh"
#include "concurrent/threads/runnable.hh"
#include "concurrent/threads/thread.hh"
#include "lang/array.hh"
#include "lang/pointers/auto_ptr.hh"
#include "mlearning/clustering/clusterers/abstract/weighted_clusterer.hh"

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

/*
 * Abstract centroid clusterer.
 */
template <typename T>
class centroid_clusterer : public weighted_clusterer<T> {
public:
   /*
    * Destructor.
    */
   virtual ~centroid_clusterer() = 0;
   
   /*
    * Clustering.
    *
    * Return the cluster assignments of the items in the
    * order in which they appear in the given collection.
    *
    * Optionally return the cluster centroids.
    *
    * The default implementation calls the weighted version
    * of the cluster(...) method with all weights set to one.
    */
   virtual array<unsigned long> cluster(
      const collection<T>&                   /* items to cluster */
   ) const;

   virtual array<unsigned long> cluster(
      const collection<T>&,                  /* items to cluster */
      auto_collection< T, array_list<T> >&   /* returned cluster centroids */
   ) const;

   /*
    * Weighted clustering.
    * Return the cluster assignments and (optionally) the centroids.
    */
   virtual array<unsigned long> cluster(
      const collection<T>&,                  /* items to cluster */
      const collection<double>&              /* weight of each item */
   ) const;
   
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
   virtual unsigned long assign_cluster(
      const T&,                              /* item */
      const array_list<T>&                   /* centroids */
   ) const = 0;

   /*
    * Cluster assignment.
    * For each item, return the id of the closest cluster centroid.
    * Note that the centroid array must not be empty.
    */
   virtual array<unsigned long> assign_cluster(
      const collection<T>&,                  /* items */
      const array_list<T>&                   /* centroids */
   ) const;
      
protected:
   /*
    * Runnable object for assigning items to clusters in parallel.
    */
   class cluster_assigner : public runnable {
   public:
      /*
       * Constructor.
       */
      explicit cluster_assigner(
         const centroid_clusterer<T>&, /* clusterer */
         unsigned long,                /* first item to process */
         unsigned long,                /* last item to process  */
         const array_list<T>&,         /* items */
         const array_list<T>&,         /* centroids */
         array<unsigned long>&         /* assignments */
      );

      /*
       * Destructor.
       */
      virtual ~cluster_assigner();
      
      /*
       * Compute assignments.
       */
      virtual void run();

   protected:
      const centroid_clusterer<T>&     _clusterer;
      const unsigned long              _start;
      const unsigned long              _end;
      const array_list<T>&             _items;
      const array_list<T>&             _centroids;
      array<unsigned long>&            _assignments;
   };
};

/***************************************************************************
 * Centroid clusterer implementation.
 ***************************************************************************/

/*
 * Pure virtual destructor.
 */
template <typename T>
centroid_clusterer<T>::~centroid_clusterer() { }

/*
 * Clustering (ignoring returned centroids).
 * Call the version that returns the centroids.
 */
template <typename T>
array<unsigned long> centroid_clusterer<T>::cluster(
   const collection<T>& items) const 
{
   auto_collection< T, array_list<T> > centroids;
   return this->cluster(items, centroids);
}

/*
 * Clustering.
 * Call the weighted version with all weights set to one.
 */
template <typename T>
array<unsigned long> centroid_clusterer<T>::cluster(
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
   return this->cluster(items, *weights, centroids);
}

/*
 * Weighted clustering (ignoring returned centroids).
 * Call the version that returns the centroids.
 */
template <typename T>
array<unsigned long> centroid_clusterer<T>::cluster(
   const collection<T>&      items,
   const collection<double>& weights) const 
{
   auto_collection< T, array_list<T> > centroids;
   return this->cluster(items, weights, centroids);
}

/*
 * Cluster assignment.
 * For each item, return the id of the closest cluster centroid.
 * Note that the centroid array must not be empty.
 */
template <typename T>
array<unsigned long> centroid_clusterer<T>::assign_cluster(
   const collection<T>& items,
   const array_list<T>& centroids) const
{
   array_list<T> items_array(items);
   unsigned long n_items = items_array.size();
   array<unsigned long> assignments(n_items);
   if (n_items > 0) {
      cluster_assigner clstr_assigner(
         *this, 0, n_items - 1, items_array, centroids, assignments
      );
      clstr_assigner.run();
   }
   return assignments;
}

/***************************************************************************
 * Cluster assigner implementation.
 ***************************************************************************/

/*
 * Constructor.
 */
template <typename T>
centroid_clusterer<T>::cluster_assigner::cluster_assigner(
   const centroid_clusterer<T>& c,
   unsigned long                start,
   unsigned long                end,
   const array_list<T>&         items,
   const array_list<T>&         centroids,
   array<unsigned long>&        assignments)
 : _clusterer(c),
   _start(start),
   _end(end),
   _items(items),
   _centroids(centroids),
   _assignments(assignments)
{ }

/*
 * Destructor.
 */
template <typename T>
centroid_clusterer<T>::cluster_assigner::~cluster_assigner() {
   /* do nothing */
}

/*
 * Compute assignments.
 */
template <typename T>
void centroid_clusterer<T>::cluster_assigner::run() {
   if ((thread::processors() > 1) && (_start < _end)) {
      /* recursively compute assignments in parallel */
      unsigned long mid = (_start + _end)/2;
      cluster_assigner clstr_assigner_left(
         _clusterer, _start,  mid,  _items, _centroids, _assignments
      );
      cluster_assigner clstr_assigner_right(
         _clusterer, mid + 1, _end, _items, _centroids, _assignments
      );
      child_thread::start(clstr_assigner_left, clstr_assigner_right);
   } else {
      /* compute assignments sequentially */
      for (unsigned long n = _start; n <= _end; n++)
         _assignments[n] = _clusterer.assign_cluster(_items[n], _centroids);
   }
}

} /* namespace abstract */
} /* namesapce clusterers */
} /* namespace clustering */
} /* namespace mlearning */

#endif
