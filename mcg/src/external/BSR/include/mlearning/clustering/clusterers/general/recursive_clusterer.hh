/*
 * General purpose recursive clusterer.
 *
 * This clusterer produces an approximate clustering by recursively merging 
 * the results of clustering a fixed maximum number of items.  The original 
 * dataset is recursively split into subsets containing no more than a 
 * specified maximum number of items (the base problem size).  A user-specified
 * clusterer is applied to these subsets to produce centroids.  Results are 
 * merged by applying the same clusterer to the centroids, weighted according 
 * to cluster membership.
 *
 * The independent clustering subproblems are solved in parallel when possible.
 *
 * Note that the maximum number of items in any clustering subproblem solved 
 * in this recursive process is max(B, 2K), where B is the user-specified 
 * base problem size, and K is the maximum number of clusters the user-
 * specified clusterer can return.
 */
#ifndef MLEARNING__CLUSTERING__CLUSTERERS__GENERAL__RECURSIVE_CLUSTERER_HH
#define MLEARNING__CLUSTERING__CLUSTERERS__GENERAL__RECURSIVE_CLUSTERER_HH

#include "collections/abstract/collection.hh"
#include "collections/array_list.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "concurrent/threads/child_thread.hh"
#include "concurrent/threads/runnable.hh"
#include "lang/array.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"
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
using collections::list;
using collections::pointers::auto_collection;
using concurrent::threads::child_thread;
using concurrent::threads::runnable;
using lang::array;
using lang::exceptions::ex_invalid_argument;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;
using mlearning::clustering::clusterers::abstract::centroid_clusterer;

/*
 * Recursive clusterer.
 */
template <typename T>
class recursive_clusterer : public centroid_clusterer<T> {
public:
   /*
    * Constructor.
    */
   explicit recursive_clusterer(
      const centroid_clusterer<T>&, /* underlying clusterer to employ */
      unsigned long                 /* maximum # items in base problem */
   );

   /*
    * Copy constructor.
    * Create a recursive cluster with the same parameters.
    */
   recursive_clusterer(const recursive_clusterer<T>&);

   /*
    * Destructor.
    */
   virtual ~recursive_clusterer();
  
   /*
    * Clustering.
    */
   using centroid_clusterer<T>::cluster;
   
   /*
    * Weighted clustering.
    * Return the cluster assignments and cluster centroids.
    */
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
   /************************************************************************
    * Recursive clusterer paremeters.
    ************************************************************************/

   const centroid_clusterer<T>& _clusterer;
   unsigned long                _max_problem_size;

   /************************************************************************
    * Recursive clusterer helpers.
    ************************************************************************/

   /*
    * Runnable object for recursively reducing the size of a clustering 
    * problem by merging the results of clustering subsets.  The number 
    * of cluster centroids returned is at most max(B, 2K), where B is the 
    * maximum problem size specified by the recursive clusterer and K is 
    * the maximum number of clusters the underlying clusterer can return.
    *
    * Results of the reducing clusterer are useful as input into a final 
    * application of the underlying clusterer.
    */
   class reducing_clusterer : public runnable {
   public:
      /*
       * Constructor.
       * Note that the order of the items should be randomized prior to 
       * calling this constructor.
       */
      explicit reducing_clusterer(
         const recursive_clusterer<T>&,                  /* recursive clusterer */
         const collection<T>&,                           /* items to cluster (in random order) */
         const collection<double>&,                      /* weight of each item */
         auto_collection< T, array_list<T> >&,           /* returned cluster centroids */
         auto_collection< double, array_list<double> >&  /* returned centroid weights */
      );

      /*
       * Destructor.
       */
      virtual ~reducing_clusterer();

      /*
       * Run the reducing clusterer.
       */
      virtual void run();

   protected:
      const recursive_clusterer<T>&                      _r_clusterer;
      const collection<T>&                               _items;
      const collection<double>&                          _weights;
      auto_collection< T, array_list<T> >&               _centroids;
      auto_collection< double, array_list<double> >&     _centroid_weights;
   };

   /*
    * Runnable object for applying a centroid clusterer (and retrieving 
    * centroids and centroid weights, while disregarding assignments).
    * This runnable object wraps the cluster(...) method.
    */
   class runnable_clusterer : public runnable {
   public:
      /*
       * Constructor.
       */
      explicit runnable_clusterer(
         const centroid_clusterer<T>&,                   /* centroid clusterer */
         const collection<T>&,                           /* items to cluster */
         const collection<double>&,                      /* weight of each item */
         auto_collection< T, array_list<T> >&,           /* returned cluster centroids */
         auto_collection< double, array_list<double> >&  /* returned centroid weights */
      );

      /*
       * Destructor.
       */
      virtual ~runnable_clusterer();

      /*
       * Run the clusterer.
       */
      virtual void run();

   protected:
      const centroid_clusterer<T>&                       _clusterer;
      const collection<T>&                               _items;
      const collection<double>&                          _weights;
      auto_collection< T, array_list<T> >&               _centroids;
      auto_collection< double, array_list<double> >&     _centroid_weights;
   };
};

/***************************************************************************
 * Recursive clusterer implementation.
 ***************************************************************************/

/*
 * Constructor.
 */
template <typename T>
recursive_clusterer<T>::recursive_clusterer(
   const centroid_clusterer<T>& c,
   unsigned long                n)
 : _clusterer(c),
   _max_problem_size(n)
{ }

/*
 * Copy constructor.
 * Create a recursive cluster with the same parameters.
 */
template <typename T>
recursive_clusterer<T>::recursive_clusterer(const recursive_clusterer<T>& c)
 : _clusterer(c._clusterer),
   _max_problem_size(c._max_problem_size)
{ }

/*
 * Destructor.
 */
template <typename T>
recursive_clusterer<T>::~recursive_clusterer() {
   /* do nothing */
}
   
/*
 * Weighted clustering.
 * Return the cluster assignments and cluster centroids.
 */
template <typename T>
array<unsigned long> recursive_clusterer<T>::cluster(
   const collection<T>&                 items,
   const collection<double>&            weights,
   auto_collection< T, array_list<T> >& centroids) const
{
   /* check problem size */
   unsigned long n_items = items.size();
   if (n_items <= _max_problem_size) {
      /* cluster small problem */
      return _clusterer.cluster(items, weights, centroids);
   } else if (_max_problem_size == 0) {
      /* cannot enforce a max problem size of zero */
      throw ex_invalid_argument(
         "max_problem_size must be nonzero when clustering nonempty collection"
      );
   } else {
      /* randomize item order (and put weights in same order as items) */
      array_list<T> items_array(items);
      array_list<double> weights_array;
      {
         array<unsigned long> idx_map = items_array.randperm();
         array_list<double> temp(weights);
         temp.subarray(idx_map, weights_array);
      }
      /* perform recursive clustering to reduce problem size */
      auto_collection< T, array_list<T> >           cntrds;
      auto_collection< double, array_list<double> > cntrd_weights;
      reducing_clusterer c_reduce(
         *this,
         items_array,
         weights_array,
         cntrds,
         cntrd_weights
      );
      c_reduce.run();
      /* apply final clustering step */
      _clusterer.cluster(*cntrds, *cntrd_weights, centroids);
      /* compute assignments to the final centroids */
      return _clusterer.assign_cluster(items, *centroids);
   }
}

/*
 * Cluster assignment.
 * Return the id of the cluster centroid closest to the item.
 * Note that the centroid array must not be empty.
 */
template <typename T>
unsigned long recursive_clusterer<T>::assign_cluster(
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
array<unsigned long> recursive_clusterer<T>::assign_cluster(
   const collection<T>& items,
   const array_list<T>& centroids) const
{
   return _clusterer.assign_cluster(items, centroids);
}

/***************************************************************************
 * Reducing clusterer implementation.
 ***************************************************************************/

/*
 * Constructor.
 * Note that the order of the items should be randomized prior to 
 * calling this constructor.
 */
template <typename T>
recursive_clusterer<T>::reducing_clusterer::reducing_clusterer(
   const recursive_clusterer<T>&                  r_clusterer,
   const collection<T>&                           items,
   const collection<double>&                      weights,
   auto_collection< T, array_list<T> >&           centroids,
   auto_collection< double, array_list<double> >& centroid_weights)
 : _r_clusterer(r_clusterer),
   _items(items),
   _weights(weights),
   _centroids(centroids),
   _centroid_weights(centroid_weights)
{ }

/*
 * Destructor.
 */
template <typename T>
recursive_clusterer<T>::reducing_clusterer::~reducing_clusterer() {
   /* do nothing */
}

/*
 * Run the reducing clusterer.
 */
template <typename T>
void recursive_clusterer<T>::reducing_clusterer::run() {
   /* split items into subsets */
   unsigned long n_items       = _items.size();
   unsigned long n_items_left  = n_items/2;
   unsigned long n_items_right = n_items - n_items_left;
   list<T> items_left;
   list<T> items_right;
   list<double> weights_left;
   list<double> weights_right;
   auto_ptr< iterator<T> >      i_items   = _items.iter_create();
   auto_ptr< iterator<double> > i_weights = _weights.iter_create();
   for (unsigned long n = 0; n < n_items_left; n++) {
      items_left.add(i_items->next());
      weights_left.add(i_weights->next());
   }
   for (unsigned long n = 0; n < n_items_right; n++) {
      items_right.add(i_items->next());
      weights_right.add(i_weights->next());
   }
   /* create containers for results of clustering each subset */
   auto_collection< T, array_list<T> > cntrds_left;
   auto_collection< T, array_list<T> > cntrds_right;
   auto_collection< double, array_list<double> > cntrd_weights_left;
   auto_collection< double, array_list<double> > cntrd_weights_right;
   /* perform clustering */
   if ((n_items_left  <= _r_clusterer._max_problem_size) &&
       (n_items_right <= _r_clusterer._max_problem_size))
   {
      /* cluster each subset */
      runnable_clusterer c_run_left(
         _r_clusterer._clusterer,
         items_left,
         weights_left,
         cntrds_left,
         cntrd_weights_left
      );
      runnable_clusterer c_run_right(
         _r_clusterer._clusterer,
         items_right,
         weights_right,
         cntrds_right,
         cntrd_weights_right
      );
      child_thread::run(c_run_left, c_run_right);
   } else {
      /* reduce the size of each subset */
      reducing_clusterer c_reduce_left( 
         _r_clusterer,
         items_left,
         weights_left,
         cntrds_left,
         cntrd_weights_left
      );
      reducing_clusterer c_reduce_right(
         _r_clusterer,
         items_right,
         weights_right,
         cntrds_right,
         cntrd_weights_right
      );
      child_thread::run(c_reduce_left, c_reduce_right);
      /* cluster the centroids if merged set will be too large */
      unsigned long n_cntrds_left  = cntrds_left->size();
      unsigned long n_cntrds_right = cntrds_right->size();
      if ((n_cntrds_left + n_cntrds_right) > _r_clusterer._max_problem_size) {
         /* create containers for results of clustering centroids */
         auto_collection< T, array_list<T> > cntrd_cntrds_left;
         auto_collection< T, array_list<T> > cntrd_cntrds_right;
         auto_collection< double, array_list<double> > cntrd_cntrd_weights_left;
         auto_collection< double, array_list<double> > cntrd_cntrd_weights_right;
         /* cluster centroids */
         runnable_clusterer c_run_left(
            _r_clusterer._clusterer,
            *cntrds_left,
            *cntrd_weights_left,
            cntrd_cntrds_left,
            cntrd_cntrd_weights_left
         );
         runnable_clusterer c_run_right(
            _r_clusterer._clusterer,
            *cntrds_right,
            *cntrd_weights_right,
            cntrd_cntrds_right,
            cntrd_cntrd_weights_right
         );
         child_thread::run(c_run_left, c_run_right);
         /* update centroids and centroid weights */
         cntrds_left  = cntrd_cntrds_left;
         cntrds_right = cntrd_cntrds_right;
         cntrd_weights_left  = cntrd_cntrd_weights_left;
         cntrd_weights_right = cntrd_cntrd_weights_right;
      }
   }
   /* place centroids into a single array */
   auto_ptr< array_list<T> > ptr_cntrds_left;
   auto_ptr< array_list<T> > ptr_cntrds_right;
   _centroids.reset(new array_list<T>());
   _centroids->add(*cntrds_left);
   ptr_cntrds_left.reset(cntrds_left.release());
   _centroids->add(*cntrds_right);
   ptr_cntrds_right.reset(cntrds_right.release());
   /* place centroid weights into a single array */
   auto_ptr< array_list<double> > ptr_cntrd_weights_left;
   auto_ptr< array_list<double> > ptr_cntrd_weights_right;
   _centroid_weights.reset(new array_list<double>());
   _centroid_weights->add(*cntrd_weights_left);
   ptr_cntrd_weights_left.reset(cntrd_weights_left.release());
   _centroid_weights->add(*cntrd_weights_right);
   ptr_cntrd_weights_right.reset(cntrd_weights_right.release());
}

/***************************************************************************
 * Runnable clusterer implementation.
 ***************************************************************************/

/*
 * Constructor.
 */
template <typename T>
recursive_clusterer<T>::runnable_clusterer::runnable_clusterer(
   const centroid_clusterer<T>&                   c,
   const collection<T>&                           items,
   const collection<double>&                      weights,
   auto_collection< T, array_list<T> >&           centroids,
   auto_collection< double, array_list<double> >& centroid_weights)
 : _clusterer(c),
   _items(items),
   _weights(weights),
   _centroids(centroids),
   _centroid_weights(centroid_weights)
{ }

/*
 * Destructor.
 */
template <typename T>
recursive_clusterer<T>::runnable_clusterer::~runnable_clusterer() {
   /* do nothing */
}

/*
 * Run the clusterer.
 */
template <typename T>
void recursive_clusterer<T>::runnable_clusterer::run() {
   /* compute centroids and assignments */
   array<unsigned long> assignments = _clusterer.cluster(
      _items,
      _weights,
      _centroids
   );
   /* compute weights from assignments */
   unsigned long n_items = assignments.size();
   unsigned long n_clusters = _centroids->size();
   _centroid_weights.reset(new array_list<double>());
   for (unsigned long n = 0; n < n_clusters; n++) {
      auto_ptr<double> w(new double(0));
      _centroid_weights->add(*w);
      w.release();
   }
   auto_ptr< iterator<double> > i = _weights.iter_create();
   for (unsigned long n = 0; n < n_items; n++)
      (*_centroid_weights)[assignments[n]] += i->next();
}

} /* namespace general */
} /* namespace clusterers */
} /* namespace clustering */
} /* namespace mlearning */

#endif
