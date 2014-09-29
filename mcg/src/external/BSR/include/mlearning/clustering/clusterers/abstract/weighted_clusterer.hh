/*
 * Abstract weighted clusterer.
 *
 * A weighted cluster takes into account real-valued weights for data items 
 * when computing their assignment to clusters.
 */
#ifndef MLEARNING__CLUSTERING__CLUSTERERS__ABSTRACT__WEIGHTED_CLUSTERER_HH
#define MLEARNING__CLUSTERING__CLUSTERERS__ABSTRACT__WEIGHTED_CLUSTERER_HH

#include "collections/abstract/collection.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "lang/array.hh"
#include "lang/pointers/auto_ptr.hh"
#include "mlearning/clustering/clusterers/abstract/clusterer.hh"

namespace mlearning {
namespace clustering {
namespace clusterers {
namespace abstract {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::list;
using collections::pointers::auto_collection;
using lang::array;
using lang::pointers::auto_ptr;

/*
 * Base class containing code common to all weighted_clusterer<T> templates.
 */
class weighted_clusterer_base {
public:
   /*
    * Destructor.
    */
   virtual ~weighted_clusterer_base() = 0;

   /*
    * Compute the total weight assigned to each cluster (the sum 
    * of the weights of its members) given the item weights and 
    * assignment map.
    */
   static array<double> cluster_weights(
      const collection<double>&,    /* weight of each item */
      const array<unsigned long>&   /* assignments */
   );
};

/*
 * Abstract interface for weighted clusterers.
 */
template <typename T>
class weighted_clusterer : public clusterer<T>,
                           public weighted_clusterer_base {
public:
   /*
    * Destructor.
    */
   virtual ~weighted_clusterer() = 0;
   
   /*
    * Clustering.
    *
    * Return the cluster assignments of the items in the 
    * order in which they appear in the given collection.
    * 
    * The default implementation calls the weighted version
    * of the cluster(...) method with all weights set to one.
    */
   virtual array<unsigned long> cluster(
      const collection<T>&          /* items to cluster */
   ) const;

   /*
    * Weighted clustering.
    */
   virtual array<unsigned long> cluster(
      const collection<T>&,         /* items to cluster */
      const collection<double>&     /* weight of each item */
   ) const = 0;
};

/*
 * Pure virtual destructor.
 */
template <typename T>
weighted_clusterer<T>::~weighted_clusterer() { }

/*
 * Clustering.
 * Call the weighted version with all weights set to one.
 */
template <typename T>
array<unsigned long> weighted_clusterer<T>::cluster(
   const collection<T>& items) const
{
   auto_collection< double, list<double> > weights(new list<double>());
   unsigned long n_items = items.size();
   for (unsigned long n = 0; n < n_items; n++) {
      auto_ptr<double> w(new double(1));
      weights->add(*w);
      w.release();
   }
   return this->cluster(items, *weights);
}

} /* namespace abstract */
} /* namesapce clusterers */
} /* namespace clustering */
} /* namespace mlearning */

#endif
