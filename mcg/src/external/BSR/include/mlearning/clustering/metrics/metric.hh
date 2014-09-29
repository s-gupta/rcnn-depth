/*
 * Clustering metric.
 *
 * This is an abstract template class for a distance metric that can be used
 * in a clustering problem.
 */
#ifndef MLEARNING__CLUSTERING__METRICS__METRIC_HH
#define MLEARNING__CLUSTERING__METRICS__METRIC_HH

#include "collections/abstract/collection.hh"
#include "lang/pointers/auto_ptr.hh"

namespace mlearning {
namespace clustering {
namespace metrics {
/*
 * Imports.
 */
using collections::abstract::collection;
using lang::pointers::auto_ptr;

/*
 * Clustering metric.
 * T is the type of object being clustered.
 * V is the numeric type used to represent distance.
 */
template <typename T, typename V = double>
class metric {
public:
   /*
    * Destructor.
    */
   virtual ~metric() { }
   
   /*
    * Distance computation.
    * Return the distance from the item to the centroid.
    */
   virtual V distance(
      const T&,                  /* item */
      const T&                   /* centroid */
   ) const = 0;
   
   /*
    * Centroid computation.
    * Return the cluster centroid given the items in the cluster.
    * Note that the cluster must not be empty.
    */
   virtual auto_ptr<T> centroid(
      const collection<T>&,      /* items in cluster */
      const collection<double>&  /* weight of each item */
   ) const = 0;
};

} /* namespace metrics */
} /* namespace clustering */
} /* namespace mlearning */

#endif
