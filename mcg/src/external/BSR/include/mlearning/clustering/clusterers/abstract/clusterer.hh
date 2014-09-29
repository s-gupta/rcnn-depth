/*
 * Abstract clusterer.
 *
 * A clusterer computes an assignment of data items to clusters.
 */
#ifndef MLEARNING__CLUSTERING__CLUSTERERS__ABSTRACT__CLUSTERER_HH
#define MLEARNING__CLUSTERING__CLUSTERERS__ABSTRACT__CLUSTERER_HH

#include "collections/abstract/collection.hh"
#include "lang/array.hh"

namespace mlearning {
namespace clustering {
namespace clusterers {
namespace abstract {
/*
 * Imports.
 */
using collections::abstract::collection;
using lang::array;

/*
 * Base class containing code common to all clusterer<T> templates.
 */
class clusterer_base {
public:
   /*
    * Destructor.
    */
   virtual ~clusterer_base() = 0;
   
   /*
    * Compute the number of clusters given the cluster assignments.
    *
    * As cluster ids must be used in order, this is the same as one 
    * more than the maximum value of an entry in the assignments array.
    */
   static unsigned long number_of_clusters(
      const array<unsigned long>&   /* assignments */
   );

   /*
    * Compute the number of items assigned to each cluster 
    * given the assignments (map of item -> cluster id). 
    */
   static array<unsigned long> cluster_sizes(
      const array<unsigned long>&   /* assignments */
   );

protected:
   /*
    * Return the size of each problem in a series of coarse to fine
    * clustering problems.
    */
   static array<unsigned long> coarse_problem_sizes(
      unsigned long,                /* number of items */
      unsigned long,                /* minimum problem size */
      double = 0.5                  /* coarsening factor (in [0,1)) */
   );
};

/*
 * Abstract interface for clusterers.
 */
template <typename T>
class clusterer : public clusterer_base {
public:
   /*
    * Destructor.
    */
   virtual ~clusterer() = 0;
   
   /*
    * Clustering.
    *
    * Return the cluster assignments of the items in the 
    * order in which they appear in the given collection.
    *
    * Note that cluster ids must be used in order, 
    * starting from zero.
    */
   virtual array<unsigned long> cluster(
      const collection<T>&          /* items to cluster */
   ) const = 0;
};

/*
 * Pure virtual destructor.
 */
template <typename T>
clusterer<T>::~clusterer() { }

} /* namespace abstract */
} /* namesapce clusterers */
} /* namespace clustering */
} /* namespace mlearning */

#endif
