/*
 * Abstract weighted clusterer.
 */
#include "collections/abstract/collection.hh"
#include "lang/array.hh"
#include "lang/iterators/iterator.hh"
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
using lang::array;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

/*
 * Pure virtual destructor.
 */
weighted_clusterer_base::~weighted_clusterer_base() { }

/*
 * Compute the total weight assigned to each cluster (the sum 
 * of the weights of its members) given the item weights and 
 * assignment map.
 */
array<double> weighted_clusterer_base::cluster_weights(
   const collection<double>&   weights,
   const array<unsigned long>& assignments)
{
   unsigned long n_items = assignments.size();
   unsigned long n_clusters = clusterer_base::number_of_clusters(assignments);
   array<double> total_weights(n_clusters);
   auto_ptr< iterator<double> > i = weights.iter_create();
   for (unsigned long n = 0; n < n_items; n++)
      total_weights[assignments[n]] += i->next();
   return total_weights;
}

} /* namespace abstract */
} /* namesapce clusterers */
} /* namespace clustering */
} /* namespace mlearning */
