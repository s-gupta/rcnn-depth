/*
 * Abstract clusterer.
 */
#include "lang/array.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "math/math.hh"
#include "mlearning/clustering/clusterers/abstract/clusterer.hh"

namespace mlearning {
namespace clustering {
namespace clusterers {
namespace abstract {
/*
 * Imports.
 */
using lang::array;
using lang::exceptions::ex_invalid_argument;

/*
 * Pure virtual destructor.
 */
clusterer_base::~clusterer_base() { }

/*
 * Compute the number of clusters given the cluster assignments.
 *
 * As cluster ids must be used in order, this is the same as one 
 * more than the maximum value of an entry in the assignments array.
 */
unsigned long clusterer_base::number_of_clusters(
   const array<unsigned long>& assignments)
{
   unsigned long n_items = assignments.size();
   if (n_items > 0) {
      unsigned long max_cluster_num = 0;
      for (unsigned long n = 0; n < n_items; n++) {
         if (assignments[n] > max_cluster_num)
            max_cluster_num = assignments[n];
      }
      return max_cluster_num + 1;
   } else {
      return 0;
   }
}

/*
 * Compute the number of items assigned to each cluster 
 * given the assignments (map of item -> cluster id). 
 */
array<unsigned long> clusterer_base::cluster_sizes(
   const array<unsigned long>& assignments)
{
   unsigned long n_items = assignments.size();
   unsigned long n_clusters = clusterer_base::number_of_clusters(assignments);
   array<unsigned long> sizes(n_clusters);
   for (unsigned long n = 0; n < n_items; n++)
      sizes[assignments[n]]++;
   return sizes;
}

/*
 * Return the size of each problem in a series of coarse to fine
 * clustering problems.
 */
array<unsigned long> clusterer_base::coarse_problem_sizes(
   unsigned long n_items, unsigned long min_size, double factor)
{
   /* check arguments */
   if ((factor < 0) || (factor >= 1))
      throw ex_invalid_argument("coarsening factor must be in [0,1)");
   /* compute number of problems */
   unsigned long n_problems = 0;
   unsigned long curr_size = n_items;
   do {
      n_problems++;
      curr_size = static_cast<unsigned long>(
         math::floor(static_cast<double>(curr_size) * factor)
      );
   } while (curr_size >= min_size);
   /* store size of each problem */
   array<unsigned long> coarse_sizes(n_problems);
   curr_size = n_items;
   do {
      n_problems--;
      coarse_sizes[n_problems] = curr_size;
      curr_size = static_cast<unsigned long>(
         math::floor(static_cast<double>(curr_size) * factor)
      );
   } while (curr_size >= min_size);
   return coarse_sizes;
}

} /* namespace abstract */
} /* namesapce clusterers */
} /* namespace clustering */
} /* namespace mlearning */
