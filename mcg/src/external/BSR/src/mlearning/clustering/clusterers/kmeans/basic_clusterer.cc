/*
 * Basic K-means clusterer.
 */
#include "lang/array.hh"
#include "mlearning/clustering/clusterers/kmeans/basic_clusterer.hh"

namespace mlearning {
namespace clustering {
namespace clusterers {
namespace kmeans {
/*
 * Imports.
 */
using lang::array;

/*
 * Pure virtual destructor.
 */
basic_clusterer_base::~basic_clusterer_base() { }

/*
 * Compute ids of changed centroids given change flags.
 */
array<unsigned long> basic_clusterer_base::compute_changed_ids(
   const array<bool>& has_changed)
{
   unsigned long n_clusters = has_changed.size();
   return basic_clusterer_base::compute_changed_ids(has_changed, n_clusters);
}

/*
 * Compute ids of changed centroids given change flags.
 * Specify the number of centroids to consider.
 */
array<unsigned long> basic_clusterer_base::compute_changed_ids(
   const array<bool>& has_changed, unsigned long n_clusters)
{
   /* compute how many centroids have changed */
   unsigned long n_changed = 0;
   for (unsigned long n = 0; n < n_clusters; n++) {
      if (has_changed[n])
         n_changed++;
   }
   /* get array of changed centroid ids */
   array<unsigned long> changed_ids(n_changed);
   for (unsigned long n = 0, chngd = 0; chngd < n_changed; n++) {
      if (has_changed[n]) {
         changed_ids[chngd] = n;
         chngd++;
      }
   }
   return changed_ids;
}

} /* namespace kmeans */
} /* namesapce clusterers */
} /* namespace clustering */
} /* namespace mlearning */
