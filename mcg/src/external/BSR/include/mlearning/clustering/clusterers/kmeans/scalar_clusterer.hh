/*
 * K-means scalar clusterer.
 */
#ifndef MLEARNING__CLUSTERING__CLUSTERERS__KMEANS__SCALAR_CLUSTERER_HH
#define MLEARNING__CLUSTERING__CLUSTERERS__KMEANS__SCALAR_CLUSTERER_HH

#include "mlearning/clustering/clusterers/kmeans/basic_clusterer.hh"
#include "mlearning/clustering/metrics/scalar_metrics.hh"
#include "mlearning/clustering/metrics/metric.hh"

namespace mlearning {
namespace clustering {
namespace clusterers {
namespace kmeans {
/*
 * Imports.
 */
using mlearning::clustering::metrics::scalar_metrics;
using mlearning::clustering::metrics::metric;

/*
 * K-means scalar clusterer.
 */
template <typename T = double>
class scalar_clusterer : public basic_clusterer<T,T> {
public:
   /*
    * Constructor.
    */
   explicit scalar_clusterer(
      unsigned long,       /* K (desired number of clusters) */
      unsigned long = 0,   /* maximum # of iterations (0 for unlimited) */
      const metric<T,T>& = scalar_metrics<T>::L2_metric()
   );

   /*
    * Copy constructor.
    * Create a clusterer with the same parameters.
    */
   scalar_clusterer(const scalar_clusterer<T>&);

   /*
    * Destructor.
    */
   virtual ~scalar_clusterer();
};

/*
 * Constructor.
 */
template <typename T>
scalar_clusterer<T>::scalar_clusterer(
   unsigned long      K,
   unsigned long      max_iterations,
   const metric<T,T>& m)
 : basic_clusterer<T,T>(K, max_iterations, m)
{ }

/*
 * Copy constructor.
 * Create a clusterer with the same parameters.
 */
template <typename T>
scalar_clusterer<T>::scalar_clusterer(const scalar_clusterer<T>& c)
 : basic_clusterer<T,T>(c)
{ }

/*
 * Destructor.
 */
template <typename T>
scalar_clusterer<T>::~scalar_clusterer() {
   /* do nothing */
}

} /* namespace kmeans */
} /* namespace clusterers */
} /* namespace clustering */
} /* namespace mlearning */

#endif
