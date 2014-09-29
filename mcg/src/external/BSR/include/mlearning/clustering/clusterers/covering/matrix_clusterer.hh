/*
 * Covering matrix clusterer.
 *
 * Using the default clustering metric (L2 distance), each matrix is treated 
 * as an n-dimensional vector, where n is the total number of matrix elements.
 * All matrices to cluster must be of the same size and shape.
 */
#ifndef MLEARNING__CLUSTERING__CLUSTERERS__COVERING__MATRIX_CLUSTERER_HH
#define MLEARNING__CLUSTERING__CLUSTERERS__COVERING__MATRIX_CLUSTERER_HH

#include "math/matrices/matrix.hh"
#include "mlearning/clustering/clusterers/covering/basic_clusterer.hh"
#include "mlearning/clustering/metrics/matrix_metrics.hh"
#include "mlearning/clustering/metrics/metric.hh"

namespace mlearning {
namespace clustering {
namespace clusterers {
namespace covering {
/*
 * Imports.
 */
using math::matrices::matrix;
using mlearning::clustering::metrics::matrix_metrics;
using mlearning::clustering::metrics::metric;

/*
 * Covering matrix clusterer.
 */
template <typename T = double>
class matrix_clusterer : public basic_clusterer<matrix<T>,T> {
public:
   /*
    * Constructor.
    */
   explicit matrix_clusterer(
      T,                   /* cluster radius */
      unsigned long = 0,   /* maximum # of samples (0 for unlimited) */
      const metric<matrix<T>,T>& = matrix_metrics<T>::L2_metric()
   );

   /*
    * Copy constructor.
    * Create a clusterer with the same parameters.
    */
   matrix_clusterer(const matrix_clusterer<T>&);

   /*
    * Destructor.
    */
   virtual ~matrix_clusterer();
};

/*
 * Constructor.
 */
template <typename T>
matrix_clusterer<T>::matrix_clusterer(
   T             radius,
   unsigned long max_samples,
   const metric<matrix<T>,T>& m)
 : basic_clusterer<matrix<T>,T>(radius, max_samples, m)
{ }

/*
 * Copy constructor.
 * Create a clusterer with the same parameters.
 */
template <typename T>
matrix_clusterer<T>::matrix_clusterer(const matrix_clusterer<T>& c)
 : basic_clusterer<matrix<T>,T>(c)
{ }

/*
 * Destructor.
 */
template <typename T>
matrix_clusterer<T>::~matrix_clusterer() {
   /* do nothing */
}

} /* namespace covering */
} /* namespace clusterers */
} /* namespace clustering */
} /* namespace mlearning */

#endif
