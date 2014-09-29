/*
 * Segmentation.
 */
#ifndef VISION__SEGMENTATION__SEGMENTATION_HH
#define VISION__SEGMENTATION__SEGMENTATION_HH

#include "collections/array_list.hh"
#include "collections/pointers/auto_collection.hh"
#include "lang/array.hh"
#include "math/matrices/matrix.hh"
#include "vision/segmentation/boundary.hh"
#include "vision/segmentation/region.hh"

namespace vision {
namespace segmentation {
/*
 * Imports.
 */
using collections::array_list;
using collections::pointers::auto_collection;
using lang::array;
using math::matrices::matrix;

/*
 * Compute region features given an oversegmentation.
 * Return the regions and boundaries.
 */
void compute_regions(
   const matrix<unsigned long>&, /* oversegmentation pixel assignments */
   const matrix<>&,              /* image L channel */
   const matrix<>&,              /* image a channel */
   const matrix<>&,              /* image b channel */
   const matrix<>&,              /* local contrast at each pixel */
   const matrix<>&,              /* pb at each pixel */
   auto_collection< region,   array_list<region> >&,  /* returned regions */
   auto_collection< boundary, array_list<boundary> >& /* returned boundaries */
);

/*
 * Create a segmentation given initial regions and boundaries.
 * Return the assignment of initial regions -> final regions.
 */
array<unsigned long> segment(
   auto_collection< region,   array_list<region> >,  /* regions */
   auto_collection< boundary, array_list<boundary> > /* boundaries */
);

} /* namespace segmentation */
} /* namespace vision */

#endif
