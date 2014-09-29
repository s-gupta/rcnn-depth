/*
 * Region.
 */
#ifndef VISION__SEGMENTATION__REGION_HH
#define VISION__SEGMENTATION__REGION_HH

#include "collections/map.hh"

namespace vision {
namespace segmentation {
/*
 * Imports.
 */
using collections::map;

/*
 * Region.
 */
class region {
public:
   /*
    * Constructor.
    */
   region();

   /*
    * Copy constructor.
    */
   region(const region&);

   /*
    * Destructor.
    */
   virtual ~region();

   /*
    * Region id and connectivity.
    * Map of region id -> id of shared boundary.
    */
   unsigned long _id;
   map<unsigned long, unsigned long> _boundary_map;

   /*
    * Region data.
    */
   unsigned long _size;    /* number of pixels in region */
   double        _sum_L;   /* sums for L channel */
   double        _sum_L2;
   double        _sum_a;   /* sums for a channel */
   double        _sum_a2;
   double        _sum_b;   /* sums for b channel */
   double        _sum_b2;
};

} /* namespace segmentation */
} /* namespace vision */

#endif
