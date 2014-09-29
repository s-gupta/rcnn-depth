/*
 * Boundary.
 *
 * A boundary separates two distinct regions.
 */
#ifndef VISION__SEGMENTATION__BOUNDARY_HH
#define VISION__SEGMENTATION__BOUNDARY_HH

namespace vision {
namespace segmentation {

/*
 * Boundary of a region.
 */
class boundary {
public:
   /*
    * Constructor.
    */
   boundary();

   /*
    * Destructor.
    */
   virtual ~boundary();

   /*
    * Boundary id.
    */
   unsigned long _id;
   
   /*
    * Boundary data.
    */
   unsigned long _size;    /* number of pixels along boundary */
   double _sum_contrast;   /* sum of local contrast measures */
   double _sum_pb;         /* sum of pb along boundary */
};

} /* namespace segmentation */
} /* namespace vision */

#endif
