/*
 * Feature identity.
 *
 * This abstract base class serves as an identity tag that can be attached 
 * to feature descriptors.  It can be extended to store useful identity 
 * information.
 */
#ifndef VISION__FEATURES__FEATURE_ID_HH
#define VISION__FEATURES__FEATURE_ID_HH

namespace vision {
namespace features {

/*
 * Abstract base class for feature identity.
 */
class feature_id {
public:
   /* 
    * Destructor.
    */
   virtual ~feature_id() = 0;
};

} /* namespace features */
} /* namespace vision */

#endif
