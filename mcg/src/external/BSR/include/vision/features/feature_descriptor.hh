/*
 * Feature descriptor.
 *
 * A feature descriptor is a feature that has an associated appearance 
 * descriptor.
 */
#ifndef VISION__FEATURES__FEATURE_DESCRIPTOR_HH
#define VISION__FEATURES__FEATURE_DESCRIPTOR_HH

#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/matrices/matrix.hh"
#include "vision/features/feature.hh"
#include "vision/features/feature_id.hh"

namespace vision {
namespace features {
/*
 * Imports.
 */
using lang::pointers::auto_ptr;
using math::matrices::matrix;

/*
 * Feature descriptor.
 */
class feature_descriptor : public feature {
public:
   /*
    * Constructor.
    * Create a feature descriptor at the given location in scale space.
    * Optionally set its identity.
    */
   explicit feature_descriptor(
      double,                       /* x-coordinate */
      double,                       /* y-coordinate */
      double,                       /* orientation (in radians) */ 
      double,                       /* scale */
      const matrix<>&,              /* descriptor */
      auto_ptr<feature_id> = auto_ptr<feature_id>(NULL)
   );
   
   /*
    * Copy constructor.
    * The copy does not have an id assigned unless specified.
    */
   feature_descriptor(
      const feature_descriptor&, 
      auto_ptr<feature_id> = auto_ptr<feature_id>(NULL)
   );

   /*
    * Destructor.
    */
   virtual ~feature_descriptor();
   
   /*
    * Get the descriptor.
    */
   const matrix<>& descriptor() const;
   
protected:
   matrix<> _descriptor;
};

} /* namespace features */
} /* namespace vision */

#endif
