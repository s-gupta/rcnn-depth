/*
 * Feature.
 *
 * The feature class serves as a base class for image features that have a 
 * local orientation and scale.
 */
#ifndef VISION__FEATURES__FEATURE_HH
#define VISION__FEATURES__FEATURE_HH

#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"
#include "lang/pointers/safe_ptr.hh"
#include "vision/features/feature_id.hh"

namespace vision {
namespace features {
/*
 * Imports.
 */
using lang::pointers::auto_ptr;
using lang::pointers::safe_ptr;

/*
 * Base class for local image features.
 */
class feature {
public:
   /*
    * Constructor.
    * Create a feature at the given location in scale space.
    * Optionally set its identity.
    */
   explicit feature(
      double,                       /* x-coordinate */
      double,                       /* y-coordinate */
      double,                       /* orientation (in radians) */ 
      double,                       /* scale */
      auto_ptr<feature_id> = auto_ptr<feature_id>(NULL)
   );

   /*
    * Copy constructor.
    * The copy does not have an id assigned unless specified.
    */
   feature(
      const feature&, 
      auto_ptr<feature_id> = auto_ptr<feature_id>(NULL)
   );
   
   /*
    * Destructor.
    */
   virtual ~feature();

   /*
    * Get feature location.
    */
   double x() const;                /* x-coordinate */
   double y() const;                /* y-coordinate */
   double orientation() const;      /* orientation (in radians) */
   double scale() const;            /* scale */
   
   /*
    * Get feature identity.
    * Returns NULL if no identity has been assigned.
    */
   safe_ptr<const feature_id> feature_identity() const;
  
   /*
    * Set feature identity.
    * Return the new identity.
    */
   safe_ptr<const feature_id> feature_identity(auto_ptr<feature_id>);

protected:
   /*
    * Feature data.
    */
   double _x;                       /* x-coordinate */
   double _y;                       /* y-coordinate */
   double _ori;                     /* orientation (in radians) */
   double _scale;                   /* scale of feature */
   auto_ptr<feature_id> _f_id;      /* feature id */
};

} /* namespace features */
} /* namespace vision */

#endif
