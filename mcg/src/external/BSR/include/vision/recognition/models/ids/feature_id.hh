/*
 * Feature identity.
 *
 * This class provides identity information for a feature within the object
 * recognition framework.  It extends the abstract base feature_id class 
 * defined in the vision::features namespace.
 *
 * In the recognition framework, a feature_id contains a reference to its 
 * parent exemplar and indicates which feature it is within that exemplar.
 */
#ifndef VISION__RECOGNITION__MODELS__IDS__FEATURE_ID_HH
#define VISION__RECOGNITION__MODELS__IDS__FEATURE_ID_HH

#include "interfaces/comparable.hh"
#include "vision/features/feature_id.hh"

/*
 * Declare exemplar class.
 */
namespace vision {
namespace recognition {
namespace models {
   class exemplar;
} /* namespace models */
} /* namespace recognition */
} /* namespace vision */

namespace vision {
namespace recognition {
namespace models {
namespace ids {
/*
 * Imports.
 */
using interfaces::comparable;
using vision::recognition::models::exemplar;

/*
 * Feature identity.
 */
class feature_id : public vision::features::feature_id,
                   public comparable<feature_id> {
public:
   /*
    * Constructor.
    */
   feature_id(
      exemplar&,           /* exemplar to which feature belongs */
      unsigned long        /* feature number within exemplar */
   );

   /*
    * Copy constructor.
    */
   feature_id(const feature_id&);
   
   /*
    * Destructor.
    */
   virtual ~feature_id();

   /*
    * Get the exemplar to which the feature belongs.
    */
   exemplar& parent_exemplar() const;

   /*
    * Get the feature number within the exemplar.
    */
   unsigned long id() const;
   
   /*
    * Compare to another feature identity.
    */   
   virtual int compare_to(const feature_id&) const;

protected:
   exemplar& _exemplar;    /* exemplar to which feature belongs */
   unsigned long _id;      /* feature number within exemplar */
};

} /* namespace ids */
} /* namespace models */
} /* namespace recognition */
} /* namespace vision */

#endif
