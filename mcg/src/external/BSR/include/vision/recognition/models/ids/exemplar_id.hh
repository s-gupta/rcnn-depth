/*
 * Exemplar identity.
 *
 * This class provides identity information for an exemplar within the object
 * recognition framework.  An exemplar_id contains a reference to the category 
 * to which it belongs and indicates which exemplar it is within that category.
 */
#ifndef VISION__RECOGNITION__MODELS__IDS__EXEMPLAR_ID_HH
#define VISION__RECOGNITION__MODELS__IDS__EXEMPLAR_ID_HH

#include "interfaces/comparable.hh"

/*
 * Declare category class.
 */
namespace vision {
namespace recognition {
namespace models {
   class category;
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
using vision::recognition::models::category;

/*
 * Exemplar identity.
 */
class exemplar_id : public comparable<exemplar_id> {
public:
   /*
    * Constructor.
    */
   exemplar_id(
      category&,           /* category to which exemplar belongs */
      unsigned long        /* exemplar number within category */
   );

   /*
    * Copy constructor.
    */
   exemplar_id(const exemplar_id&);
   
   /*
    * Destructor.
    */
   virtual ~exemplar_id();

   /*
    * Get the category to which the exemplar belongs.
    */
   category& parent_category() const;

   /*
    * Get the exemplar number within the category.
    */
   unsigned long id() const;
   
   /*
    * Compare to another exemplar identity.
    */   
   virtual int compare_to(const exemplar_id&) const;
  
protected:
   category& _category;    /* category to which exemplar belongs */
   unsigned long _id;      /* exemplar number within category */
};

} /* namespace ids */
} /* namespace models */
} /* namespace recognition */
} /* namespace vision */

#endif
