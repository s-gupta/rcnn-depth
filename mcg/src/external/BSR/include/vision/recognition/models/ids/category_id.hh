/*
 * Category identity.
 *
 * This class provides identity information for an category within the object
 * recognition framework.  An category_id contains a reference to the category 
 * to which it belongs and indicates which subcategory it is within that 
 * parent category.
 */
#ifndef VISION__RECOGNITION__MODELS__IDS__CATEGORY_ID_HH
#define VISION__RECOGNITION__MODELS__IDS__CATEGORY_ID_HH

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
 * Category identity.
 */
class category_id : public comparable<category_id> {
public:
   /*
    * Constructor.
    */
   category_id(
      category&,           /* parent category to which subcategory belongs */
      unsigned long        /* subcategory number within parent category */
   );

   /*
    * Copy constructor.
    */
   category_id(const category_id&);
   
   /*
    * Destructor.
    */
   virtual ~category_id();

   /*
    * Get the parent category to which the subcategory belongs.
    */
   category& parent_category() const;

   /*
    * Get the subcategory number within the parent category.
    */
   unsigned long id() const;
   
   /*
    * Compare to another category identity.
    */   
   virtual int compare_to(const category_id&) const;
   
protected:
   category& _category;    /* parent category to which sub category belongs */
   unsigned long _id;      /* subcategory number within parent category */
};

} /* namespace ids */
} /* namespace models */
} /* namespace recognition */
} /* namespace vision */

#endif
