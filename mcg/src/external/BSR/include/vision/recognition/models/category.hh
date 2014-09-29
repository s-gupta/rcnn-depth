/*
 * Category.
 *
 * A category is a collection of exemplars sharing a common visual theme.
 * Categories may be further divided into subcategories.
 */
#ifndef VISION__RECOGNITION__MODELS__CATEGORY_HH
#define VISION__RECOGNITION__MODELS__CATEGORY_HH

#include "collections/array_list.hh"
#include "collections/pointers/auto_collection.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"
#include "lang/pointers/safe_ptr.hh"
#include "vision/recognition/models/exemplar.hh"
#include "vision/recognition/models/ids/category_id.hh"

namespace vision {
namespace recognition {
namespace models {
/*
 * Imports.
 */
using collections::array_list;
using collections::pointers::auto_collection;
using lang::pointers::auto_ptr;
using lang::pointers::safe_ptr;
using vision::recognition::models::ids::category_id;

/*
 * Category.
 */
class category {
public:
   /*
    * Constructor.
    * Create an empty category.
    * Optionally set its identity.
    */
   explicit category(auto_ptr<category_id> = auto_ptr<category_id>(NULL));

   /*
    * Copy constructor.
    * Perform a deep copy of the category.
    * The copy does not have an id assigned unless specified.
    */
   category(
      const category&, 
      auto_ptr<category_id> = auto_ptr<category_id>(NULL)
   );

   /*
    * Destructor.
    */
   virtual ~category();

   /*
    * Add an exemplar to the category.
    * Return a reference to the category.
    */
   category& add_exemplar(auto_ptr<exemplar>);

   /*
    * Add a subcategory to the category.
    * Return a reference to the category.
    */
   category& add_subcategory(auto_ptr<category>);
  
   /*
    * Get exemplars within the category.
    */
   const array_list<exemplar>& exemplars() const;

   /*
    * Get subcategories.
    */
   const array_list<category>& subcategories() const;
    
   /*
    * Get category identity.
    */
   safe_ptr<const category_id> category_identity() const;

   /*
    * Set category identity.
    * Return the new identity.
    */
   safe_ptr<const category_id> category_identity(auto_ptr<category_id>);

protected:
   /*
    * Exemplars within category.
    */
   auto_collection< exemplar, array_list<exemplar> > _exemplars;

   /*
    * Subcategories.
    */
   auto_collection< category, array_list<category> > _subcategories;

   /*
    * Category id.
    */
   auto_ptr<category_id> _c_id;
};

} /* namespace models */
} /* namespace recognition */
} /* namespace vision */

#endif
