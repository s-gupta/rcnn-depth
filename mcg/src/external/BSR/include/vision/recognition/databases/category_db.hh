/*
 * Category DB.
 *
 * A category database is a searchable database of all exemplars 
 * within a category and any of its descendent subcategories.
 */
#ifndef VISION__RECOGNITION__DATABASES__CATEGORY_DB_HH
#define VISION__RECOGNITION__DATABASES__CATEGORY_DB_HH

#include "collections/kd_tree.hh"
#include "lang/pointers/auto_ptr.hh"
#include "vision/features/geometric_blur.hh"
#include "vision/recognition/models/category.hh"

namespace vision {
namespace recognition {
namespace databases {
/*
 * Imports.
 */
using collections::kd_tree;
using lang::pointers::auto_ptr;
using vision::features::geometric_blur;
using vision::recognition::models::category;

/*
 * Category DB.
 */
class category_db {
public:
   /*
    * Constructor.
    * Build a database for the given category.
    */
   explicit category_db(auto_ptr<category>);

   /*
    * Copy constructor.
    * Make a deep copy of the database.
    */
   category_db(const category_db&);

   /*
    * Destructor.
    */
   virtual ~category_db();

   /*
    * Database search structure access.
    */
   const kd_tree<geometric_blur>& gb_features() const;
   
protected:
   /************************************************************************
    * Category DB data structures.
    ************************************************************************/
    
   /*
    * Root category.
    */
   auto_ptr<category> _root_category;

   /*
    * Search structures.
    */
   auto_ptr< kd_tree<geometric_blur> > _gb_db;  /* geometric blur feature database */
   
   /************************************************************************
    * Category DB helper functions.
    ************************************************************************/

   /*
    * Build the database search structures from the root category.
    */
   void build_db();
};
   
} /* namespace databases */
} /* namespace recognition */
} /* namespace vision */

#endif
