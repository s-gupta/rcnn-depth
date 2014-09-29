/*
 * Classifier.
 *
 * A classifier takes an exemplar and determines the category (if any) to which 
 * it belongs.
 */
#ifndef VISION__RECOGNITION__CLASSIFIERS__ABSTRACT__CLASSIFIER_HH
#define VISION__RECOGNITION__CLASSIFIERS__ABSTRACT__CLASSIFIER_HH

#include "lang/pointers/safe_ptr.hh"
#include "vision/recognition/databases/category_db.hh"
#include "vision/recognition/models/exemplar.hh"
#include "vision/recognition/models/ids/category_id.hh"

namespace vision {
namespace recognition {
namespace classifiers {
namespace abstract {
/*
 * Imports.
 */
using lang::pointers::safe_ptr;
using vision::recognition::databases::category_db;
using vision::recognition::models::exemplar;
using vision::recognition::models::ids::category_id;

/*
 * Abstract base class for classifiers.
 */
class classifier {
public:
   /*
    * Destructor.
    */
   virtual ~classifier() = 0;

   /*
    * Classify the exemplar.
    * Return the identity of the category to which the exemplar belongs.
    * The returned id is NULL if the exemplar does not belong to any category.
    */
   virtual safe_ptr<const category_id> classify(
      const category_db&, 
      const exemplar&
   ) = 0;
};

} /* namespace abstract */
} /* namespace classifiers */
} /* namespace recognition */
} /* namespace vision */

#endif
