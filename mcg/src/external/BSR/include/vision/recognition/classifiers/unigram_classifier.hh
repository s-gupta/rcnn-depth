/*
 * Unigram classifier.
 *
 * Classify exemplars by treating them as a bag-of-features.
 * Each feature votes for categories and the exemplar is assigned to the 
 * most popular category selected by the features.
 */
#ifndef VISION__RECOGNITION__CLASSIFIERS__UNIGRAM_CLASSIFIER_HH
#define VISION__RECOGNITION__CLASSIFIERS__UNIGRAM_CLASSIFIER_HH

#include "lang/pointers/safe_ptr.hh"
#include "vision/recognition/databases/category_db.hh"
#include "vision/recognition/classifiers/abstract/classifier.hh"
#include "vision/recognition/models/exemplar.hh"
#include "vision/recognition/models/ids/category_id.hh"

namespace vision {
namespace recognition {
namespace classifiers {
/*
 * Imports.
 */
using lang::pointers::safe_ptr;
using vision::recognition::databases::category_db;
using vision::recognition::models::exemplar;
using vision::recognition::models::ids::category_id;

/*
 * Unigram classifier.
 */
class unigram_classifier : public abstract::classifier {
public:
   /*
    * Constructor.
    * Specify parameters for the classifier.
    */
   unigram_classifier(
      unsigned long,          /* number of approx nearest neighbors to vote */
      unsigned long           /* maximum neighbors to consider for each feature */
   );

   /*
    * Copy constructor.
    */
   unigram_classifier(const unigram_classifier&);

   /*
    * Destructor.
    */
   virtual ~unigram_classifier();

   /*
    * Classify exemplar.
    */
   virtual safe_ptr<const category_id> classify(
      const category_db&, 
      const exemplar&
   );

protected:
   unsigned long _num_nn;     /* number of approx nearest neighbors to vote */
   unsigned long _item_limit; /* maximum neighbors to consider for each feature */
};

} /* namespace classifiers */
} /* namespace recognition */
} /* namespace vision */

#endif
