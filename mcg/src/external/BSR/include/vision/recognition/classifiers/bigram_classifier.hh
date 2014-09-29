/*
 * Bigram classifier.
 *
 * Classify an exemplar by matching pairs of features in the exemplar to 
 * corresponding pairs of features in the training dataset.
 */
#ifndef VISION__RECOGNITION__CLASSIFIERS__BIGRAM_CLASSIFIER_HH
#define VISION__RECOGNITION__CLASSIFIERS__BIGRAM_CLASSIFIER_HH

#include "lang/pointers/safe_ptr.hh"
#include "vision/features/geometric_blur.hh"
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
using vision::features::geometric_blur;
using vision::recognition::databases::category_db;
using vision::recognition::models::exemplar;
using vision::recognition::models::ids::category_id;

/*
 * Bigram classifier.
 */
class bigram_classifier : public abstract::classifier {
public:
   /*
    * Constructor.
    * Specify parameters for the classifier.
    * FIXME: correct default parameters?
    */
   bigram_classifier(
      unsigned long,          /* number of approx nearest neighbors to use */
      unsigned long,          /* maximum neighbors to consider for each feature */
      double = 0.2,           /* max fractional change in bigram distance */
      double = 0.39270,       /* max absolute change in bigram angle */
      double = 50.0           /* minimum distance between training features in bigram */
   );

   /*
    * Copy constructor.
    */
   bigram_classifier(const bigram_classifier&);

   /*
    * Destructor.
    */
   virtual ~bigram_classifier();

   /*
    * Classify exemplar.
    */
   virtual safe_ptr<const category_id> classify(
      const category_db&, 
      const exemplar&
   );

protected:
   unsigned long _num_nn;     /* number of approx nearest neighbors to use */
   unsigned long _item_limit; /* maximum neighbors to consider for each feature */
   double _max_dist_change;
   double _max_angle_change;
   double _min_dist;

   /*
    * Bigram verification.
    */
   double verify_match(
      const geometric_blur&, 
      const geometric_blur&, 
      const geometric_blur&, 
      const geometric_blur&
   );
};

} /* namespace classifiers */
} /* namespace recognition */
} /* namespace vision */

#endif
