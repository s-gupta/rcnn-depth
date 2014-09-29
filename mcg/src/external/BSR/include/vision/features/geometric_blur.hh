/*
 * Geometric Blur.
 */
#ifndef VISION__FEATURES__GEOMETRIC_BLUR_HH
#define VISION__FEATURES__GEOMETRIC_BLUR_HH

#include "collections/abstract/collection.hh"
#include "interfaces/equalable.hh"
#include "interfaces/kd_treeable.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/matrices/matrix.hh"
#include "vision/features/feature_descriptor.hh"
#include "vision/features/feature_id.hh"
#include "vision/features/parameters/geometric_blur_params.hh"

namespace vision {
namespace features {
/*
 * Imports.
 */
using collections::abstract::collection;
using interfaces::equalable;
using interfaces::kd_tree_key;
using interfaces::kd_treeable;
using math::matrices::matrix;
using vision::features::parameters::geometric_blur_params;

/*
 * Geometric Blur.
 * This class associates a geometric blur feature with the parameter values 
 * used to compute it.  It is only meaningful to compare geometric blur 
 * features created using the same parameter settings.
 */
class geometric_blur : public feature_descriptor, 
                       public equalable<const geometric_blur>, 
                       public kd_treeable<geometric_blur> {
public:
   /*
    * Constructor.
    * Create a geometric blur feature as specified.
    * FIXME: this should be removed in the future 
    * (once c++ implementation of computing features is done).
    */
   explicit geometric_blur(
      double,                       /* x-coordinate */
      double,                       /* y-coordinate */
      double,                       /* orientation (in radians) */ 
      double,                       /* scale */
      const matrix<>&,              /* descriptor */
      auto_ptr<feature_id> = auto_ptr<feature_id>(NULL),
      const geometric_blur_params& = geometric_blur_params::default_parameters
   );
   
   /*
    * Copy constructor.
    * The copy does not have an id assigned unless specified.
    */
   geometric_blur(
      const geometric_blur&, 
      auto_ptr<feature_id> = auto_ptr<feature_id>(NULL)
   );

   /*
    * Destructor.
    */
   virtual ~geometric_blur();

   /*
    * Get the parameters associated with the feature.
    */
   const geometric_blur_params& parameters() const;
  
   /*
    * Distance from feature vector to half-space.
    */
   double distance_to(const kd_tree_key<>&) const;

   /*
    * Distance between feature vectors.
    * Throw an invalid argument exception if the features are not comparable
    * (do not share the same parameters).
    */
   double distance_to(const geometric_blur&) const;

   /*
    * Determine side of half-space in which feature vector lies.
    */
   int compare_to(const kd_tree_key<>&) const;

   /*
    * Check equality of feature vectors.
    * Throw an invalid argument exception if the features are not comparable
    * (do not share the same parameters).
    */
   bool is_equal_to(const geometric_blur&) const;

   /*
    * Return dimensionality of geometric blur feature vectors in collection.
    * Throw an invalid argument exception if the items in the collection do 
    * not share the same parameters (and hence have the same dimensionality).
    * An invalid argument exception is also thrown if the collection is empty.
    */
   static unsigned long dimensionality(const collection<geometric_blur>&);

   /*
    * Compute key for splitting a collection of geometric blur features 
    * along the dimension of the feature vector with highest variance.  
    * Throw an invalid argument exception if the features do not share 
    * the same parameters (or the collection is empty).
    */
   static kd_tree_key<> key_split(const collection<geometric_blur>&);
   
protected:
   /*
    * Parameters used when creating descriptor.
    */
   const geometric_blur_params& _params;

   /*
    * Assert that the given collection of geometric blur features is nonempty 
    * and that all features in the collection share the same parameters.
    * Throw an invalid argument exception if this assertation fails.
    */
   static void assert_nonempty_compatible(const collection<geometric_blur>&);
};

} /* namespace features */
} /* namespace vision */

#endif
