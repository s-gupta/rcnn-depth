/*
 * Geometric blur parameters.
 */
#ifndef VISION__FEATURES__PARAMETERS__GEOMETRIC_BLUR_PARAMS_HH
#define VISION__FEATURES__PARAMETERS__GEOMETRIC_BLUR_PARAMS_HH

#include "lang/array.hh"

namespace vision {
namespace features {
namespace parameters {
/*
 * Imports.
 */
using lang::array;

/*
 * Geometric Blur parameters.
 */
class geometric_blur_params {
public:
   /*
    * Constructor.
    * Return the default set of parameters.
    */
   geometric_blur_params();

   /*
    * Constructor.
    * Specify parameters.
    */
   explicit geometric_blur_params(
      double,                       /* alpha */
      double,                       /* beta */
      const array<double>&,         /* radii */
      const array<unsigned long>&,  /* number of orientations at each radius */
      unsigned long                 /* number of signal channels */
   );

   /*
    * Copy constructor.
    */
   geometric_blur_params(const geometric_blur_params&);
   
   /*
    * Destructor.
    */
   virtual ~geometric_blur_params();

   /*
    * Get parameters used for controlling variation of blur with distance.
    * The blur amount at distance |x| is alpha*|x| + beta.
    */
   double alpha() const;
   double beta() const;

   /*
    * Get the set of radii at which geometric blur is sampled.
    */
   const array<double>& radii() const;

   /*
    * Get the number of orientations to sample at each radius.
    */
   const array<unsigned long>& orientations() const;

   /*
    * Get the number of signal channels.
    */
   unsigned long channels() const;

   /*
    * Return the number of bins in the geometric blur feature vector.
    */
   unsigned long descriptor_size() const;
   
   /*
    * Compute blur amount at the given distance from the feature center.
    */
   double blur_sigma(double) const;
  
   /*
    * Default parameter set.
    */
   static const geometric_blur_params default_parameters;
   
protected:
   /*
    * Parameters.
    */
   double               _alpha;     /* blur parameters */
   double               _beta;
   array<double>        _radii;     /* radii at which to sample */
   array<unsigned long> _oris;      /* number of orientation sampled at each radius */
   unsigned long        _channels;  /* number of signal channels */

   /*
    * Dimensionality of feature vector.
    * (determined by parameters above)
    */
   unsigned long _descriptor_size;  /* feature vector dimensionality */

   /*
    * Compute feature vector size from radii and orientation arrays.
    */
   void compute_descriptor_size();
};

} /* namespace parameters */
} /* namespace features */
} /* namespace vision */

#endif
