/*
 * Exemplar.
 *
 * An exemplar stores the visual information for a single example of an object.
 */
#ifndef VISION__RECOGNITION__MODELS__EXEMPLAR_HH
#define VISION__RECOGNITION__MODELS__EXEMPLAR_HH

#include "collections/array_list.hh"
#include "collections/pointers/auto_collection.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"
#include "lang/pointers/safe_ptr.hh"
#include "vision/features/geometric_blur.hh"
#include "vision/recognition/models/ids/exemplar_id.hh"

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
using vision::features::geometric_blur;
using vision::recognition::models::ids::exemplar_id;

/*
 * Exemplar.
 */
class exemplar {
public:
   /*
    * Constructor.
    * Create an exemplar from the given features.
    * Optionally set its identity.
    */
   explicit exemplar(
      auto_collection<geometric_blur>, /* geometric blur features */
      auto_ptr<exemplar_id> = auto_ptr<exemplar_id>(NULL)
   );

   /*
    * Copy constructor.
    * Perform a deep copy of the exemplar.
    * The copy does not have an id assigned unless specified.
    */
   exemplar(
      const exemplar&,
      auto_ptr<exemplar_id> = auto_ptr<exemplar_id>(NULL)
   );

   /*
    * Destructor.
    */
   virtual ~exemplar();

   /*
    * Get geometric blur features in exemplar.
    */
   const array_list<geometric_blur>& gb_features() const;
   
   /*
    * Get exemplar identity.
    */
   safe_ptr<const exemplar_id> exemplar_identity() const;

   /*
    * Set exemplar identity.
    * Return the new identity.
    */
   safe_ptr<const exemplar_id> exemplar_identity(auto_ptr<exemplar_id>);
   
protected:
   /*
    * Features.
    */
   auto_collection< geometric_blur, array_list<geometric_blur> > _gb_features;
   
   /*
    * Exemplar id.
    */
   auto_ptr<exemplar_id> _e_id;
};

} /* namespace models */
} /* namespace recognition */
} /* namespace vision */

#endif
