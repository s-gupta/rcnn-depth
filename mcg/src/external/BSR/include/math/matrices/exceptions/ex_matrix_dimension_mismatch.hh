/*
 * Matrix dimension mismatch exception.
 */
#ifndef MATH__MATRICES__EXCEPTIONS__MATRIX_DIMENSION_MISMATCH_HH
#define MATH__MATRICES__EXCEPTIONS__MATRIX_DIMENSION_MISMATCH_HH

#include "lang/exceptions/exception.hh"
#include "lang/null.hh"

namespace math {
namespace matrices {
namespace exceptions {
/*
 * Imports.
 */
using lang::exceptions::exception;

/*
 * Matrix dimension mismatch exception.
 */
class ex_matrix_dimension_mismatch : public exception {
public:
   /*
    * Constructor.
    */
   explicit ex_matrix_dimension_mismatch(
      const char* = NULL   /* message (use default if NULL) */
   );

   /*
    * Copy constructor.
    */
   ex_matrix_dimension_mismatch(const ex_matrix_dimension_mismatch&);

   /*
    * Destructor.
    */
   virtual ~ex_matrix_dimension_mismatch();

   /*
    * Clone the exception.
    */
   virtual ex_matrix_dimension_mismatch* clone() const;

   /*
    * Throw the exception.
    */
   virtual void raise() const;
};

} /* namespace exceptions */
} /* namespace matices */
} /* namespace math */

#endif
