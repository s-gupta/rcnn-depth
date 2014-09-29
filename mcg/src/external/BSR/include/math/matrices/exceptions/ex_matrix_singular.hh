/*
 * Singular matrix exception.
 */
#ifndef MATH__MATRICES__EXCEPTIONS__MATRIX_SINGULAR_HH
#define MATH__MATRICES__EXCEPTIONS__MATRIX_SINGULAR_HH

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
 * Singular matrix exception.
 */
class ex_matrix_singular : public exception {
public:
   /*
    * Constructor.
    */
   explicit ex_matrix_singular(
      const char* = NULL   /* message (use default if NULL) */
   );

   /*
    * Copy constructor.
    */
   ex_matrix_singular(const ex_matrix_singular&);

   /*
    * Destructor.
    */
   virtual ~ex_matrix_singular();

   /*
    * Clone the exception.
    */
   virtual ex_matrix_singular* clone() const;

   /*
    * Throw the exception.
    */
   virtual void raise() const;
};

} /* namespace exceptions */
} /* namespace matrices */
} /* namespace math */

#endif
