/*
 * Matrix dimension mismatch exception.
 */
#include "lang/exceptions/exception.hh"
#include "lang/null.hh"
#include "math/matrices/exceptions/ex_matrix_dimension_mismatch.hh"

namespace math {
namespace matrices {
namespace exceptions {
/*
 * Imports.
 */
using lang::exceptions::exception;

/*
 * Constructor.
 */
ex_matrix_dimension_mismatch::ex_matrix_dimension_mismatch(const char* msg)
 : exception(
      ((msg != NULL) ? msg : "matrix dimension mismatch")
   )
{ }

/*
 * Copy constructor.
 */
ex_matrix_dimension_mismatch::ex_matrix_dimension_mismatch(
   const ex_matrix_dimension_mismatch& e)
 : exception(e)
{ }

/*
 * Destructor.
 */
ex_matrix_dimension_mismatch::~ex_matrix_dimension_mismatch() {
   /* do nothing */
}

/*
 * Clone the exception.
 */
ex_matrix_dimension_mismatch* ex_matrix_dimension_mismatch::clone() const {
   return (new ex_matrix_dimension_mismatch(*this));
}

/*
 * Throw the exception.
 */
void ex_matrix_dimension_mismatch::raise() const {
   throw *this;
}

} /* namespace exceptions */
} /* namespace matrices */
} /* namespace math */
