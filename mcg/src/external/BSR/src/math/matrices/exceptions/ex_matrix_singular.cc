/*
 * Singular matrix exception.
 */
#include "lang/exceptions/exception.hh"
#include "lang/null.hh"
#include "math/matrices/exceptions/ex_matrix_singular.hh"

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
ex_matrix_singular::ex_matrix_singular(const char* msg)
 : exception(
      ((msg != NULL) ? msg : "singular matrix")
   )
{ }

/*
 * Copy constructor.
 */
ex_matrix_singular::ex_matrix_singular(const ex_matrix_singular& e)
 : exception(e)
{ }

/*
 * Destructor.
 */
ex_matrix_singular::~ex_matrix_singular() {
   /* do nothing */
}

/*
 * Clone the exception.
 */
ex_matrix_singular* ex_matrix_singular::clone() const {
   return (new ex_matrix_singular(*this));
}

/*
 * Throw the exception.
 */
void ex_matrix_singular::raise() const {
   throw *this;
}

} /* namespace exceptions */
} /* namespace matrices */
} /* namespace math */
