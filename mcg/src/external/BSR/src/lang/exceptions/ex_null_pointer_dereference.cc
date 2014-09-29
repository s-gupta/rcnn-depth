/*
 * Null pointer dereference exception.
 */
#include "lang/exceptions/exception.hh"
#include "lang/exceptions/ex_null_pointer_dereference.hh"
#include "lang/null.hh"

namespace lang {
namespace exceptions {

/*
 * Constructor.
 */
ex_null_pointer_dereference::ex_null_pointer_dereference(const char* msg)
 : exception(
      ((msg != NULL) ? msg : "null pointer dereference")
   )
{ }

/*
 * Copy constructor.
 */
ex_null_pointer_dereference::ex_null_pointer_dereference(
   const ex_null_pointer_dereference& e)
 : exception(e)
{ }

/*
 * Destructor.
 */
ex_null_pointer_dereference::~ex_null_pointer_dereference() {
   /* do nothing */
}

/*
 * Clone the exception.
 */
ex_null_pointer_dereference* ex_null_pointer_dereference::clone() const {
   return (new ex_null_pointer_dereference(*this));
}

/*
 * Throw the exception.
 */
void ex_null_pointer_dereference::raise() const {
   throw *this;
}

} /* namespace exceptions */
} /* namespace lang */
