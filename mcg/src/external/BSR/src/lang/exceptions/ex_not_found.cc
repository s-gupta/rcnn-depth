/*
 * Not found exception.
 */
#include "lang/exceptions/exception.hh"
#include "lang/exceptions/ex_not_found.hh"
#include "lang/null.hh"

namespace lang {
namespace exceptions {

/*
 * Constructor.
 */
ex_not_found::ex_not_found(const char* msg)
 : exception(
      ((msg != NULL) ? msg : "not found")
   )
{ }

/*
 * Copy constructor.
 */
ex_not_found::ex_not_found(const ex_not_found& e)
 : exception(e)
{ }

/*
 * Destructor.
 */
ex_not_found::~ex_not_found() {
   /* do nothing */
}

/*
 * Clone the exception.
 */
ex_not_found* ex_not_found::clone() const {
   return (new ex_not_found(*this));
}

/*
 * Throw the exception.
 */
void ex_not_found::raise() const {
   throw *this;
}

} /* namespace exceptions */
} /* namespace lang */
