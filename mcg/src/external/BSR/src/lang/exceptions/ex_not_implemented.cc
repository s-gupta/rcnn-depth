/*
 * Not implemented exception.
 */
#include "lang/exceptions/exception.hh"
#include "lang/exceptions/ex_not_implemented.hh"
#include "lang/null.hh"

namespace lang {
namespace exceptions {

/*
 * Constructor.
 */
ex_not_implemented::ex_not_implemented(const char* msg)
 : exception(
      ((msg != NULL) ? msg : "not implemented")
   )
{ }

/*
 * Copy constructor.
 */
ex_not_implemented::ex_not_implemented(const ex_not_implemented& e)
 : exception(e)
{ }

/*
 * Destructor.
 */
ex_not_implemented::~ex_not_implemented() {
   /* do nothing */
}

/*
 * Clone the exception.
 */
ex_not_implemented* ex_not_implemented::clone() const {
   return (new ex_not_implemented(*this));
}

/*
 * Throw the exception.
 */
void ex_not_implemented::raise() const {
   throw *this;
}

} /* namespace exceptions */
} /* namespace lang */
