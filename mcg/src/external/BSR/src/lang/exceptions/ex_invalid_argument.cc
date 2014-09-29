/*
 * Invalid argument exception.
 */
#include "lang/exceptions/exception.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "lang/null.hh"

namespace lang {
namespace exceptions {

/*
 * Constructor.
 */
ex_invalid_argument::ex_invalid_argument(const char* msg)
 : exception(
      ((msg != NULL) ? msg : "invalid argument")
   )
{ }

/*
 * Copy constructor.
 */
ex_invalid_argument::ex_invalid_argument(const ex_invalid_argument& e)
 : exception(e)
{ }

/*
 * Destructor.
 */
ex_invalid_argument::~ex_invalid_argument() {
   /* do nothing */
}

/*
 * Clone the exception.
 */
ex_invalid_argument* ex_invalid_argument::clone() const {
   return (new ex_invalid_argument(*this));
}

/*
 * Throw the exception.
 */
void ex_invalid_argument::raise() const {
   throw *this;
}

} /* namespace exceptions */
} /* namespace lang */
