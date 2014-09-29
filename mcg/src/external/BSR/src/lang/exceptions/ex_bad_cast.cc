/*
 * Bad cast exception.
 */
#include "lang/exceptions/exception.hh"
#include "lang/exceptions/ex_bad_cast.hh"
#include "lang/null.hh"

namespace lang {
namespace exceptions {

/*
 * Constructor.
 */
ex_bad_cast::ex_bad_cast(const char* msg)
 : exception(
      ((msg != NULL) ? msg : "bad cast")
   )
{ }

/*
 * Copy constructor.
 */
ex_bad_cast::ex_bad_cast(const ex_bad_cast& e)
 : exception(e)
{ }

/*
 * Destructor.
 */
ex_bad_cast::~ex_bad_cast() {
   /* do nothing */
}

/*
 * Clone the exception.
 */
ex_bad_cast* ex_bad_cast::clone() const {
   return (new ex_bad_cast(*this));
}

/*
 * Throw the exception.
 */
void ex_bad_cast::raise() const {
   throw *this;
}

} /* namespace exceptions */
} /* namespace lang */
