/*
 * File error exception.
 */
#include "lang/exceptions/exception.hh"
#include "lang/exceptions/ex_file_error.hh"
#include "lang/null.hh"

namespace lang {
namespace exceptions {

/*
 * Constructor.
 */
ex_file_error::ex_file_error(const char* msg)
 : exception(
      ((msg != NULL) ? msg : "file error")
   )
{ }

/*
 * Copy constructor.
 */
ex_file_error::ex_file_error(const ex_file_error& e)
 : exception(e)
{ }

/*
 * Destructor.
 */
ex_file_error::~ex_file_error() {
   /* do nothing */
}

/*
 * Clone the exception.
 */
ex_file_error* ex_file_error::clone() const {
   return (new ex_file_error(*this));
}

/*
 * Throw the exception.
 */
void ex_file_error::raise() const {
   throw *this;
}

} /* namespace exceptions */
} /* namespace lang */
