/*
 * File read error exception.
 */
#include "lang/exceptions/ex_file_error.hh"
#include "lang/exceptions/ex_file_read_error.hh"
#include "lang/null.hh"
#include "lang/string.hh"

namespace lang {
namespace exceptions {
/*
 * Imports.
 */
using lang::string;

/*
 * Constructor.
 */
ex_file_read_error::ex_file_read_error(const char* msg)
 : ex_file_error(
      ((msg != NULL) ? msg : "file read error")
   )
{ }

/*
 * Constructor.
 * Generate an error message reporting the filename and error description.
 */
ex_file_read_error::ex_file_read_error(const char* filename, const char* err)
 : ex_file_error(
      string<>("error reading ").concat(filename).concat(": ").concat(err)
   )
{ }

/*
 * Copy constructor.
 */
ex_file_read_error::ex_file_read_error(const ex_file_read_error& e)
 : ex_file_error(e)
{ }

/*
 * Destructor.
 */
ex_file_read_error::~ex_file_read_error() {
   /* do nothing */
}

/*
 * Clone the exception.
 */
ex_file_read_error* ex_file_read_error::clone() const {
   return (new ex_file_read_error(*this));
}

/*
 * Throw the exception.
 */
void ex_file_read_error::raise() const {
   throw *this;
}

} /* namespace exceptions */
} /* namespace lang */
