/*
 * File open error exception.
 */
#include "lang/exceptions/ex_file_error.hh"
#include "lang/exceptions/ex_file_open_error.hh"
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
ex_file_open_error::ex_file_open_error(const char* msg)
 : ex_file_error(
      ((msg != NULL) ? msg : "file open error")
   )
{ }

/*
 * Constructor.
 * Generate an error message reporting the filename and mode.
 */
ex_file_open_error::ex_file_open_error(const char* filename, const char* mode)
 : ex_file_error(
      string<>("unable to open ").concat(filename).concat(" for ").concat(mode)
   )
{ }

/*
 * Copy constructor.
 */
ex_file_open_error::ex_file_open_error(const ex_file_open_error& e)
 : ex_file_error(e)
{ }

/*
 * Destructor.
 */
ex_file_open_error::~ex_file_open_error() {
   /* do nothing */
}

/*
 * Clone the exception.
 */
ex_file_open_error* ex_file_open_error::clone() const {
   return (new ex_file_open_error(*this));
}

/*
 * Throw the exception.
 */
void ex_file_open_error::raise() const {
   throw *this;
}

} /* namespace exceptions */
} /* namespace lang */
