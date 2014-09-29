/*
 * File read error exception.
 */
#ifndef LANG__EXCEPTIONS__EX_FILE_READ_ERROR_HH
#define LANG__EXCEPTIONS__EX_FILE_READ_ERROR_HH

#include "lang/exceptions/ex_file_error.hh"
#include "lang/null.hh"

namespace lang {
namespace exceptions {

/*
 * File read error exception.
 */
class ex_file_read_error : public ex_file_error {
public:
   /*
    * Constructor.
    */
   explicit ex_file_read_error(
      const char* = NULL   /* message (use default if NULL) */
   );

   /*
    * Constructor.
    * Generate an error message reporting the filename and error description.
    */
   explicit ex_file_read_error(
      const char*,         /* filename */
      const char*          /* error description */
   );
   
   /*
    * Copy constructor.
    */
   ex_file_read_error(const ex_file_read_error&);

   /*
    * Destructor.
    */
   virtual ~ex_file_read_error();

   /*
    * Clone the exception.
    */
   virtual ex_file_read_error* clone() const;

   /*
    * Throw the exception.
    */
   virtual void raise() const;
};

} /* namespace exceptions */
} /* namespace lang */

#endif
