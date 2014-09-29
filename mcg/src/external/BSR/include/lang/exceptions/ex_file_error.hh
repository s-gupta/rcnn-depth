/*
 * File error exception.
 */
#ifndef LANG__EXCEPTIONS__EX_FILE_ERROR_HH
#define LANG__EXCEPTIONS__EX_FILE_ERROR_HH

#include "lang/exceptions/exception.hh"
#include "lang/null.hh"

namespace lang {
namespace exceptions {

/*
 * File error exception.
 */
class ex_file_error : public exception {
public:
   /*
    * Constructor.
    */
   explicit ex_file_error(
      const char* = NULL   /* message (use default if NULL) */
   );

   /*
    * Copy constructor.
    */
   ex_file_error(const ex_file_error&);

   /*
    * Destructor.
    */
   virtual ~ex_file_error();

   /*
    * Clone the exception.
    */
   virtual ex_file_error* clone() const;

   /*
    * Throw the exception.
    */
   virtual void raise() const;
};

} /* namespace exceptions */
} /* namespace lang */

#endif
