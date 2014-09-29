/*
 * Invalid argument exception.
 */
#ifndef LANG__EXCEPTIONS__EX_INVALID_ARGUMENT_HH
#define LANG__EXCEPTIONS__EX_INVALID_ARGUMENT_HH

#include "lang/exceptions/exception.hh"
#include "lang/null.hh"

namespace lang {
namespace exceptions {

/*
 * Invalid argument exception.
 */
class ex_invalid_argument : public exception {
public:
   /*
    * Constructor.
    */
   explicit ex_invalid_argument(
      const char* = NULL   /* message (use default if NULL) */
   );

   /*
    * Copy constructor.
    */
   ex_invalid_argument(const ex_invalid_argument&);

   /*
    * Destructor.
    */
   virtual ~ex_invalid_argument();

   /*
    * Clone the exception.
    */
   virtual ex_invalid_argument* clone() const;

   /*
    * Throw the exception.
    */
   virtual void raise() const;
};

} /* namespace exceptions */
} /* namespace lang */

#endif
