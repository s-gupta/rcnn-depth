/*
 * Not found exception.
 */
#ifndef LANG__EXCEPTIONS__EX_NOT_FOUND_HH
#define LANG__EXCEPTIONS__EX_NOT_FOUND_HH

#include "lang/exceptions/exception.hh"
#include "lang/null.hh"

namespace lang {
namespace exceptions {

/*
 * Not found exception.
 */
class ex_not_found : public exception {
public:
   /*
    * Constructor.
    */
   explicit ex_not_found(
      const char* = NULL   /* message (use default if NULL) */
   );

   /*
    * Copy constructor.
    */
   ex_not_found(const ex_not_found&);

   /*
    * Destructor.
    */
   virtual ~ex_not_found();

   /*
    * Clone the exception.
    */
   virtual ex_not_found* clone() const;

   /*
    * Throw the exception.
    */
   virtual void raise() const;
};

} /* namespace exceptions */
} /* namespace lang */

#endif
