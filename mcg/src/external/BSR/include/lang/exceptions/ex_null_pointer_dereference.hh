/*
 * Null pointer dereference exception.
 */
#ifndef LANG__EXCEPTIONS__EX_NULL_POINTER_DEREFERENCE_HH
#define LANG__EXCEPTIONS__EX_NULL_POINTER_DEREFERENCE_HH

#include "lang/exceptions/exception.hh"
#include "lang/null.hh"

namespace lang {
namespace exceptions {

/*
 * Null pointer dereference exception.
 */
class ex_null_pointer_dereference : public exception {
public:
   /*
    * Constructor.
    */
   explicit ex_null_pointer_dereference(
      const char* = NULL   /* message (use default if NULL) */
   );

   /*
    * Copy constructor.
    */
   ex_null_pointer_dereference(const ex_null_pointer_dereference&);

   /*
    * Destructor.
    */
   virtual ~ex_null_pointer_dereference();

   /*
    * Clone the exception.
    */
   virtual ex_null_pointer_dereference* clone() const;

   /*
    * Throw the exception.
    */
   virtual void raise() const;
};

} /* namespace exceptions */
} /* namespace lang */

#endif
