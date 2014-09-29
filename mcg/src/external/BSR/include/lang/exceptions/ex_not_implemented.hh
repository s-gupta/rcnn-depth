/*
 * Not implemented exception.
 */
#ifndef LANG__EXCEPTIONS__EX_NOT_IMPLEMENTED_HH
#define LANG__EXCEPTIONS__EX_NOT_IMPLEMENTED_HH

#include "lang/exceptions/exception.hh"
#include "lang/null.hh"

namespace lang {
namespace exceptions {

/*
 * Not implemented exception.
 */
class ex_not_implemented : public exception {
public:
   /*
    * Constructor.
    */
   explicit ex_not_implemented(
      const char* = NULL   /* message (use default if NULL) */
   );

   /*
    * Copy constructor.
    */
   ex_not_implemented(const ex_not_implemented&);

   /*
    * Destructor.
    */
   virtual ~ex_not_implemented();

   /*
    * Clone the exception.
    */
   virtual ex_not_implemented* clone() const;

   /*
    * Throw the exception.
    */
   virtual void raise() const;
};

} /* namespace exceptions */
} /* namespace lang */

#endif
