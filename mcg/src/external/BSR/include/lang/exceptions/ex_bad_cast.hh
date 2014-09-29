/*
 * Bad cast exception.
 */
#ifndef LANG__EXCEPTIONS__EX_BAD_CAST_HH
#define LANG__EXCEPTIONS__EX_BAD_CAST_HH

#include "lang/exceptions/exception.hh"
#include "lang/null.hh"

namespace lang {
namespace exceptions {

/*
 * Bad cast exception.
 */
class ex_bad_cast : public exception {
public:
   /*
    * Constructor.
    */
   explicit ex_bad_cast(
      const char* = NULL   /* message (use default if NULL) */
   );

   /*
    * Copy constructor.
    */
   ex_bad_cast(const ex_bad_cast&);

   /*
    * Destructor.
    */
   virtual ~ex_bad_cast();

   /*
    * Clone the exception.
    */
   virtual ex_bad_cast* clone() const;

   /*
    * Throw the exception.
    */
   virtual void raise() const;
};

} /* namespace exceptions */
} /* namespace lang */

#endif
