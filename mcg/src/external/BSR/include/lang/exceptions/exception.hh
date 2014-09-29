/*
 * Exception.
 */
#ifndef LANG__EXCEPTIONS__EXCEPTION_HH
#define LANG__EXCEPTIONS__EXCEPTION_HH

#include "io/streams/ostream.hh"
#include "lang/exceptions/throwable.hh"
#include "lang/null.hh"

namespace lang {
namespace exceptions {
/*
 * Imports.
 */
using io::streams::ostream;

/*
 * Exception base class.
 */
class exception : public throwable {
public:
   /*
    * Constructor.
    */
   explicit exception(
      const char* = NULL   /* message */
   );

   /*
    * Copy constructor.
    */
   exception(const exception&);

   /*
    * Destructor.
    */
   virtual ~exception();
   
   /*
    * Clone the exception.
    */
   virtual exception* clone() const;

   /*
    * Throw the exception.
    */
   virtual void raise() const;
    
   /*
    * Get the exception's message.
    */
   virtual const char* what() const;

   /*
    * Get the stack trace at the time the exception was thrown.
    */
   const char* const* stack_trace() const;

   /*
    * Output message and stack trace to stream.
    */
   friend ostream& operator<<(ostream&, const exception&);

protected:
   char*         _message;
   char**        _stack_trace;
   unsigned long _stack_trace_size;
};

} /* namespace exceptions */
} /* namespace lang */

#endif
