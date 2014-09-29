/*
 * Throwable.
 *
 * Any object which can be thrown as an exception should derive from throwable.
 */
#ifndef LANG__EXCEPTIONS__THROWABLE_HH
#define LANG__EXCEPTIONS__THROWABLE_HH

namespace lang {
namespace exceptions {

/*
 * Throwable interface.
 */
class throwable {
public:
   /*
    * Destructor.
    */
   virtual ~throwable() = 0;

   /*
    * Clone the object.
    */
   virtual throwable* clone() const = 0;

   /*
    * Throw the object.
    */
   virtual void raise() const = 0;
};

} /* namespace exceptions */
} /* namespace lang */

#endif
