/*
 * Counter.
 * This class implements a thread-safe counter.
 */
#ifndef CONCURRENT__THREADS__SYNCHRONIZATION__COUNTER_HH
#define CONCURRENT__THREADS__SYNCHRONIZATION__COUNTER_HH

#include "concurrent/threads/synchronization/mutex.hh"

namespace concurrent {
namespace threads {
namespace synchronization {

/*
 * Counter class (thread-safe).
 */
class counter {
public:
   /*
    * Constructor.
    * Initialize the counter with value zero.
    */
   counter();

   /*
    * Constructor.
    * Initialize the counter with the specified value.
    */
   explicit counter(unsigned long);
   
   /*
    * Copy constructor.
    * Create a counter with the same value as the original.
    */
   explicit counter(const counter&);

   /*
    * Destructor.
    */
   virtual ~counter();

   /*
    * Get the counter value.
    */
   unsigned long get() const;

   /*
    * Set the counter value.
    */
   void set(unsigned long);
   
   /*
    * Increment the counter and return the new value.
    */
   unsigned long increment(unsigned long = 1 /* increment amount */);

   /*
    * Decrement the counter and return the new value.
    * A decrement by more than the counter value sets the counter to zero.
    */
   unsigned long decrement(unsigned long = 1 /* decrement amount */);

protected:
   unsigned long _value;         /* counter value */
   mutable mutex _value_mutex;   /* mutex for reading/writing counter value */
};

} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */

#endif
