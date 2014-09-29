/*
 * Counter.
 * This class implements a thread-safe counter.
 */
#include "concurrent/threads/synchronization/counter.hh"

namespace concurrent {
namespace threads {
namespace synchronization {

/*
 * Constructor.
 * Initialize the counter with value zero.
 */
counter::counter()
 : _value(0), 
   _value_mutex()
{ }

/*
 * Constructor.
 * Initialize the counter with the specified value.
 */
counter::counter(unsigned long n)
 : _value(n),
   _value_mutex()
{ }

/*
 * Copy constructor.
 * Create a new counter with the same value.
 */
counter::counter(const counter& c) 
 : _value(c.get()), 
   _value_mutex()
{ }

/*
 * Destructor.
 */
counter::~counter() {
   /* do nothing */
}

/*
 * Get the counter value.
 */
unsigned long counter::get() const {
   _value_mutex.lock();
   unsigned long n = _value;
   _value_mutex.unlock();
   return n;
}

/*
 * Set the counter value.
 */
void counter::set(unsigned long n) {
   _value_mutex.lock();
   _value = n;
   _value_mutex.unlock();
}

/*
 * Increment the counter and return the new value.
 */
unsigned long counter::increment(unsigned long inc) {
   _value_mutex.lock();
   _value += inc;
   unsigned long n = _value;
   _value_mutex.unlock();
   return n;
}

/*
 * Decrement the counter and return the new value.
 * A decrement by more than the counter value sets the counter to zero.
 */
unsigned long counter::decrement(unsigned long dec) {
   _value_mutex.lock();
   if (_value > dec)
      _value -= dec;
   else
      _value = 0;
   unsigned long n = _value;
   _value_mutex.unlock();
   return n;
}

} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */
