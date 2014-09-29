/*
 * Semaphore.
 */
#include "concurrent/threads/synchronization/mutex.hh"
#include "concurrent/threads/synchronization/semaphore.hh"

namespace concurrent {
namespace threads {
namespace synchronization {

/*
 * Default constructor.
 * Initialize the semaphore with value 1.
 */
semaphore::semaphore() 
 : _value(1), 
   _value_mutex(), 
   _access_mutex() { }

/*
 * Constructor.
 * Initialize the semaphore with the specified value.
 * Semaphores initialized to zero are initially locked.
 */
semaphore::semaphore(unsigned long n) 
 : _value(n), 
   _value_mutex(), 
   _access_mutex() 
{
   /* lock semaphore if needed */
   if (n == 0)
      _access_mutex.lock();
}

/*
 * Copy constructor.
 * Create a new semaphore with the same value.
 */
semaphore::semaphore(const semaphore& s)
 : _value(0), 
   _value_mutex(), 
   _access_mutex() 
{
   /* copy semaphore value */
   s._value_mutex.lock();
   _value = s._value;
   s._value_mutex.unlock();
   /* lock semaphore if needed */
   if (_value == 0)
      _access_mutex.lock();
}

/*
 * Destructor.
 */
semaphore::~semaphore() {
   /* do nothing  */
}

/*
 * Acquire the semaphore.
 * Wait until the semaphore value is greater than 0, then decrement it by 1.
 */
void semaphore::acquire() {
   _access_mutex.lock();
   _value_mutex.lock();
   _value--;
   bool unlock = (_value > 0);
   _value_mutex.unlock();
   if (unlock)
      _access_mutex.unlock();
}

/*
 * Acquire the semaphore n times (decrementing its value by 1 each time).
 */
void semaphore::acquire(unsigned long n) {
   for (unsigned long i = 0; i < n; i++)
      this->acquire();
}

/*
 * Release the semaphore.
 * Increase its value by one.
 */
void semaphore::release() {
   _value_mutex.lock();
   bool unlock = (_value == 0);
   _value++;
   _value_mutex.unlock();
   if (unlock)
      _access_mutex.unlock();
}

/*
 * Release the semaphore.
 * Increase its value by n.
 */
void semaphore::release(unsigned long n) {
   _value_mutex.lock();
   bool unlock = ((_value == 0) && (n > 0));
   _value += n;
   _value_mutex.unlock();
   if (unlock)
      _access_mutex.unlock();
}

} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */
