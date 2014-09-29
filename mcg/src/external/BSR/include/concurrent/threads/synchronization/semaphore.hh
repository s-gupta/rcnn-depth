/*
 * Semaphore.
 */
#ifndef CONCURRENT__THREADS__SYNCHRONIZATION__SEMAPHORE_HH
#define CONCURRENT__THREADS__SYNCHRONIZATION__SEMAPHORE_HH

#include "concurrent/threads/synchronization/mutex.hh"

namespace concurrent {
namespace threads {
namespace synchronization {

class semaphore {
public:
   /*
    * Constructor.
    * Initialize the semaphore with value 1.
    */
   semaphore();

   /*
    * Constructor.
    * Initialize the semaphore with the specified value.
    * Semaphores initialized to zero are initially locked.
    */
   explicit semaphore(unsigned long);

   /*
    * Copy constructor.
    * Create a new semaphore with the same value.
    */
   explicit semaphore(const semaphore&);

   /*
    * Destructor.
    */
   virtual ~semaphore();

   /*
    * Acquire the semaphore.
    */
   void acquire();
   void acquire(unsigned long);

   /*
    * Release the semaphore.
    */
   void release();
   void release(unsigned long);

protected:
   unsigned long _value;         /* number of remaining acquires allowed */
   mutable mutex _value_mutex;   /* mutex for reading/writing semaphore value */
   mutex         _access_mutex;  /* mutex for taking semaphore lock */
};

} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */

#endif
