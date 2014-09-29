/*
 * The synchronized class implements access control using a binary semaphore.
 * Only exclusive access is supported.  
 *
 * A request for read or write access is equivalent to a request for exclusive
 * access.
 */
#ifndef CONCURRENT__THREADS__SYNCHRONIZATION__SYNCHRONIZABLES__SYNCHRONIZED_HH
#define CONCURRENT__THREADS__SYNCHRONIZATION__SYNCHRONIZABLES__SYNCHRONIZED_HH

#include "concurrent/threads/synchronization/mutex.hh"
#include "concurrent/threads/synchronization/synchronizables/synchronizable.hh"

namespace concurrent {
namespace threads {
namespace synchronization {
namespace synchronizables {
/*
 * Imports.
 */
using concurrent::threads::synchronization::mutex;

class synchronized : public synchronizable {
public:
   /*
    * Constructor.
    */
   synchronized();
  
   /*
    * Destructor.
    */
   virtual ~synchronized();

   /*
    * Claim exclusive access.
    */
   virtual void abstract_lock() const;
   inline void lock() const;
   
   /*
    * Release exclusive access.
    */
   virtual void abstract_unlock() const;
   inline void unlock() const;

   /*
    * Claim read access.
    */
   virtual void abstract_read_lock() const;
   inline void read_lock() const;

   /*
    * Release read access.
    */
   virtual void abstract_read_unlock() const;
   inline void read_unlock() const;

   /*
    * Claim write access.
    */
   virtual void abstract_write_lock() const;
   inline void write_lock() const;

   /*
    * Release write access.
    */
   virtual void abstract_write_unlock() const;
   inline void write_unlock() const;

private:
   /*
    * Private copy constructor.
    * (synchronized objects should not be copied)
    */
   explicit synchronized(const synchronized&);
   
   /* binary semaphore mutex */
   mutable mutex _mutex;
};

/*
 * Claim exclusive access (inline version).
 */
inline void synchronized::lock() const {
   _mutex.lock();
}

/*
 * Release exclusive access (inline version).
 */
inline void synchronized::unlock() const {
   _mutex.unlock();
}

/*
 * Claim read access (inline version).
 */
inline void synchronized::read_lock() const {
   _mutex.lock();
}

/*
 * Release read access (inline version).
 */
inline void synchronized::read_unlock() const {
   _mutex.unlock();
}

/*
 * Claim write access (inline version).
 */
inline void synchronized::write_lock() const {
   _mutex.lock();
}

/*
 * Release write access (inline version).
 */
inline void synchronized::write_unlock() const {
   _mutex.unlock();
}

} /* namespace synchronizables */
} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */

#endif
