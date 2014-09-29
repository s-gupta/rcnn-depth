/*
 * Synchronized class.
 */
#include "concurrent/threads/synchronization/mutex.hh"
#include "concurrent/threads/synchronization/synchronizables/synchronizable.hh"
#include "concurrent/threads/synchronization/synchronizables/synchronized.hh"

namespace concurrent {
namespace threads {
namespace synchronization {
namespace synchronizables {
/*
 * Imports.
 */
using concurrent::threads::synchronization::mutex;

/*
 * Constructor.
 * Initialize mutex.
 */
synchronized::synchronized() 
 : synchronizable(),
   _mutex()
{ }

/*
 * Private copy constructor.
 * (synchronized objects should not be copied)
 */
synchronized::synchronized(const synchronized&)
 : synchronizable(),
   _mutex()
{ }

/*
 * Destructor.
 */
synchronized::~synchronized() {
   /* do nothing */
}

/*
 * Claim exclusive access (abstract interface).
 */
void synchronized::abstract_lock() const {
   this->lock();
}

/*
 * Release exclusive access (abstract interface).
 */
void synchronized::abstract_unlock() const {
   this->unlock();
}

/*
 * Claim read access (abstract interface).
 */
void synchronized::abstract_read_lock() const {
   this->read_lock();
}

/*
 * Release read access (abstract interface).
 */
void synchronized::abstract_read_unlock() const {
   this->read_unlock();
}

/*
 * Claim write access (abstract interface).
 */
void synchronized::abstract_write_lock() const {
   this->write_lock();
}

/*
 * Release write access (abstract interface).
 */
void synchronized::abstract_write_unlock() const {
   this->write_unlock();
}

} /* namespace synchronizables */
} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */
