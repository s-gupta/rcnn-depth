/*
 * Unsynchronized class.
 */
#include "concurrent/threads/synchronization/synchronizables/synchronizable.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"

namespace concurrent {
namespace threads {
namespace synchronization {
namespace synchronizables {

/*
 * Constructor.
 */
unsynchronized::unsynchronized()
 : synchronizable(0)
{ }

/*
 * Private copy constructor.
 * (unsynchronized objects should not be copied)
 */
unsynchronized::unsynchronized(const unsynchronized&)
 : synchronizable(0)
{ }

/*
 * Destructor.
 */
unsynchronized::~unsynchronized() {
   /* do nothing */
}

/*
 * Claim exclusive access (abstract interface).
 */
void unsynchronized::abstract_lock() const {
   this->lock();
}

/*
 * Release exclusive access (abstract interface).
 */
void unsynchronized::abstract_unlock() const {
   this->unlock();
}

/*
 * Claim read access (abstract interface).
 */
void unsynchronized::abstract_read_lock() const {
   this->read_lock();
}

/*
 * Release read access (abstract interface).
 */
void unsynchronized::abstract_read_unlock() const {
   this->read_unlock();
}

/*
 * Claim write access (abstract interface).
 */
void unsynchronized::abstract_write_lock() const {
   this->write_lock();
}

/*
 * Release write access (abstract interface).
 */
void unsynchronized::abstract_write_unlock() const {
   this->write_unlock();
}

} /* namespace synchronizables */
} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */
