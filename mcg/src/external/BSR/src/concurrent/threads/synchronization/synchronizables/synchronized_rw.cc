/*
 * Synchronized_rw class.
 */
#include "concurrent/threads/synchronization/mutex.hh"
#include "concurrent/threads/synchronization/synchronizables/synchronizable.hh"
#include "concurrent/threads/synchronization/synchronizables/synchronized_rw.hh"
#include "lang/null.hh"

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
 * Initialize locks.
 */
synchronized_rw::synchronized_rw()
 : synchronizable(),
   _n_readers(0), 
   _rcnt_mutex(),
   _read_mutex(),
   _write_mutex() 
{ }

/*
 * Private copy constructor.
 * (synchronized_rw objects should not be copied)
 */
synchronized_rw::synchronized_rw(const synchronized_rw&)
 : synchronizable(),
   _n_readers(0), 
   _rcnt_mutex(),
   _read_mutex(),
   _write_mutex() 
{ }

/*
 * Destructor.
 */
synchronized_rw::~synchronized_rw() {
   /* do nothing */
}

/*
 * Claim exclusive access (abstract interface).
 */
void synchronized_rw::abstract_lock() const {
   this->lock();
}

/*
 * Release exclusive access (abstract interface).
 */
void synchronized_rw::abstract_unlock() const {
   this->unlock();
}

/*
 * Claim read access (abstract interface).
 */
void synchronized_rw::abstract_read_lock() const {
   this->read_lock();
}

/*
 * Release read access (abstract interface).
 */
void synchronized_rw::abstract_read_unlock() const {
   this->read_unlock();
}

/*
 * Claim write access (abstract interface).
 */
void synchronized_rw::abstract_write_lock() const {
   this->write_lock();
}

/*
 * Release write access (abstract interface).
 */
void synchronized_rw::abstract_write_unlock() const {
   this->write_unlock();
}

/*
 * Claim exclusive access. 
 */
void synchronized_rw::lock() const {
   _read_mutex.lock();
   _write_mutex.lock();
}

/*
 * Release exclusive access.
 */
void synchronized_rw::unlock() const {
   _write_mutex.unlock();
   _read_mutex.unlock();
}

/*
 * Claim read access.
 */
void synchronized_rw::read_lock()  const {
   _read_mutex.lock();
   _rcnt_mutex.lock();
   if (_n_readers == 0)
      _write_mutex.lock();
   _n_readers++;
   _rcnt_mutex.unlock();
   _read_mutex.unlock();
}

/*
 * Release read access.
 */
void synchronized_rw::read_unlock() const {
   _rcnt_mutex.lock();
   _n_readers--;
   if (_n_readers == 0)
      _write_mutex.unlock();
   _rcnt_mutex.unlock();
}

/*
 * Claim write access.
 */
void synchronized_rw::write_lock() const {
   _read_mutex.lock();
   _write_mutex.lock();
}

/*
 * Release write access.
 */
void synchronized_rw::write_unlock() const {
   _write_mutex.unlock();
   _read_mutex.unlock();
}

} /* namespace synchronizables */
} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */
