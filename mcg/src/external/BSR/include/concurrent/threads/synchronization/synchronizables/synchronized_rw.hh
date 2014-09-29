/*
 * The synchronized_rw class implements access control locks according to a 
 * fair readers/writers algorithm.  
 *
 * Any number of readers may be active at once, but no readers and no other 
 * writers are permitted when one writer is active.
 *
 * Access is granted on a first-come first-serve basis.
 */
#ifndef CONCURRENT__THREADS__SYNCHRONIZATION__SYNCHRONIZABLES__SYNCHRONIZED_RW_HH
#define CONCURRENT__THREADS__SYNCHRONIZATION__SYNCHRONIZABLES__SYNCHRONIZED_RW_HH

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

class synchronized_rw : public synchronizable {
public:
   /*
    * Constructor.
    */
   synchronized_rw();
  
   /*
    * Destructor.
    */
   virtual ~synchronized_rw();

   /*
    * Claim exclusive access.
    */
   virtual void abstract_lock() const;
   void lock() const;
   
   /*
    * Release exclusive access.
    */
   virtual void abstract_unlock() const;
   void unlock() const;

   /*
    * Claim read access.
    */
   virtual void abstract_read_lock() const;
   void read_lock() const;

   /*
    * Release read access.
    */
   virtual void abstract_read_unlock() const;
   void read_unlock() const;

   /*
    * Claim write access.
    */
   virtual void abstract_write_lock() const;
   void write_lock() const;

   /*
    * Release write access.
    */
   virtual void abstract_write_unlock() const;
   void write_unlock() const;

private:
   /*
    * Private copy constructor.
    * (synchronized_rw objects should not be copied)
    */
   explicit synchronized_rw(const synchronized_rw&);
   
   /*
    * Access control data.
    */
   mutable unsigned long _n_readers;   /* number of active readers */
   mutable mutex _rcnt_mutex;          /* lock for counting # readers */
   mutable mutex _read_mutex;          /* lock for read capability */
   mutable mutex _write_mutex;         /* lock for write capability */
};

} /* namespace synchronizables */
} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */

#endif
