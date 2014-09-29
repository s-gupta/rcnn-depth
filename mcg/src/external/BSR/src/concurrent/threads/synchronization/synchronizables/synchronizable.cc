/*
 * Synchronizable interface.
 */
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

/*
 * Issue the next available identity for a synchronizable object.
 * The zero id is reserved for unsynchronized objects.
 */
unsigned long synchronizable::id_next() {
   static unsigned long id = 1;
   static mutex id_mutex;
   id_mutex.lock();
   if (id == 0) { id = 1; }
   unsigned long issue_id = id;
   id++;
   id_mutex.unlock();
   return issue_id;
}

/*
 * Constructor.
 * Initialize the identity of the synchronizable object.
 */
synchronizable::synchronizable()
 : _id(synchronizable::id_next())
{ }

/*
 * Constructor.
 * Initialize the synchronizable object to have the given id.
 */
synchronizable::synchronizable(unsigned long id)
 : _id(id)
{ }

/*
 * Destructor.
 */
synchronizable::~synchronizable() {
   /* do nothing */
}

/*
 * Claim exclusive access to two synchronizable objects.
 */
void synchronizable::lock(
   const synchronizable& s0,
   const synchronizable& s1)
{
   if (s0.syn_id() < s1.syn_id()) {
      s0.abstract_lock();
      s1.abstract_lock();
   } else if (s0.syn_id() > s1.syn_id()) {
      s1.abstract_lock();
      s0.abstract_lock();
   } else {
      s0.abstract_lock();
   }
}
   
/*
 * Release exclusive access to two synchronizable objects.
 */
void synchronizable::unlock(
   const synchronizable& s0,
   const synchronizable& s1)
{
   s0.abstract_unlock();
   if (s0.syn_id() != s1.syn_id())
      s1.abstract_unlock();
}

/*
 * Claim read access to two synchronizable objects.
 */
void synchronizable::read_lock(
   const synchronizable& s0,
   const synchronizable& s1)
{
   if (s0.syn_id() < s1.syn_id()) {
      s0.abstract_read_lock();
      s1.abstract_read_lock();
   } else if (s0.syn_id() > s1.syn_id()) {
      s1.abstract_read_lock();
      s0.abstract_read_lock();
   } else {
      s0.abstract_read_lock();
   }
}

/*
 * Release read access to two synchronizable objects.
 */
void synchronizable::read_unlock(
   const synchronizable& s0,
   const synchronizable& s1)
{
   s0.abstract_read_unlock();
   if (s0.syn_id() != s1.syn_id())
      s1.abstract_read_unlock();
}

/*
 * Claim write access to two synchronizable objects.
 */
void synchronizable::write_lock(
   const synchronizable& s0,
   const synchronizable& s1)
{
   if (s0.syn_id() < s1.syn_id()) {
      s0.abstract_write_lock();
      s1.abstract_write_lock();
   } else if (s0.syn_id() > s1.syn_id()) {
      s1.abstract_write_lock();
      s0.abstract_write_lock();
   } else {
      s0.abstract_write_lock();
   }
}

/*
 * Release write access to two synchronizable objects.
 */
void synchronizable::write_unlock(
   const synchronizable& s0,
   const synchronizable& s1)
{
   s0.abstract_write_unlock();
   if (s0.syn_id() != s1.syn_id())
      s1.abstract_write_unlock();
}

/*
 * Claim read access to the first object and write access to the second.
 */
void synchronizable::read_write_lock(
   const synchronizable& s0,
   const synchronizable& s1)
{
   if (s0.syn_id() < s1.syn_id()) {
      s0.abstract_read_lock();
      s1.abstract_write_lock();
   } else if (s0.syn_id() > s1.syn_id()) {
      s1.abstract_write_lock();
      s0.abstract_read_lock();
   } else {
      s0.abstract_write_lock();
   }
   
}

/*
 * Release read access to the first object and write access to the second.
 */
void synchronizable::read_write_unlock(
   const synchronizable& s0,
   const synchronizable& s1)
{
   if (s0.syn_id() != s1.syn_id()) {
      s0.abstract_read_unlock();
      s1.abstract_write_unlock();
   } else {
      s0.abstract_write_unlock();
   }  
}

} /* namespace synchronizables */
} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */
