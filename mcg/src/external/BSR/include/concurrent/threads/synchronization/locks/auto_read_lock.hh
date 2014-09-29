/*
 * Auto read lock.
 *
 * An auto_read_lock takes a read lock on a synchronizable object when
 * constructed and releases the read lock when destructed.
 */
#ifndef CONCURRENT__THREADS__SYNCHRONIZATION__LOCKS__AUTO_READ_LOCK_HH
#define CONCURRENT__THREADS__SYNCHRONIZATION__LOCKS__AUTO_READ_LOCK_HH

namespace concurrent {
namespace threads {
namespace synchronization {
namespace locks {

/*
 * Auto read lock.
 */
template <typename Syn>
class auto_read_lock {
public:
   /*
    * Constructor.
    * Acquire read lock on the given object.
    */
   explicit inline auto_read_lock(Syn&);
   
   /*
    * Copy constructor.
    * Copying an auto read lock invokes another read lock on the object.
    */
   explicit inline auto_read_lock(const auto_read_lock&);

   /*
    * Destructor.
    * Release read lock.
    */
   inline ~auto_read_lock();

protected:
   Syn& _s;    /* synchronizable object */
};

/*
 * Constructor.
 * Acquire read lock on the given object.
 */
template <typename Syn>
inline auto_read_lock<Syn>::auto_read_lock(Syn& s)
 : _s(s) 
{
   _s.read_lock();
}

/*
 * Copy constructor.
 * Copying an auto read lock invokes another read lock on the object.
 */
template <typename Syn>
inline auto_read_lock<Syn>::auto_read_lock(const auto_read_lock& l)
 : _s(l._s) 
{
   _s.read_lock();
}

/*
 * Destructor.
 * Release read lock.
 */
template <typename Syn>
inline auto_read_lock<Syn>::~auto_read_lock() {
   _s.read_unlock();
}

} /* namespace locks */
} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */

#endif
