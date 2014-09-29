/*
 * Auto lock.
 *
 * An auto_lock takes a lock on a mutex or a synchronizable object when
 * constructed and releases the lock when destructed.
 */
#ifndef CONCURRENT__THREADS__SYNCHRONIZATION__LOCKS__AUTO_LOCK_HH
#define CONCURRENT__THREADS__SYNCHRONIZATION__LOCKS__AUTO_LOCK_HH

namespace concurrent {
namespace threads {
namespace synchronization {
namespace locks {

/*
 * Auto lock.
 */
template <typename Syn>
class auto_lock {
public:
   /*
    * Constructor.
    * Acquire lock on the given object.
    */
   explicit inline auto_lock(Syn&);

   /*
    * Destructor.
    * Release lock.
    */
   inline ~auto_lock();

protected:
   Syn& _s;    /* synchronizable object */
   
private:
   /*
    * Private copy constructor.
    * Auto locks should not be copied as they provide exclusive access to the
    * locked object.
    */
   explicit inline auto_lock(const auto_lock&);
};

/*
 * Constructor.
 * Acquire lock on the given object.
 */
template <typename Syn>
inline auto_lock<Syn>::auto_lock(Syn& s)
 : _s(s) 
{
   _s.lock();
}

/*
 * Private copy constructor.
 * Auto locks should not be copied as they provide exclusive access to the
 * locked object.
 */
template <typename Syn>
inline auto_lock<Syn>::auto_lock(const auto_lock& l)
 : _s(l._s) 
{
   _s.lock();
}

/*
 * Destructor.
 * Release lock.
 */
template <typename Syn>
inline auto_lock<Syn>::~auto_lock() {
   _s.unlock();
}

} /* namespace locks */
} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */

#endif
