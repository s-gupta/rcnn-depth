/*
 * An auto_read_read_lock simultaneously takes read locks on two synchronizable
 * objects when constructed and releases the locks when destructed.  The order
 * in which the locks are taken is handled automatically.
 */
#ifndef CONCURRENT__THREADS__SYNCHRONIZATION__LOCKS__AUTO_READ_READ_LOCK_HH
#define CONCURRENT__THREADS__SYNCHRONIZATION__LOCKS__AUTO_READ_READ_LOCK_HH

namespace concurrent {
namespace threads {
namespace synchronization {
namespace locks {

/*
 * Auto read read lock.
 */
template <typename Syn>
class auto_read_read_lock {
public:
   /*
    * Constructor.
    * Acquire read locks on the given objects.
    */
   explicit inline auto_read_read_lock(Syn&, Syn&);
   
   /*
    * Copy constructor.
    * Copying the lock invokes another set of read locks on the objects.
    */
   explicit inline auto_read_read_lock(const auto_read_read_lock&);

   /*
    * Destructor.
    * Release read locks.
    */
   inline ~auto_read_read_lock();

protected:
   Syn& _s0;   /* synchronizable objects */
   Syn& _s1;
};

/*
 * Constructor.
 * Acquire read locks on the given objects.
 */
template <typename Syn>
inline auto_read_read_lock<Syn>::auto_read_read_lock(Syn& s0, Syn& s1)
 : _s0(s0),
   _s1(s1)
{
   Syn::read_lock(_s0, _s1);
}

/*
 * Copy constructor.
 * Copying the lock invokes another set of read locks on the objects.
 */
template <typename Syn>
inline auto_read_read_lock<Syn>::auto_read_read_lock(
   const auto_read_read_lock& l)
 : _s0(l._s0),
   _s1(l._s1)
{
   Syn::read_lock(_s0, _s1);
}

/*
 * Destructor.
 * Release read locks.
 */
template <typename Syn>
inline auto_read_read_lock<Syn>::~auto_read_read_lock() {
   Syn::read_unlock(_s0, _s1);
}

} /* namespace locks */
} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */

#endif
