/*
 * An auto_read_write_lock simultaneously takes a read lock on one
 * synchronizable object and a write lock on another when constructed 
 * and releases the locks when destructed.  The order in which the 
 * locks are taken is handled automatically.
 */
#ifndef CONCURRENT__THREADS__SYNCHRONIZATION__LOCKS__AUTO_READ_WRITE_LOCK_HH
#define CONCURRENT__THREADS__SYNCHRONIZATION__LOCKS__AUTO_READ_WRITE_LOCK_HH

namespace concurrent {
namespace threads {
namespace synchronization {
namespace locks {

/*
 * Auto read write lock.
 */
template <typename Syn>
class auto_read_write_lock {
public:
   /*
    * Constructor.
    * Acquire read lock on the first object and write lock on the second.
    */
   explicit inline auto_read_write_lock(Syn&, Syn&);
   

   /*
    * Destructor.
    * Release locks.
    */
   inline ~auto_read_write_lock();

protected:
   Syn& _sr;   /* synchronizable object on which to claim read lock  */
   Syn& _sw;   /* synchronizable object on which to claim write lock */

private:
   /*
    * Private copy constructor.
    * Auto read write locks should not be copied as they provide
    * exclusive access to the object locked for writing.
    */
   explicit inline auto_read_write_lock(const auto_read_write_lock&);
};

/*
 * Constructor.
 * Acquire read lock on the first object and write lock on the second.
 */
template <typename Syn>
inline auto_read_write_lock<Syn>::auto_read_write_lock(Syn& sr, Syn& sw)
 : _sr(sr),
   _sw(sw)
{
   Syn::read_write_lock(_sr, _sw);
}

/*
 * Private copy constructor.
 * Auto read write locks should not be copied as they provide
 * exclusive access to the object locked for writing.
 */
template <typename Syn>
inline auto_read_write_lock<Syn>::auto_read_write_lock(
   const auto_read_write_lock& l)
 : _sr(l._sr),
   _sw(l._sw)
{
   Syn::read_write_lock(_sr, _sw);
}

/*
 * Destructor.
 * Release locks.
 */
template <typename Syn>
inline auto_read_write_lock<Syn>::~auto_read_write_lock() {
   Syn::read_write_unlock(_sr, _sw);
}

} /* namespace locks */
} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */

#endif
