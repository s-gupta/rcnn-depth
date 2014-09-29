/*
 * Auto write lock.
 *
 * An auto_write_lock takes a write lock on a synchronizable object when
 * constructed and releases the write lock when destructed.
 */
#ifndef CONCURRENT__THREADS__SYNCHRONIZATION__LOCKS__AUTO_WRITE_LOCK_HH
#define CONCURRENT__THREADS__SYNCHRONIZATION__LOCKS__AUTO_WRITE_LOCK_HH

namespace concurrent {
namespace threads {
namespace synchronization {
namespace locks {

/*
 * Auto write lock.
 */
template <typename Syn>
class auto_write_lock {
public:
   /*
    * Constructor.
    * Acquire write lock on the given object.
    */
   explicit inline auto_write_lock(Syn&);

   /*
    * Destructor.
    * Release write lock.
    */
   inline ~auto_write_lock();

protected:
   Syn& _s;    /* synchronizable object */
   
private:
   /*
    * Private copy constructor.
    * Auto write locks should not be copied as they provide exclusive access
    * to the locked object.
    */
   explicit inline auto_write_lock(const auto_write_lock&);
};

/*
 * Constructor.
 * Acquire write lock on the given object.
 */
template <typename Syn>
inline auto_write_lock<Syn>::auto_write_lock(Syn& s)
 : _s(s) 
{
   _s.write_lock();
}

/*
 * Private copy constructor.
 * Auto write locks should not be copied as they provide exclusive access
 * to the locked object.
 */
template <typename Syn>
inline auto_write_lock<Syn>::auto_write_lock(const auto_write_lock& l)
 : _s(l._s) 
{
   _s.write_lock();
}

/*
 * Destructor.
 * Release write lock.
 */
template <typename Syn>
inline auto_write_lock<Syn>::~auto_write_lock() {
   _s.write_unlock();
}

} /* namespace locks */
} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */

#endif
