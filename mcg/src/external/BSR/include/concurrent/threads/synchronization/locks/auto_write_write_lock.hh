/*
 * An auto_write_write_lock simultaneously takes write locks on two 
 * synchronizable objects when constructed and releases the locks 
 * when destructed.  The order in which the locks are taken is 
 * handled automatically.
 */
#ifndef CONCURRENT__THREADS__SYNCHRONIZATION__LOCKS__AUTO_WRITE_WRITE_LOCK_HH
#define CONCURRENT__THREADS__SYNCHRONIZATION__LOCKS__AUTO_WRITE_WRITE_LOCK_HH

namespace concurrent {
namespace threads {
namespace synchronization {
namespace locks {

/*
 * Auto write write lock.
 */
template <typename Syn>
class auto_write_write_lock {
public:
   /*
    * Constructor.
    * Acquire write locks on the given objects.
    */
   explicit inline auto_write_write_lock(Syn&, Syn&);

   /*
    * Destructor.
    * Release write locks.
    */
   inline ~auto_write_write_lock();

protected:
   Syn& _s0;   /* synchronizable objects */
   Syn& _s1;   

private:
   /*
    * Private copy constructor.
    * Auto write write locks should not be copied as they provide exclusive
    * access to the objects locked for writing.
    */
   explicit inline auto_write_write_lock(const auto_write_write_lock&);
};

/*
 * Constructor.
 * Acquire write locks on the given objects.
 */
template <typename Syn>
inline auto_write_write_lock<Syn>::auto_write_write_lock(Syn& s0, Syn& s1)
 : _s0(s0),
   _s1(s1)
{
   Syn::write_lock(_s0, _s1);
}

/*
 * Private copy constructor.
 * Auto write write locks should not be copied as they provide exclusive
 * access to the objects locked for writing.
 */
template <typename Syn>
inline auto_write_write_lock<Syn>::auto_write_write_lock(
   const auto_write_write_lock& l)
 : _s0(l._s0),
   _s1(l._s1)
{
   Syn::write_lock(_s0, _s1);
}

/*
 * Destructor.
 * Release write locks.
 */
template <typename Syn>
inline auto_write_write_lock<Syn>::~auto_write_write_lock() {
   Syn::write_unlock(_s0, _s1);
}

} /* namespace locks */
} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */

#endif
