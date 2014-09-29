/*
 * Atomic values.
 * An atomic<T> is an atomic access control wrapper around a value of type T.
 */
#ifndef CONCURRENT__THREADS__SYNCHRONIZATION__ATOMIC_HH
#define CONCURRENT__THREADS__SYNCHRONIZATION__ATOMIC_HH

#include "concurrent/threads/synchronization/mutex.hh"

namespace concurrent {
namespace threads {
namespace synchronization {

/*
 * Atomic class.
 */
template <typename T>
class atomic {
public:
   /*
    * Constructors.
    */
   atomic();
   explicit atomic(const T&);
   
   /*
    * Copy constructor.
    */
   atomic(const atomic<T>&);

   /*
    * Destructor.
    */
   virtual ~atomic();

   /*
    * Get the value.
    */
   T get() const;

   /*
    * Set the value.
    * Return the prior value.
    */
   T set(const T&);

protected:
   T             _value;         /* contained value */
   mutable mutex _value_mutex;   /* mutex for reading/writing value */
};

/***************************************************************************
 * Atomic implementation.
 ***************************************************************************/

/*
 * Constructor.
 * Default initialize the contained value.
 */
template <typename T>
atomic<T>::atomic()
 : _value(),
   _value_mutex()
{ }

/*
 * Constructor.
 * Initialize with the specified value.
 */
template <typename T>
atomic<T>::atomic(const T& t)
 : _value(t),
   _value_mutex()
{ }

/*
 * Copy constructor.
 */
template <typename T>
atomic<T>::atomic(const atomic<T>& a)
 : _value(a.get()),
   _value_mutex()
{ }

/*
 * Destructor.
 */
template <typename T>
atomic<T>::~atomic() {
   /* do nothing */
}

/*
 * Get the value.
 */
template <typename T>
T atomic<T>::get() const {
   _value_mutex.lock();
   T t(_value);
   _value_mutex.unlock();
   return t;
}

/*
 * Set the value.
 * Return the prior value.
 */
template <typename T>
T atomic<T>::set(const T& t) {
   _value_mutex.lock();
   T t_prev(_value);
   _value = t;
   _value_mutex.unlock();
   return t_prev;
}

} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */

#endif
