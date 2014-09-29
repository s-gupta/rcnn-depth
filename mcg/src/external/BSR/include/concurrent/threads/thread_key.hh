/*
 * Thread key.
 * 
 * A thead key can be used to create variables that have different values in 
 * different threads.
 */
#ifndef CONCURRENT__THREADS__THREAD_KEY_HH
#define CONCURRENT__THREADS__THREAD_KEY_HH

#include "lang/exceptions/ex_not_found.hh"
#include "lang/null.hh"

#include <pthread.h>

namespace concurrent {
namespace threads {
/*
 * Imports.
 */
using lang::exceptions::ex_not_found;

/*
 * Thread key.
 */
template <typename T>
class thread_key {
public:
   /*
    * Constructor.
    * Allocate a new key that can be used to associate each thread with an 
    * object of type T.
    */
   thread_key();
   
   /*
    * Constructor.
    * Allocate a new key that can be used to associate each thread with an 
    * object of type T and set the associated object for the current thread.
    */
   explicit thread_key(T&);

   /*
    * Copy constructor.
    * Allocate a new key and set the associated object for the current thread 
    * to the same object associated with the given key.
    */
   explicit thread_key(const thread_key<T>&);

   /*
    * Destructor.
    * Deallocate the key.
    */
   virtual ~thread_key();
   
   /*
    * Set the object associated with the current thread.
    */
   void set(T&);

   /*
    * Get the object associated with the current thread.
    * Throw an exception (ex_not_found) if no object is associated with the 
    * current thread.
    */
   T& get() const;

   /*
    * Remove the object (if any) associated with the current thread.
    * This does not delete the object, only removes the association.
    */
   void clear();

protected:
   pthread_key_t _key;  /* key for thread-specific data */
};

/*
 * Constructor.
 * Allocate a new key that can be used to associate each thread with an 
 * object of type T.
 */
template <typename T>
thread_key<T>::thread_key() {
   pthread_key_create(&_key, NULL);
}

/*
 * Constructor.
 * Allocate a new key that can be used to associate each thread with an 
 * object of type T and set the associated object for the current thread.
 */
template <typename T>
thread_key<T>::thread_key(T& t) {
   pthread_key_create(&_key, NULL);
   this->set(t);
}

/*
 * Copy constructor.
 * Allocate a new key and set the associated object for the current thread 
 * to the same object associated with the given key.
 */
template <typename T>
thread_key<T>::thread_key(const thread_key<T>& t_key) {
   T& t = t_key.get();
   pthread_key_create(&_key, NULL);
   this->set(t);
}

/*
 * Destructor.
 * Deallocate the key.
 */
template <typename T>
thread_key<T>::~thread_key() {
   pthread_key_delete(_key);
}

/*
 * Set the object associated with the current thread.
 */
template <typename T>
void thread_key<T>::set(T& t) {
   pthread_setspecific(_key, &t);
}

/*
 * Get the object associated with the current thread.
 * Throw an exception (ex_not_found) if no object is associated with the 
 * current thread.
 */
template <typename T>
T& thread_key<T>::get() const {
   T* p = static_cast<T*>(pthread_getspecific(_key));
   if (p == NULL)
      throw ex_not_found(
         "no object associated with key in current thread"
      );
   return *p;
}

/*
 * Remove the object (if any) associated with the current thread.
 * This does not delete the object, only removes the association.
 */
template <typename T>
void thread_key<T>::clear() {
   pthread_setspecific(_key, NULL);
}

} /* namespace threads */
} /* namespace concurrent */

#endif
