/*
 * Mutex.
 */
#ifndef CONCURRENT__THREADS__SYNCHRONIZATION__MUTEX_HH
#define CONCURRENT__THREADS__SYNCHRONIZATION__MUTEX_HH

#include <pthread.h>

namespace concurrent { 
namespace threads {
namespace synchronization {

class mutex {
public:
   /*
    * Constructor.
    * Initialize the mutex to an unlocked state.
    */
   inline mutex();

   /*
    * Destructor.
    */
   inline ~mutex();

   /*
    * Lock the mutex.
    */
   inline void lock();

   /*
    * Unlock the mutex.
    */
   inline void unlock();

protected:
   /*
    * Private copy constructor.
    * A mutex should not be copied.
    */
   explicit inline mutex(const mutex&);

   pthread_mutex_t _m;  /* mutex */
};

/*
 * Constructor.
 * Initialize the mutex to an unlocked state.
 */
inline mutex::mutex() {
   pthread_mutex_init(&_m, NULL);
}

/*
 * Private copy constructor.
 * A mutex should not be copied.
 */
inline mutex::mutex(const mutex&) {
   pthread_mutex_init(&_m, NULL);
}

/*
 * Destructor.
 */
inline mutex::~mutex() {
   pthread_mutex_destroy(&_m);
}

/*
 * Lock the mutex.
 */
inline void mutex::lock() {
   pthread_mutex_lock(&_m);
}

/*
 * Unlock the mutex.
 */
inline void mutex::unlock() {
   pthread_mutex_unlock(&_m);
}

} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */

#endif
