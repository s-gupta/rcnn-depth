/*
 * Thread.
 */
#include "concurrent/threads/runnable.hh"
#include "concurrent/threads/thread.hh"
#include "lang/null.hh"

#include <pthread.h>
#include <time.h>
#include <unistd.h>

namespace concurrent {
namespace threads {

/*
 * Default constructor.
 * The thread itself provides the run method.
 */
thread::thread()
 : _runnable_obj(*this),
   _join_lock(),
   _n_processors(0)
{ }

/*
 * Constructor.
 * The specified runnable object provides the run method.
 */
thread::thread(runnable& r)
 : _runnable_obj(r),
   _join_lock(),
   _n_processors(0)
{ }

/*
 * Private copy constructor.
 * Threads should not be copied.
 */
thread::thread(const thread& t)
 : _runnable_obj(t._runnable_obj),
   _join_lock(),
   _n_processors(0)
{ }

/*
 * Destructor.
 */
thread::~thread() {
   /* do nothing */
}

/*
 * Default run method (does nothing).
 */
void thread::run() { }

/*
 * Start thread execution.
 */
void thread::start() {
   pthread_attr_t thread_attr;
   pthread_t thread_id;
   /* initialize thread attributes */
   pthread_attr_init(&thread_attr);
   pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_DETACHED);
   /* acquire the join lock */
   _join_lock.lock();
   /* inherit all available processors from the parent thread */
   _n_processors = thread::processors();
   /* create new pthread */
   pthread_create(
      &thread_id,
      &thread_attr,
      thread::thread_run,
      static_cast<void*>(this)
   );
   /* destroy attributes structure */
   pthread_attr_destroy(&thread_attr);
}

/*
 * Join the thread (wait until it terminates).
 * Return immediately if the thread is not running.
 */
void thread::join() {
   _join_lock.lock();
   _join_lock.unlock();
}

/*
 * Helper function for thread execution.
 */ 
void* thread::thread_run(void* t) {
   thread* thread_obj = static_cast<thread*>(t);
   /* set key for # processors available to the current thread */
   thread::n_processors_key().set(thread_obj->_n_processors);
   /* run the thread */
   thread_obj->_runnable_obj.run();
   /* clear key for # processors available to the current thread */
   thread::n_processors_key().clear();
   /* allow waiting threads to join */
   thread_obj->_join_lock.unlock();
   return NULL;
}

/*
 * Sleep the current thread for the specified number of seconds.
 * The sleep may be interrupted by a signal.
 */
void thread::sleep_sec(unsigned long n) {
   struct timespec req;
   struct timespec rem;
   req.tv_sec  = static_cast<time_t>(n);
   req.tv_nsec = 0;
   nanosleep(&req, &rem);
}

/*
 * Sleep the current thread for the specified number of microseconds.
 * The sleep may be interrupted by a signal.
 */
void thread::sleep_usec(unsigned long n) {
   struct timespec req;
   struct timespec rem;
   long sec    = (static_cast<long>(n))/(long(1000000));
   req.tv_sec  = static_cast<time_t>(sec);
   req.tv_nsec = (static_cast<long>(n) - (sec*1000000))*1000;
   nanosleep(&req, &rem);
}

/*
 * Sleep the current thread for the specified number of nanoseconds.
 * The sleep may be interrupted by a signal.
 */
void thread::sleep_nsec(unsigned long n) {
   struct timespec req;
   struct timespec rem;
   long sec    = (static_cast<long>(n))/(long(1000000000));
   req.tv_sec  = static_cast<time_t>(sec);
   req.tv_nsec = (static_cast<long>(n) - (sec*1000000000));
   nanosleep(&req, &rem);
}

/*
 * Return the number of processors available on the system.
 * Programs should treat this number as an estimate of the amount 
 * of thread-level parallelism they may attempt to utilize.
 */
unsigned long thread::processors() {
   return thread::n_processors_key().get();
}

/*
 * Set the number of processors available on the system.
 * This may be more or less than the actual number of processors.
 * An argument of zero resets the processor count to the actual 
 * number of processors on the machine.  Return the new count.
 */
unsigned long thread::processors(unsigned long n) {
   if (n == 0)
      n = 1;
   thread::n_processors_key().get() = n;
   return n;
}

/*
 * Return the key for the available processor count.
 */
thread_key<unsigned long>& thread::n_processors_key() {
   static thread_key<unsigned long>* n_procs_key 
      = new thread_key<unsigned long>(thread::n_processors_main());
   return *n_procs_key;
}

/*
 * Return the number of processors available to the main program.
 * Initialize the processor count to the actual number of processors.
 */
unsigned long& thread::n_processors_main() {
   static unsigned long n_procs_main 
      = 1;
   return n_procs_main;
}

} /* namespace threads */
} /* namespace concurrent */
