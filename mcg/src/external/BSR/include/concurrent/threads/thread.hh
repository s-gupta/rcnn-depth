/*
 * Thread.
 */
#ifndef CONCURRENT__THREADS__THREAD_HH
#define CONCURRENT__THREADS__THREAD_HH

#include "concurrent/threads/synchronization/mutex.hh"
#include "concurrent/threads/runnable.hh"
#include "concurrent/threads/thread_key.hh"

namespace concurrent {
namespace threads {
/*
 * Imports.
 */
using concurrent::threads::synchronization::mutex;

/*
 * Thread.
 */
class thread : public runnable {
public:
   /*
    * Constructors.
    * Associate a runnable object (defaults to the thread itself) with 
    * the thread.
    */
   thread();
   explicit thread(runnable&);

   /*
    * Destructor.
    */
   virtual ~thread();
    
   /*
    * Default run method.
    */
   virtual void run();

   /*
    * Start thread execution.
    * The run() method of the associated runnable object is executed in a 
    * new thread.
    */
   void start();
   
   /*
    * Join the thread (wait until it terminates).
    * Return immediately if the thread is not running.
    */
   void join();

   /*
    * Sleep the current thread for the specified number of seconds.
    * The sleep may be interrupted by a signal.
    */
   static void sleep_sec(unsigned long);

   /*
    * Sleep the current thread for the specified number of microseconds.
    * The sleep may be interrupted by a signal.
    */
   static void sleep_usec(unsigned long);

   /*
    * Sleep the current thread for the specified number of nanoseconds.
    * The sleep may be interrupted by a signal.
    */
   static void sleep_nsec(unsigned long);

   /*
    * Return the number of processors available to the current thread.
    * Programs should treat this number as an estimate of the amount 
    * of thread-level parallelism they may attempt to utilize.
    */
   static unsigned long processors();

   /*
    * Set the number of processors available to the current thread.
    * This may be more or less than the actual number of processors.
    * An argument of zero resets the processor count to the actual 
    * number of processors on the machine.  Return the new count.
    */
   static unsigned long processors(unsigned long);

private:
   /*
    * Private copy constructor.
    * Threads should not be copied.
    */
   explicit thread(const thread&);

   /*
    * Thread data.
    */
   runnable&     _runnable_obj;  /* runnable interface */
   mutex         _join_lock;     /* lock for join() */
   unsigned long _n_processors;  /* # processors available to thread */

   /*
    * Number of processors available to each thread.
    */
   static thread_key<unsigned long>& n_processors_key(); /* current thread */
   static unsigned long& n_processors_main();            /* main thread */
   
   /*
    * Helper function for thread execution.
    */ 
   static void* thread_run(void*);
};

} /* namespace threads */
} /* namespace concurrent */

#endif
