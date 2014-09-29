/*
 * Child thread.
 *
 * The child thread class provides means for creating threads in which an 
 * uncaught exception does not terminate the program, but is instead rethrown 
 * when the child thread is joined.
 *
 * A parent thread does not execute while it has running child threads.  
 * Once starting a collection of child threads, the parent thread waits 
 * for all of them to complete before continuing.
 *
 * In addition, by default, the parent thread's available processors are 
 * split evenly among its child threads (with each thread getting one 
 * processor if there are more threads than processors).  This behavior 
 * can be changed by calling thread::processors(...) during execution of 
 * a child thread.
 *
 * In the case that mutiple child threads throw uncaught exceptions, 
 * the exception rethrown is that belonging to the child thread appearing 
 * earliest in the collection of child threads.
 *
 * NOTE: Child threads will only prevent uncaught exceptions derived from 
 * the throwable class from terminating the program.
 */
#ifndef CONCURRENT__THREADS__CHILD_THREAD_HH
#define CONCURRENT__THREADS__CHILD_THREAD_HH

#include "collections/abstract/collection.hh"
#include "concurrent/threads/runnable.hh"
#include "lang/exceptions/throwable.hh"
#include "lang/pointers/auto_ptr.hh"

namespace concurrent {
namespace threads {
/*
 * Imports.
 */
using collections::abstract::collection;
using lang::exceptions::throwable;
using lang::pointers::auto_ptr;

/*
 * Child thread.
 */
class child_thread : public runnable {
public:
   /*
    * Constructors.
    * Associate a runnable object (defaults to the child thread itself) with 
    * the child thread.
    */
   child_thread();
   explicit child_thread(runnable&);

   /*
    * Destructor.
    */
   virtual ~child_thread();

   /*
    * Default run method.
    */
   virtual void run();
   
   /*
    * Start (in new execution threads) and run to completion both child 
    * threads.  The child threads execute in parallel.
    *
    * This function call blocks until both child threads have finished.
    *
    * An exception thrown but not handled by one child thread is rethrown
    * upon completion of both child threads.  In the case that both 
    * child threads throw uncaught exceptions, the exception rethrown 
    * is that belonging to the first child thread.
    */
   static void start(child_thread&, child_thread&);

   /*
    * Start (in new execution threads) and run to completion all child 
    * threads in the collection.  The child threads execute in parallel.
    *
    * This function call blocks until all child threads have finished.
    *
    * An exception thrown but not handled by a child thread is rethrown
    * upon completion of all child threads.  In the case that multiple 
    * child threads throw uncaught exceptions, the exception rethrown 
    * is that belonging to the child thread appearing earliest in the 
    * collection.
    */
   static void start(const collection<child_thread>&);

   /*
    * Create and start child threads for the given runnable objects.
    */
   static void start(runnable&, runnable&);
   static void start(const collection<runnable>&);
   
   /*
    * Run (either sequentially or in parallel) both child threads.
    *
    * If thread::processors() indicates that more than one processor is 
    * available, the start(...) function above is called to execute the child 
    * threads in parallel and rethrow any uncaught exception upon completion 
    * of both child threads.
    *
    * If only one processor is available, the child threads are executed 
    * sequentially.  Note that in this case exceptions are thrown immediately 
    * rather than held until completion of both child threads.
    */
   static void run(child_thread&, child_thread&);

   /*
    * Run (either sequentially or in parallel) all child threads.
    *
    * If thread::processors() indicates that at least as many processors are 
    * available as child threads in the collection, then the start(...) 
    * function above is called to execute all child threads in parallel and 
    * and rethrow any uncaught exception upon completion of all child threads.
    *
    * Otherwise, the threads are divided into exactly thread::processors() 
    * groups, the groups are run in parallel, and the threads within each 
    * group are run sequentially.  Note that in this case exceptions are 
    * thrown immediately within each group rather than held until completion 
    * of all child threads within the group.
    */
   static void run(const collection<child_thread>&);

   /*
    * Create and run child threads for the given runnable objects.
    */
   static void run(runnable&, runnable&);
   static void run(const collection<runnable>&);
   
private:
   /*
    * Private copy constructor.
    * Child threads should not be copied.
    */
   explicit child_thread(const child_thread&);
   
   /*
    * Child thread data.
    */
   runnable& _runnable_obj;  /* runnable interface */

   /*
    * Runnable object for packaging together multiple runnable objects 
    * to be executed sequentially upon calling the group's run() method.
    */
   class runnable_group : public runnable {
   public:
      /*
       * Constructor.
       * Create an empty runnable group.
       */
      runnable_group();

      /*
       * Copy constructor.
       */
      runnable_group(const runnable_group&);
      
      /*
       * Destructor.
       */
      virtual ~runnable_group();

      /*
       * Run method.
       * Run each member of the group sequentially.
       */
      virtual void run();
      
      /*
       * Group members.
       */
      collection<runnable>& runnables;
   };
   
   /*
    * Runnable object which wraps the run() method of the child thread, 
    * setting the default number of processors, and catching exceptions.
    */
   class child_runnable : public runnable {
   public:
      /*
       * Constructor.
       */
      explicit child_runnable(
         runnable&,     /* object to run */
         unsigned long  /* # processors to allocate to it */
      );

      /*
       * Destructor.
       */
      virtual ~child_runnable();

      /*
       * Run method.
       * Wrap the run method of the thread.
       */
      virtual void run();
     
      /*
       * Raise pending exception (if any).
       */
      void raise();
      
   private:
      runnable&           _runnable_obj;  /* object to run */
      unsigned long       _n_processors;  /* # processors to allocate to child */
      auto_ptr<throwable> _exception;     /* pending exception (if any) */
   };
};

} /* namespace threads */
} /* namespace concurrent */

#endif
