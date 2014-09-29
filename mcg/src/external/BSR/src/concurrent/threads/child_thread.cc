/*
 * Child thread.
 */
#include "collections/abstract/collection.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "concurrent/threads/child_thread.hh"
#include "concurrent/threads/runnable.hh"
#include "concurrent/threads/thread.hh"
#include "lang/exceptions/throwable.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"

namespace concurrent {
namespace threads {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::list;
using collections::pointers::auto_collection;
using lang::exceptions::throwable;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

/*
 * Default constructor.
 * The child thread itself provides the run method.
 */
child_thread::child_thread()
 : _runnable_obj(*this)
{ }

/*
 * Constructor.
 * The specified runnable object provides the run method.
 */
child_thread::child_thread(runnable& r)
 : _runnable_obj(r)
{ }

/*
 * Private copy constructor.
 * Child threads should not be copied.
 */
child_thread::child_thread(const child_thread& t)
 : _runnable_obj(t._runnable_obj)
{ }

/*
 * Destructor.
 */
child_thread::~child_thread() {
   /* do nothing */
}

/*
 * Default run method (does nothing).
 */
void child_thread::run() { }

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
void child_thread::start(child_thread& t0, child_thread& t1) {
   /* split the available processors between child threads */
   unsigned long n_proc = thread::processors();
   unsigned long n_proc0 = 1;
   unsigned long n_proc1 = 1;
   if (n_proc > 2) {
      n_proc0 = n_proc/2;
      n_proc1 = n_proc - n_proc0;
   }
   /* wrap the run methods of the child threads */
   child_runnable r0(t0._runnable_obj, n_proc0);
   child_runnable r1(t1._runnable_obj, n_proc1);
   /* create thread objects */
   thread thrd0(r0);
   thread thrd1(r1);
   /* start both threads */
   thrd0.start();
   thrd1.start();
   /* join both threads */
   thrd0.join();
   thrd1.join();
   /* rethrow first uncaught exception */
   r0.raise();
   r1.raise();
}

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
void child_thread::start(const collection<child_thread>& c) {
   /* check that collection is nonempty */
   unsigned long n_threads = c.size();
   if (n_threads > 0) {
      /* split the available processors between child threads */
      unsigned long n_procs = thread::processors();
      unsigned long n_procs_per_thread = 1;
      unsigned long n_procs_extra = 0;
      if (n_procs > n_threads) {
         n_procs_per_thread = n_procs/n_threads;
         n_procs_extra = n_procs - (n_procs_per_thread*n_threads);
      }
      /* wrap run methods of the child threads */
      auto_collection< child_runnable, list<child_runnable> > r_list(
         new list<child_runnable>()
      );
      auto_ptr< iterator<child_thread> > i_thread = c.iter_create();
      for (unsigned long n = 0; n < n_procs_extra; n++) {
         auto_ptr<child_runnable> r(
            new child_runnable(
               i_thread->next()._runnable_obj, n_procs_per_thread + 1
            )
         );
         r_list->add(*r);
         r.release();
      }
      while (i_thread->has_next()) {
         auto_ptr<child_runnable> r(
            new child_runnable(
               i_thread->next()._runnable_obj, n_procs_per_thread
            )
         );
         r_list->add(*r);
         r.release();
      }
      /* create thread objects */
      auto_collection< thread, list<thread> > t_list(new list<thread>());
      for (list<child_runnable>::iterator_t i(*r_list); i.has_next(); ) {
         auto_ptr<thread> t(new thread(i.next()));
         t_list->add(*t);
         t.release();
      }
      /* start all threads */
      for (list<thread>::iterator_t i(*t_list); i.has_next(); )
         i.next().start();
      /* join all threads */
      for (list<thread>::iterator_t i(*t_list); i.has_next(); )
         i.next().join();
      /* rethrow first uncaught exception */
      for (list<child_runnable>::iterator_t i(*r_list); i.has_next(); )
         i.next().raise();
   }
}

/*
 * Create and start child threads for the given runnable objects.
 */
void child_thread::start(runnable& r0, runnable& r1) {
   child_thread t0(r0);
   child_thread t1(r1);
   child_thread::start(t0, t1);
}

void child_thread::start(const collection<runnable>& c) {
   auto_collection< child_thread, list<child_thread> > t_list(
      new list<child_thread>()
   );
   auto_ptr< iterator<runnable> > i = c.iter_create();
   while (i->has_next()) {
      auto_ptr<child_thread> t(new child_thread(i->next()));
      t_list->add(*t);
      t.release();
   }
   child_thread::start(*t_list);
}

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
void child_thread::run(child_thread& t0, child_thread& t1) {
   if (thread::processors() > 1) {
      child_thread::start(t0, t1);
   } else {
      t0._runnable_obj.run();
      t1._runnable_obj.run();
   }
}

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
void child_thread::run(const collection<child_thread>& c) {
   unsigned long n_threads = c.size();
   unsigned long n_procs = thread::processors();
   if (n_procs >= n_threads) {
      /* start all threads in parallel */
      child_thread::start(c);
   } else {
      /* compute # of threads assigned to each processor */
      unsigned long n_threads_per_proc = n_threads/n_procs;
      unsigned long n_threads_extra = n_threads - (n_threads_per_proc*n_procs);
      /* initialize runnable groups */
      auto_collection< runnable, list<runnable> > r_group_list(
         new list<runnable>()
      );
      /* place threads into runnable groups */
      auto_ptr< iterator<child_thread> > i = c.iter_create();
      for (unsigned long n_group = 0; n_group < n_threads_extra; n_group++) {
         /* create runnable group that includes an extra thread */
         auto_ptr<runnable_group> r_group(new runnable_group());
         for (unsigned long n = 0; n <= n_threads_per_proc; n++)
            r_group->runnables.add(i->next()._runnable_obj);
         r_group_list->add(*r_group);
         r_group.release();
      }
      while (i->has_next()) {
         /* create runnable group */
         auto_ptr<runnable_group> r_group(new runnable_group());
         for (unsigned long n = 0; n < n_threads_per_proc; n++)
            r_group->runnables.add(i->next()._runnable_obj);
         r_group_list->add(*r_group);
         r_group.release();
      }
      /* start runnable groups in parallel */
      child_thread::start(*r_group_list);
   }
}

/*
 * Create and run child threads for the given runnable objects.
 */
void child_thread::run(runnable& r0, runnable& r1) {
   child_thread t0(r0);
   child_thread t1(r1);
   child_thread::run(t0, t1);
}

void child_thread::run(const collection<runnable>& c) {
   auto_collection< child_thread, list<child_thread> > t_list(
      new list<child_thread>()
   );
   auto_ptr< iterator<runnable> > i = c.iter_create();
   while (i->has_next()) {
      auto_ptr<child_thread> t(new child_thread(i->next()));
      t_list->add(*t);
      t.release();
   }
   child_thread::run(*t_list);
}

/***************************************************************************
 * Runnable group implementation.
 ***************************************************************************/

/*
 * Constructor.
 * Create an empty runnable group.
 */
child_thread::runnable_group::runnable_group()
 : runnables(*(new list<runnable>()))
{ }

/*
 * Copy constructor.
 */
child_thread::runnable_group::runnable_group(const runnable_group& r)
 : runnables(*(new list<runnable>(r.runnables)))
{ }

/*
 * Destructor.
 */
child_thread::runnable_group::~runnable_group() {
   delete &runnables;
}

/*
 * Run method.
 * Run each member of the group sequentially.
 */
void child_thread::runnable_group::run() {
   auto_ptr< iterator<runnable> > i = runnables.iter_create();
   while (i->has_next())
      i->next().run();
}

/***************************************************************************
 * Child runnable implementation.
 ***************************************************************************/

/*
 * Constructor.
 */
child_thread::child_runnable::child_runnable(
   runnable&     runnable_obj,
   unsigned long n_processors)
 : _runnable_obj(runnable_obj),
   _n_processors(n_processors),
   _exception(NULL)
{ }

/*
 * Destructor.
 */
child_thread::child_runnable::~child_runnable() {
   /* do nothing */
}

/*
 * Run method.
 * Wrap the run method of the thread.
 */
void child_thread::child_runnable::run() {
   /* set # of available processors for child thread */
   thread::processors(_n_processors);
   /* run thread, catch unhandled exceptions */
   try {
      _runnable_obj.run();
   } catch (throwable& e) {
      _exception.reset(e.clone());
   }
}

/*
 * Raise pending exception (if any).
 */
void child_thread::child_runnable::raise() {
   if (_exception.get() != NULL)
      _exception->raise();
}

} /* namespace threads */
} /* namespace concurrent */
