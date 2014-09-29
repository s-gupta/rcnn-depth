/*
 * Runnable interface.
 *
 * This interface should be extended by any object to be executed by a thread.
 */
#ifndef CONCURRENT__THREADS__RUNNABLE_HH
#define CONCURRENT__THREADS__RUNNABLE_HH

namespace concurrent {
namespace threads {

/*
 * Runnable interface.
 */
class runnable {
public:
   /*
    * Destructor.
    */
   virtual ~runnable() = 0;
   
   /*
    * Method to be executed by the thread.
    */
   virtual void run() = 0;
};

} /* namespace threads */
} /* namespace concurrent */

#endif
