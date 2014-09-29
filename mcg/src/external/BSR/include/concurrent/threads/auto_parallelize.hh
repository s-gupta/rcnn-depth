/*
 * Auto parallelize.
 *
 * An auto_parallelize object sets the number of processors available to the
 * current thread (as indicated by thread::processors()) upon creation and
 * restores the previous number of available processors upon deletion.
 *
 * Programs use the value of thread::processors() as an estimate of the amount
 * of thread-level parallelism they may attempt to utilize.
 *
 * An auto_parallelize object provides a convenient way of controlling the 
 * number of processors used within a block of code.
 */
#ifndef CONCURRENT__THREADS__AUTO_PARALLELIZE_HH
#define CONCURRENT__THREADS__AUTO_PARALLELIZE_HH

namespace concurrent {
namespace threads {

/*
 * Auto parallelize.
 */
class auto_parallelize {
public:
   /*
    * Constructor.
    * Set the number of processors available to the current thread.
    * An argument of zero resets the processor count to the actual
    * number of processors on the machine.
    */
   explicit auto_parallelize(unsigned long = 0);

   /*
    * Constructor.
    * Multiply the number of processors available to the current thread by 
    * the specified factor.  An argument of zero resets the processor count
    * to the actual number of processors on the machine.
    */
   explicit auto_parallelize(double);
  
   /*
    * Destructor.
    * Reset the number of processors available to the current thread
    * to its value before creation of the auto_parallelize object.
    */
   ~auto_parallelize();

protected:
   unsigned long _n_processors;  /* # of processors to restore upon deletion */

private:
   /*
    * Private copy constructor.
    * Copying auto_parallelize objects is not permitted.
    */
   explicit auto_parallelize(const auto_parallelize&);
};

} /* namespace threads */
} /* namespace concurrent */

#endif
