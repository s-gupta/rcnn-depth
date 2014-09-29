/*
 * Auto parallelize.
 */
#include "concurrent/threads/auto_parallelize.hh"
#include "concurrent/threads/thread.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "math/math.hh"

namespace concurrent {
namespace threads {
/*
 * Imports.
 */
using lang::exceptions::ex_invalid_argument;

/*
 * Constructor.
 * Set the number of processors available to the current thread.
 * An argument of zero resets the processor count to the actual
 * number of processors on the machine.
 */
auto_parallelize::auto_parallelize(unsigned long n)
 : _n_processors(thread::processors())
{
   thread::processors(n);
}

/*
 * Constructor.
 * Multiply the number of processors available to the current thread by 
 * the specified factor.  An argument of zero resets the processor count
 * to the actual number of processors on the machine.
 */
auto_parallelize::auto_parallelize(double factor)
 : _n_processors(thread::processors())
{
   /* check argument */
   if (factor < 0)
      throw ex_invalid_argument("factor must be >= 0");
   /* set number of processors */
   unsigned long n = static_cast<unsigned long>(
      math::ceil(factor * static_cast<double>(_n_processors))
   );
   thread::processors(n);
}
   
/*
 * Private copy constructor.
 * Copying auto_parallelize objects is not permitted.
 */
auto_parallelize::auto_parallelize(const auto_parallelize& p)
 : _n_processors(p._n_processors)
{ }

/*
 * Destructor.
 * Reset the number of processors available to the current thread
 * to its value before creation of the auto_parallelize object.
 */
auto_parallelize::~auto_parallelize() {
   thread::processors(_n_processors);
}

} /* namespace threads */
} /* namespace concurrent */
