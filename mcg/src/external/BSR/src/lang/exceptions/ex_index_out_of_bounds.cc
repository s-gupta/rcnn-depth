/*
 * Index out of bounds exception.
 */
#include "lang/exceptions/exception.hh"
#include "lang/exceptions/ex_index_out_of_bounds.hh"
#include "lang/null.hh"

namespace lang {
namespace exceptions {

/*
 * Constructor.
 */
ex_index_out_of_bounds::ex_index_out_of_bounds(
   const char* msg,
   unsigned long idx)
 : exception(
      ((msg != NULL) ? msg : "index out of bounds")
   ),
   _index(idx)
{ }

/*
 * Copy constructor.
 */
ex_index_out_of_bounds::ex_index_out_of_bounds(const ex_index_out_of_bounds& e)
 : exception(e),
   _index(e._index) 
{ }

/*
 * Destructor.
 */
ex_index_out_of_bounds::~ex_index_out_of_bounds() {
   /* do nothing */
}

/*
 * Clone the exception.
 */
ex_index_out_of_bounds* ex_index_out_of_bounds::clone() const {
   return (new ex_index_out_of_bounds(*this));
}

/*
 * Throw the exception.
 */
void ex_index_out_of_bounds::raise() const {
   throw *this;
}

/*
 * Get the out-of-bounds index.
 */
unsigned long ex_index_out_of_bounds::index() const {
   return _index;
}

} /* namespace exceptions */
} /* namespace lang */
