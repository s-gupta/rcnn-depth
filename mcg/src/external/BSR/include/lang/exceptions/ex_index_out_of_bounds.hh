/*
 * Index out of bounds exception.
 */
#ifndef LANG__EXCEPTIONS__EX_INDEX_OUT_OF_BOUNDS_HH
#define LANG__EXCEPTIONS__EX_INDEX_OUT_OF_BOUNDS_HH

#include "lang/exceptions/exception.hh"

namespace lang {
namespace exceptions {

/*
 * Index out of bounds exception.
 */
class ex_index_out_of_bounds : public exception {
public:
   /*
    * Constructor.
    */
   explicit ex_index_out_of_bounds(
      const char*,   /* message (use default if NULL) */
      unsigned long  /* index */
   );
   
   /*
    * Copy constructor.
    */
   ex_index_out_of_bounds(const ex_index_out_of_bounds&);

   /*
    * Destructor.
    */
   virtual ~ex_index_out_of_bounds();

   /*
    * Clone the exception.
    */
   virtual ex_index_out_of_bounds* clone() const;

   /*
    * Throw the exception.
    */
   virtual void raise() const;

   /*
    * Get the out-of-bounds index.
    */
   virtual unsigned long index() const;

protected:
   unsigned long _index;   /* index */
};

} /* namespace exceptions */
} /* namespace lang */

#endif
