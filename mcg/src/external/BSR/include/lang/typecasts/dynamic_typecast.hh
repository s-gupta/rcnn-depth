/*
 * Dynamic typecast.
 *
 * This is equivalent to the built-in dynamic_cast, except that it throws an
 * ex_bad_cast exception on error instead of the built-in bad_cast exception.
 */
#ifndef LANG__TYPECASTS__DYNAMIC_TYPECAST_HH
#define LANG__TYPECASTS__DYNAMIC_TYPECAST_HH

#include "lang/exceptions/ex_bad_cast.hh"
#include <typeinfo>

namespace lang {
namespace typecasts {
/*
 * Imports.
 */
using lang::exceptions::ex_bad_cast;

/*
 * Dynamic typecast.
 */
template <typename T, typename U>
T dynamic_typecast(U& u) {
   try {
      return dynamic_cast<T>(u);
   } catch (std::bad_cast&) {
      throw ex_bad_cast(
         "dynamic typecast failed"
      );
   }
}

} /* namespace typecasts */
} /* namespace lang */

#endif
