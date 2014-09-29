/*
 * Run-time type information.
 * Type info can be obtained using the typeid(...) operator.
 */
#ifndef LANG__TYPES__TYPE_INFO_HH
#define LANG__TYPES__TYPE_INFO_HH

#include <typeinfo>

namespace lang {
namespace types {
/*
 * Export std::type_info class.
 */
using std::type_info;

} /* namespace types */
} /* namespace lang */

#endif
