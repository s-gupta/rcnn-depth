/*
 * Standard stream manipulators.
 */
#ifndef IO__STREAMS__IOMANIP_HH
#define IO__STREAMS__IOMANIP_HH

#include <iomanip>

namespace io {
namespace streams {
namespace iomanip {
/*
 * Export standard manipulators.
 */
using std::setiosflags;
using std::resetiosflags;
using std::setw;
using std::setfill;
using std::setprecision;
using std::setbase;

} /* namespace iomanip */
} /* namespace streams */
} /* namespace io */

#endif
