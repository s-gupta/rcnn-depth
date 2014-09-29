/*
 * Exception.
 */
#include "io/streams/ostream.hh"
#include "lang/exceptions/exception.hh"
#include "lang/null.hh"

#include <cstdlib>
#include <cstring>

namespace lang {
namespace exceptions {
/*
 * Imports.
 */
using io::streams::ostream;

/*
 * Constructor.
 * Create an exception with the specified message.
 */
exception::exception(const char* msg) {
   /* copy message */
   unsigned long msg_size = ((msg == NULL) ? 0 : std::strlen(msg));
   _message = new char[msg_size+1]();
   for (unsigned long n = 0; n < msg_size; n++)
      _message[n] = msg[n];
   _message[msg_size] = '\0';
   /* store stack trace */
   const int strace_max_size = 512;
   void *strace[strace_max_size];
   //int strace_size = backtrace(strace, strace_max_size);
   //_stack_trace = backtrace_symbols(strace, strace_size);
   _stack_trace_size = static_cast<unsigned long>(strace_max_size);
}

/*
 * Copy constructor.
 */
exception::exception(const exception& e) {
   /* copy message */
   unsigned long msg_size = ((e._message == NULL) ? 0 : std::strlen(e._message));
   _message = new char[msg_size+1]();
   for (unsigned long n = 0; n < msg_size; n++)
      _message[n] = e._message[n];
   _message[msg_size] = '\0';
   /* copy stack trace */
   _stack_trace_size = e._stack_trace_size;
   _stack_trace = static_cast<char**>(malloc(_stack_trace_size*sizeof(char*)));
   for (unsigned long n = 0; n < _stack_trace_size; n++)
      _stack_trace[n] = e._stack_trace[n];
}

/*
 * Destructor.
 */
exception::~exception() {
   delete [] _message;
   free(_stack_trace);
}

/*
 * Clone the exception.
 */
exception* exception::clone() const {
   return (new exception(*this));
}

/*
 * Throw the exception.
 */
void exception::raise() const {
   throw *this;
}

/*
 * Get the exception's message.
 */
const char* exception::what() const {
   return _message;
}

/*
 * Get the stack trace at the time the exception was thrown.
 */
const char* const* exception::stack_trace() const {
   return _stack_trace;
}

/*
 * Output message and stack trace to stream.
 */
ostream& operator<<(ostream& os, const exception& e) {
   os << "exception: " << e.what();
   for (unsigned long n = 0; n < e._stack_trace_size; n++)
      os << "\n" << e._stack_trace[n];
   return os;
}

} /* namespace exceptions */
} /* namespace lang */
