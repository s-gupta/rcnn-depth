/*
 * Serial output stream built on top of a standard ostream.
 */
#ifndef IO__SERIALIZATION__SERIAL_OSTREAM_HH
#define IO__SERIALIZATION__SERIAL_OSTREAM_HH

#include "io/serialization/serial_output_stream.hh"
#include "io/streams/ostream.hh"

namespace io {
namespace serialization {
/*
 * Imports.
 */
using io::streams::ostream;

/*
 * Declare serial_istream class.
 */
class serial_istream;

/*
 * Serial output stream built on top of a standard ostream.
 */
class serial_ostream : public serial_output_stream {
public:
   /*
    * Friend class.
    */
   friend class serial_istream;

   /*
    * Constructor.
    * Attach the ostream for use as output.
    * The ostream must be ready for writing.
    */
   explicit serial_ostream(ostream&);

   /*
    * Copy constructor.
    */
   explicit serial_ostream(serial_ostream&);

   /*
    * Destructor.
    */
   virtual ~serial_ostream();

   /*
    * Serialization operators.
    * Write the specified data to the stream.
    */
   virtual serial_ostream& operator<<(const bool&);
   virtual serial_ostream& operator<<(const char&);
   virtual serial_ostream& operator<<(const unsigned char&);
   virtual serial_ostream& operator<<(const short&);
   virtual serial_ostream& operator<<(const unsigned short&);
   virtual serial_ostream& operator<<(const int&);
   virtual serial_ostream& operator<<(const unsigned int&);
   virtual serial_ostream& operator<<(const long&);
   virtual serial_ostream& operator<<(const unsigned long&);
   virtual serial_ostream& operator<<(const long long&);
   virtual serial_ostream& operator<<(const unsigned long long&);
   virtual serial_ostream& operator<<(const float&);
   virtual serial_ostream& operator<<(const double&);
   virtual serial_ostream& operator<<(const long double&);

protected:
   /*
    * Return a value that changes under byte-reversal.
    * This returns the same value as serial_istream::endian_signature().
    */
   static unsigned short endian_signature();

   ostream& _os;  /* attached ostream */
};

} /* namespace serialization */
} /* namespace io */

#endif
