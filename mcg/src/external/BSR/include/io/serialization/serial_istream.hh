/*
 * Serial input stream built on top of a standard istream.
 */
#ifndef IO__SERIALIZATION__SERIAL_ISTREAM_HH
#define IO__SERIALIZATION__SERIAL_ISTREAM_HH

#include "io/serialization/serial_input_stream.hh"
#include "io/streams/istream.hh"

namespace io {
namespace serialization {
/*
 * Imports.
 */
using io::streams::istream;

/*
 * Serial input stream built on top of a standard istream.
 */
class serial_istream : public serial_input_stream {
public:
   /*
    * Constructor.
    * Attach the istream for use as input.
    * The istream must be ready for reading.
    */
   explicit serial_istream(istream&);

   /*
    * Copy constructor.
    */
   explicit serial_istream(serial_istream&);

   /*
    * Destructor.
    */
   virtual ~serial_istream();

   /*
    * Deserialization operators.
    * Read the specified data from the file stream.
    */
   virtual serial_istream& operator>>(bool&);
   virtual serial_istream& operator>>(char&);
   virtual serial_istream& operator>>(unsigned char&);
   virtual serial_istream& operator>>(short&);
   virtual serial_istream& operator>>(unsigned short&);
   virtual serial_istream& operator>>(int&);
   virtual serial_istream& operator>>(unsigned int&);
   virtual serial_istream& operator>>(long&);
   virtual serial_istream& operator>>(unsigned long&);
   virtual serial_istream& operator>>(long long&);
   virtual serial_istream& operator>>(unsigned long long&);
   virtual serial_istream& operator>>(float&);
   virtual serial_istream& operator>>(double&);
   virtual serial_istream& operator>>(long double&);

protected:
   /*
    * Return a value that changes under byte-reversal.
    * This returns the same value as serial_ostream::endian_signature().
    */
   static unsigned short endian_signature();

   istream& _is;              /* attached istream */
   bool     _reverse_endian;  /* reverse endianness of input? */
};

} /* namespace serialization */
} /* namespace io */

#endif
