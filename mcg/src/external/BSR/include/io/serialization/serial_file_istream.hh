/*
 * Serial file input stream.
 */
#ifndef IO__SERIALIZATION__SERIAL_FILE_ISTREAM_HH
#define IO__SERIALIZATION__SERIAL_FILE_ISTREAM_HH

#include "io/serialization/serial_istream.hh"
#include "io/streams/ifstream.hh"
#include "lang/string.hh"

namespace io {
namespace serialization {
/*
 * Imports.
 */
using io::streams::ifstream;
using lang::string;

/*
 * Base class for serial_file_istream.
 */
class serial_file_istream_base {
public:
   /*
    * Destructor.
    * Close the file.
    */
   virtual ~serial_file_istream_base();

protected:
   /*
    * Constructor.
    * Open the specified file for input.
    */
   explicit serial_file_istream_base(const string<>& /* filename */);
   
   ifstream _f;   /* file stream */
};

/*
 * Serial file input stream.
 */
class serial_file_istream : protected serial_file_istream_base,
                            public    serial_istream {
public:
   /*
    * Constructor.
    * Open the specified file for use in deserialization.
    */
   explicit serial_file_istream(const string<>& /* filename */);

   /* 
    * Destructor.
    * Close the file.
    */
   virtual ~serial_file_istream();
   
   /*
    * Deserialization operators.
    * Read the specified data from the file stream.
    */
   virtual serial_file_istream& operator>>(bool&);
   virtual serial_file_istream& operator>>(char&);
   virtual serial_file_istream& operator>>(unsigned char&);
   virtual serial_file_istream& operator>>(short&);
   virtual serial_file_istream& operator>>(unsigned short&);
   virtual serial_file_istream& operator>>(int&);
   virtual serial_file_istream& operator>>(unsigned int&);
   virtual serial_file_istream& operator>>(long&);
   virtual serial_file_istream& operator>>(unsigned long&);
   virtual serial_file_istream& operator>>(long long&);
   virtual serial_file_istream& operator>>(unsigned long long&);
   virtual serial_file_istream& operator>>(float&);
   virtual serial_file_istream& operator>>(double&);
   virtual serial_file_istream& operator>>(long double&);
};

} /* namespace serialization */
} /* namespace io */

#endif
