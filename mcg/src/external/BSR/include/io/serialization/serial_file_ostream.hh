/*
 * Serial file output stream.
 */
#ifndef IO__SERIALIZATION__SERIAL_FILE_OSTREAM_HH
#define IO__SERIALIZATION__SERIAL_FILE_OSTREAM_HH

#include "io/serialization/serial_ostream.hh"
#include "io/streams/ofstream.hh"
#include "lang/string.hh"

namespace io {
namespace serialization {
/*
 * Imports.
 */
using io::streams::ofstream;
using lang::string;

/*
 * Base class for serial_file_ostream.
 */
class serial_file_ostream_base {
public:
   /*
    * Destructor.
    * Close the file.
    */
   virtual ~serial_file_ostream_base();

protected:
   /*
    * Constructor.
    * Open the specified file for output.
    */
   explicit serial_file_ostream_base(const string<>& /* filename */);

   ofstream _f;   /* file stream */
};

/*
 * Serial file output stream.
 */
class serial_file_ostream : protected serial_file_ostream_base,
                            public    serial_ostream {
public:
   /*
    * Constructor.
    * Open the specified file for use in serialization.
    */
   explicit serial_file_ostream(const string<>& /* filename */);

   /* 
    * Destructor.
    * Close the file.
    */
   virtual ~serial_file_ostream();

   /*
    * Serialization operators.
    * Write the specified data to the file stream.
    */
   virtual serial_file_ostream& operator<<(const bool&);
   virtual serial_file_ostream& operator<<(const char&);
   virtual serial_file_ostream& operator<<(const unsigned char&);
   virtual serial_file_ostream& operator<<(const short&);
   virtual serial_file_ostream& operator<<(const unsigned short&);
   virtual serial_file_ostream& operator<<(const int&);
   virtual serial_file_ostream& operator<<(const unsigned int&);
   virtual serial_file_ostream& operator<<(const long&);
   virtual serial_file_ostream& operator<<(const unsigned long&);
   virtual serial_file_ostream& operator<<(const long long&);
   virtual serial_file_ostream& operator<<(const unsigned long long&);
   virtual serial_file_ostream& operator<<(const float&);
   virtual serial_file_ostream& operator<<(const double&);
   virtual serial_file_ostream& operator<<(const long double&);
};

} /* namespace serialization */
} /* namespace io */

#endif
