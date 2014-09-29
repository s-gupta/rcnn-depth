/*
 * Output stream for serialization.
 */
#ifndef IO__SERIALIZATION__SERIAL_OUTPUT_STREAM_HH
#define IO__SERIALIZATION__SERIAL_OUTPUT_STREAM_HH

namespace io {
namespace serialization {

/*
 * Output stream for serialization.
 */
class serial_output_stream {
public:
   /*
    * Destructor.
    */
   virtual ~serial_output_stream() = 0;

   /*
    * Serialization operators.
    * Write the specified data to the stream.
    */
   virtual serial_output_stream& operator<<(const bool&) = 0;
   virtual serial_output_stream& operator<<(const char&) = 0;
   virtual serial_output_stream& operator<<(const unsigned char&) = 0;
   virtual serial_output_stream& operator<<(const short&) = 0;
   virtual serial_output_stream& operator<<(const unsigned short&) = 0;
   virtual serial_output_stream& operator<<(const int&) = 0;
   virtual serial_output_stream& operator<<(const unsigned int&) = 0;
   virtual serial_output_stream& operator<<(const long&) = 0;
   virtual serial_output_stream& operator<<(const unsigned long&) = 0;
   virtual serial_output_stream& operator<<(const long long&) = 0;
   virtual serial_output_stream& operator<<(const unsigned long long&) = 0;
   virtual serial_output_stream& operator<<(const float&) = 0;
   virtual serial_output_stream& operator<<(const double&) = 0;
   virtual serial_output_stream& operator<<(const long double&) = 0;
};

} /* namespace serialization */
} /* namespace io */

#endif
