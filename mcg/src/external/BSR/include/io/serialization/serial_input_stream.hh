/*
 * Input stream for deserialization.
 */
#ifndef IO__SERIALIZATION__SERIAL_INPUT_STREAM_HH
#define IO__SERIALIZATION__SERIAL_INPUT_STREAM_HH

namespace io {
namespace serialization {

/*
 * Input stream for deserialization.
 */
class serial_input_stream {
public:
   /*
    * Destructor.
    */
   virtual ~serial_input_stream() = 0;

   /*
    * Deserialization operators.
    * Read the specified data from the stream.
    */
   virtual serial_input_stream& operator>>(bool&) = 0;
   virtual serial_input_stream& operator>>(char&) = 0;
   virtual serial_input_stream& operator>>(unsigned char&) = 0;
   virtual serial_input_stream& operator>>(short&) = 0;
   virtual serial_input_stream& operator>>(unsigned short&) = 0;
   virtual serial_input_stream& operator>>(int&) = 0;
   virtual serial_input_stream& operator>>(unsigned int&) = 0;
   virtual serial_input_stream& operator>>(long&) = 0;
   virtual serial_input_stream& operator>>(unsigned long&) = 0;
   virtual serial_input_stream& operator>>(long long&) = 0;
   virtual serial_input_stream& operator>>(unsigned long long&) = 0;
   virtual serial_input_stream& operator>>(float&) = 0;
   virtual serial_input_stream& operator>>(double&) = 0;
   virtual serial_input_stream& operator>>(long double&) = 0;
};

} /* namespace serialization */
} /* namespace io */

#endif
