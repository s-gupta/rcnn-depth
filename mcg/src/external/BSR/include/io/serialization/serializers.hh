/*
 * Serializers.
 */
#ifndef IO__SERIALIZATION__SERIALIZERS_HH
#define IO__SERIALIZATION__SERIALIZERS_HH

#include "io/serialization/serial_input_stream.hh"
#include "io/serialization/serial_output_stream.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"

namespace io {
namespace serialization {
/*
 * Imports.
 */
using lang::pointers::auto_ptr;

/*
 * Abstract serialization functor.
 */
template <typename T>
class serializer {
public:
   virtual ~serializer() { }
   virtual void serialize(serial_output_stream&, const T&) const = 0;
   virtual auto_ptr<T> deserialize(serial_input_stream&) const = 0;
};

template <typename T>
class serializer<const T> : public serializer<T> { };

/*
 * Default functors for serialization.
 *
 * These functors are valid for any built-in type (via template specialization
 * below) as well as for any class T specifiying its own serialization methods
 * with the following type signatures:
 *
 *    [virtual] void serialize(serial_output_stream&, ... ) const;
 *
 *    static auto_ptr<T> deserialize(serial_input_stream&, ... );
 *
 * where the virtual keyword is optional and ... indicates optional additional
 * arguments for which default values are specified.
 */
template <typename T>
class default_serializer : public serializer<const T> {
public:
   void serialize(serial_output_stream& s, const T& t) const {
      t.serialize(s);
   }

   auto_ptr<T> deserialize(serial_input_stream& s) const {
      return T::deserialize(s);
   }
};

template <typename T>
class default_serializer<const T> : public default_serializer<T> { };

/*
 * Globally accessible set of default serializers.
 */
template <typename T>
class serializers {
public:
   static const default_serializer<T>& s_default();
};

template <typename T>
const default_serializer<T>& serializers<T>::s_default() {
   static const default_serializer<T>* s = new default_serializer<T>();
   return *s;
}

/*
 * Generic functor for serialization of built-in types.
 */
template <typename T>
class built_in_serializer : public serializer<const T> {
public:
   void serialize(serial_output_stream& s, const T& t) const {
      s << t;
   }

   auto_ptr<T> deserialize(serial_input_stream& s) const {
      auto_ptr<T> t_ptr(new T());
      s >> (*t_ptr);
      return t_ptr;
   }
};

template <typename T>
class built_in_serializer<const T> : public built_in_serializer<T> { };

/*
 * Default functor for serializing a pointer.
 */
template <typename T>
class default_serializer<T*> : public serializer<T* const> {
public:
   void serialize(serial_output_stream& s, T* const& t_ptr) const {
      bool is_null = (t_ptr == NULL);
      s << is_null;
      if (!is_null) {
         const default_serializer<T>& slzr = serializers<T>::s_default();
         slzr.serialize(s, *t_ptr);
      }
   }

   auto_ptr<T> deserialize(serial_input_stream& s) const {
      bool is_null = false;
      s >> is_null;
      if (is_null) {
         return auto_ptr<T>(NULL);
      } else {
         const default_serializer<T>& slzr = serializers<T>::s_default();
         return slzr.deserialize(s);
      }
   }
};

/*
 * Default functors for serialization of built-in types.
 */
template <>
class default_serializer<bool>
 : public built_in_serializer<const bool> { };

template <>
class default_serializer<char>
 : public built_in_serializer<const char> { };

template <>
class default_serializer<unsigned char>
 : public built_in_serializer<const unsigned char> { };

template <>
class default_serializer<short>
 : public built_in_serializer<const short> { };

template <>
class default_serializer<unsigned short>
 : public built_in_serializer<const unsigned short> { };

template <>
class default_serializer<int>
 : public built_in_serializer<const int> { };

template <>
class default_serializer<unsigned int>
 : public built_in_serializer<const unsigned int> { };

template <>
class default_serializer<long>
 : public built_in_serializer<const long> { };

template <>
class default_serializer<unsigned long>
 : public built_in_serializer<const unsigned long> { };

template <>
class default_serializer<long long>
 : public built_in_serializer<const long long> { };

template <>
class default_serializer<unsigned long long>
 : public built_in_serializer<const unsigned long long> { };

template <>
class default_serializer<float>
 : public built_in_serializer<const float> { };

template <>
class default_serializer<double>
 : public built_in_serializer<const double> { };

template <>
class default_serializer<long double>
 : public built_in_serializer<const long double> { };

/*
 * Container serialization functors.
 *
 * A container serialization functor takes a serialization functor for type T
 * and turns it into a serialization functor for a type C which is a container
 * of T.  These functors are valid for any container type specifiying its own
 * serialization methods with the following type signatures:
 *
 * [virtual] void serialize(serial_output_stream&, const serializer<T>&) const;
 *
 * static auto_ptr<C> deserialize(serial_input_stream&, const serializer<T>&);
 *
 * where the virtual keyword is optional and the methods may also take
 * additional default arguments.
 */
template <typename C, typename T>
class container_serializer : public serializer<const C> {
public:
   container_serializer(
      const serializer<T>& slzr = serializers<T>::s_default()
   ) : _slzr(slzr) { }
   
   void serialize(serial_output_stream& s, const C& c) const {
      c.serialize(s, _slzr);
   }
   
   auto_ptr<C> deserialize(serial_input_stream& s) const {
      return C::deserialize(s, _slzr);
   }

protected:
   const serializer<T>& _slzr;
};

template <typename C, typename T>
class container_serializer<const C, T>
 : public container_serializer<C,T> { };

template <typename C, typename T>
class container_serializer<C, const T>
 : public container_serializer<C,T> { };

template <typename C, typename T>
class container_serializer<const C, const T>
 : public container_serializer<C,T> { };

/*
 * Globally accessible set of default container serializers.
 */
template <typename C, typename T>
class container_serializers {
public:
   static const container_serializer<C,T>& s_default();
};

template <typename C, typename T>
const container_serializer<C,T>& container_serializers<C,T>::s_default() {
   static const container_serializer<C,T>* s = new container_serializer<C,T>();
   return *s;
}

} /* namespace serialization */
} /* namespace io */

#endif
