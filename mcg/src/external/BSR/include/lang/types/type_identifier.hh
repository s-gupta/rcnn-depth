/*
 * Type identifier.
 *
 * A type identifier is a tag that represents a type.
 */
#ifndef LANG__TYPES__TYPE_IDENTIFIER_HH
#define LANG__TYPES__TYPE_IDENTIFIER_HH

#include "interfaces/comparable.hh"
#include "io/serialization/serial_input_stream.hh"
#include "io/serialization/serial_output_stream.hh"
#include "lang/pointers/auto_ptr.hh"
#include "lang/types/type_info.hh"

namespace lang {
namespace types {
/*
 * Imports.
 */
using interfaces::comparable;
using io::serialization::serial_input_stream;
using io::serialization::serial_output_stream;
using lang::pointers::auto_ptr;

/*
 * Declare type identifier implementation class.
 */
class type_identifier_impl;

/*
 * Type identifier.
 */
class type_identifier : public comparable<type_identifier> {
public:
   /*
    * Constructor.
    * Return the type identifier corresponding to the given type information.
    */
   explicit type_identifier(const type_info&);

   /*
    * Copy constructor.
    */
   explicit type_identifier(const type_identifier&);

   /*
    * Destructor.
    */
   ~type_identifier();

   /*
    * Serialize.
    */
   void serialize(serial_output_stream&) const;

   /*
    * Deserialize.
    */
   static auto_ptr<type_identifier> deserialize(serial_input_stream&);

   /*
    * Return the unique name corresponding to the type id.
    */
   const char* name() const;

   /*
    * Comparable comparison method.
    */
   int compare_to(const type_identifier&) const;

protected:
   /*
    * Protected constructor.
    */
   type_identifier();

private:
   auto_ptr<type_identifier_impl> _impl;
};

} /* namespace types */
} /* namespace lang */

#endif
