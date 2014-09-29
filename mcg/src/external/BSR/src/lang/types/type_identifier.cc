/*
 * Type identifier.
 */
#include "io/serialization/serial_input_stream.hh"
#include "io/serialization/serial_output_stream.hh"
#include "lang/string.hh"
#include "lang/pointers/auto_ptr.hh"
#include "lang/types/type_identifier.hh"
#include "lang/types/type_info.hh"

namespace lang {
namespace types {
/*
 * Imports.
 */
using io::serialization::serial_input_stream;
using io::serialization::serial_output_stream;
using lang::pointers::auto_ptr;

/*
 * Define private data used by the type identifier implementation.
 */
class type_identifier_impl {
public:
   explicit type_identifier_impl(const char*);
   explicit type_identifier_impl(const type_info&);
   explicit type_identifier_impl(const type_identifier_impl&);
   string<> name_str;   /* name uniquely identifying the type */
};

type_identifier_impl::type_identifier_impl(const char* str)
 : name_str(str)
{ }

type_identifier_impl::type_identifier_impl(const type_info& t_info)
 : name_str(t_info.name())
{ }

type_identifier_impl::type_identifier_impl(const type_identifier_impl& t_impl)
 : name_str(t_impl.name_str)
{ }

/*
 * Protected constructor.
 */
type_identifier::type_identifier()
 : _impl(NULL)
{ }

/*
 * Constructor.
 * Return the type identifier corresponding to the given type information.
 */
type_identifier::type_identifier(const type_info& t_info)
 : _impl(new type_identifier_impl(t_info))
{ }

/*
 * Copy constructor.
 */
type_identifier::type_identifier(const type_identifier& t_id)
 : _impl(new type_identifier_impl(*(t_id._impl)))
{ }

/*
 * Destructor.
 */
type_identifier::~type_identifier() {
   /* do nothing */
}

/*
 * Serialize.
 */
void type_identifier::serialize(serial_output_stream& s) const {
   unsigned long size = _impl->name_str.size();
   s << size;
   for (unsigned long n = 0; n < size; n++)
      s << _impl->name_str[n];
}

/*
 * Deserialize.
 */
auto_ptr<type_identifier> type_identifier::deserialize(serial_input_stream& s) {
   unsigned long size = 0;
   s >> size;
   string<> name_str(size);
   for (unsigned long n = 0; n < size; n++)
      s >> name_str[n];
   auto_ptr<type_identifier> t_id(new type_identifier());
   t_id->_impl.reset(new type_identifier_impl(name_str));
   return t_id;
}

/*
 * Return the unique name corresponding to the type id.
 */
const char* type_identifier::name() const {
   return _impl->name_str;
}

/*
 * Comparable comparison method.
 */
int type_identifier::compare_to(const type_identifier& t_id) const {
   return _impl->name_str.compare_to(t_id._impl->name_str);
}

} /* namespace types */
} /* namespace lang */
