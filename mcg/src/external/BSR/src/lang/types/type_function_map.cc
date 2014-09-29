/*
 * Type function map.
 */
#include "collections/list.hh"
#include "collections/map.hh"
#include "collections/pointers/auto_collection.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "lang/pointers/auto_ptr.hh"
#include "lang/types/type_function_map.hh"
#include "lang/types/type_identifier.hh"
#include "lang/types/type_info.hh"

namespace lang {
namespace types {
/*
 * Imports.
 */
using collections::list;
using collections::map;
using collections::pointers::auto_collection;
using lang::exceptions::ex_invalid_argument;
using lang::pointers::auto_ptr;

/*
 * Define private data used by type function map implementation.
 */
class type_function_map_impl_private {
public:
   auto_collection< type_identifier, list<type_identifier> > _domain;
   auto_collection< void*, list<void*> >                     _range;
   map< type_identifier, void* >                             _map;
};

/*
 * Constructor.
 */
type_function_map_impl::type_function_map_impl()
 : _impl_private(new type_function_map_impl_private())
{
   _impl_private->_domain.reset(new list<type_identifier>());
   _impl_private->_range.reset(new list<void*>());
}

/*
 * Copy constructor.
 */
type_function_map_impl::type_function_map_impl(const type_function_map_impl& m)
 :  _impl_private(new type_function_map_impl_private())
{
   /* initialize domain and range */
   _impl_private->_domain.reset(new list<type_identifier>());
   _impl_private->_range.reset(new list<void*>());
   /* rebuild map */
   list<type_identifier>::iterator_t iter_domain(*(_impl_private->_domain));
   list<void*>::iterator_t           iter_range(*(_impl_private->_range));
   while (iter_domain.has_next())
      this->add(iter_domain.next(), iter_range.next());
}

/*
 * Destructor.
 */
type_function_map_impl::~type_function_map_impl() {
   /* do nothing */
}

/*
 * Add a type -> function mapping (given the type info).
 */
void type_function_map_impl::add(const type_info& t_info, void* f) {
   const type_identifier t_id(t_info);
   this->add(t_id, f);
}

/*
 * Add a type -> function mapping (given the type identifier).
 */
void type_function_map_impl::add(const type_identifier& t_id, void* f) {
   /* copy elements to be mapped */
   auto_ptr<type_identifier> t_id_copy(new type_identifier(t_id));
   auto_ptr<void*>           f_copy(new (void*)(f));
   type_identifier& t_id_copy_ref(*t_id_copy);
   void*&           f_copy_ref(*f_copy);
   /* check that mapping does not already exist */
   if (_impl_private->_map.contains(t_id_copy_ref))
      throw ex_invalid_argument(
         "cannot overwrite existing mapping in type function map"
      );
   /* add mapping */
   _impl_private->_domain->add(t_id_copy_ref); t_id_copy.release();
   _impl_private->_range->add(f_copy_ref);     f_copy.release();
   _impl_private->_map.add(t_id_copy_ref, f_copy_ref);
}

/*
 * Lookup the function for the given type (using the type info).
 */
void* type_function_map_impl::lookup(const type_info& t_info) const {
   const type_identifier t_id(t_info);
   return this->lookup(t_id);
}

/*
 * Lookup the function for the given type (using the type identifier).
 */
void* type_function_map_impl::lookup(const type_identifier& t_id) const {
   return _impl_private->_map.find_image(t_id);
}

} /* namespace types */
} /* namespace lang */
