/*
 * Type function map.
 *
 * A type function map associates type information with function pointers.
 */
#ifndef LANG__TYPES__TYPE_FUNCTION_MAP_HH
#define LANG__TYPES__TYPE_FUNCTION_MAP_HH

#include "lang/pointers/auto_ptr.hh"
#include "lang/types/type_identifier.hh"
#include "lang/types/type_info.hh"

namespace lang {
namespace types {
/*
 * Imports.
 */
using lang::pointers::auto_ptr;

/*
 * Declare private implementation class.
 */
class type_function_map_impl_private;

/*
 * Concrete version of a type function map using void pointers.
 * This concrete version is made type-safe by the template wrapper below.
 */
class type_function_map_impl {
public:
   /*
    * Friend classes.
    */
   template <typename F> friend class type_function_map;

   /*
    * Destructor.
    */
   virtual ~type_function_map_impl();

protected:
   /*
    * Constructor.
    */
   type_function_map_impl();

   /*
    * Copy constructor.
    */
   explicit type_function_map_impl(const type_function_map_impl&);

   /*
    * Add a type -> function mapping (given the type info).
    */
   void add(const type_info&, void*);
   
   /*
    * Add a type -> function mapping (given the type identifier).
    */
   void add(const type_identifier&, void*);

   /*
    * Lookup the function for the given type (using the type info).
    */
   void* lookup(const type_info&) const;
   
   /*
    * Lookup the function for the given type (using the type identifier).
    */
   void* lookup(const type_identifier&) const;

private:
   auto_ptr<type_function_map_impl_private> _impl_private;
};

/*
 * Type function map.
 */
template <typename F>
class type_function_map { };

template <typename F>
class type_function_map<F*> {
public:
   /*
    * Constructor.
    * Return an empty map.
    */
   type_function_map();

   /*
    * Copy constructor.
    */
   explicit type_function_map(const type_function_map<F>&);

   /*
    * Destructor.
    */
   virtual ~type_function_map();

   /*
    * Add a type -> function mapping (given the type info).
    */
   void add(const type_info&, F*);
   
   /*
    * Add a type -> function mapping (given the type identifier).
    */
   void add(const type_identifier&, F*);

   /*
    * Lookup the function for the given type (using the type info).
    */
   F* lookup(const type_info&) const;
   
   /*
    * Lookup the function for the given type (using the type identifier).
    */
   F* lookup(const type_identifier&) const;

protected:
   type_function_map_impl _impl;
};

/*
 * Constructor.
 * Return an empty map.
 */
template <typename F>
type_function_map<F*>::type_function_map()
 : _impl()
{ }

/*
 * Copy constructor.
 */
template <typename F>
type_function_map<F*>::type_function_map(const type_function_map<F>& m)
 : _impl(m._impl)
{ }

/*
 * Destructor.
 */
template <typename F>
type_function_map<F*>::~type_function_map() {
   /* do nothing */
}

/*
 * Add a type -> function mapping (given the type info).
 */
template <typename F>
void type_function_map<F*>::add(const type_info& t_info, F* f) {
   _impl.add(
      t_info, reinterpret_cast<void*>(reinterpret_cast<unsigned long>(f))
   );
}

/*
 * Add a type -> function mapping (given the type identifier).
 */
template <typename F>
void type_function_map<F*>::add(const type_identifier& t_id, F* f) {
   _impl.add(
      t_id, reinterpret_cast<void*>(reinterpret_cast<unsigned long>(f))
   );
}

/*
 * Lookup the function for the given type (using the type info).
 */
template <typename F>
F* type_function_map<F*>::lookup(const type_info& t_info) const {
   void* f = _impl.lookup(t_info);
   return reinterpret_cast<F*>(reinterpret_cast<unsigned long>(f));
}

/*
 * Lookup the function for the given type (using the type identifier).
 */
template <typename F>
F* type_function_map<F*>::lookup(const type_identifier& t_id) const {
   void* f = _impl.lookup(t_id);
   return reinterpret_cast<F*>(reinterpret_cast<unsigned long>(f));
}

/*
 * Type function map initialization.
 */
template <typename F>
class type_function_map_init { };

template <typename F>
class type_function_map_init<F*> {
public:
   /*
    * Constructor.
    * Add the specified type -> function mapping to the map.
    */
   explicit type_function_map_init(
      type_function_map<F*>& m, const type_info& t_info, F* f)
   {
      m.add(t_info, f);
   }
   
   /*
    * Constructor.
    * Add the specified type -> function mapping to the map.
    */
   explicit type_function_map_init(
      type_function_map<F*>& m, const type_identifier& t_id, F* f)
   {
      m.add(t_id, f);
   }
};

} /* namespace types */
} /* namespace lang */

#endif
