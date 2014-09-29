/*
 * Auto collection.
 *
 * An auto collection owns both the collection to which it holds a pointer and
 * the items contained within that collection.  Upon its destruction, both the 
 * collection and its contents are deleted.  Copying an auto_collection copies 
 * the pointer to the collection and transfers ownership.  For any given 
 * collection, at most one auto_collection should own it.
 *
 * If one desires to automatically manage the collection, but not its
 * contents, consider using an auto_ptr instead.
 *
 * Note that collections of const elements may safely be owned by an 
 * auto_collection.  In this case, the const elements are not deleted 
 * during destruction.
 *
 * WARNING: A collection containing multiple references to the same element 
 * should not be owned by an auto_collection.  Otherwise, multiple deletes of 
 * the same object will be attempted at destruction time.
 */
#ifndef COLLECTIONS__POINTERS__AUTO_COLLECTION_HH
#define COLLECTIONS__POINTERS__AUTO_COLLECTION_HH

#include "collections/abstract/collection.hh"
#include "config/safety.hh"
#include "functors/iterable_functors.hh"
#include "lang/exceptions/ex_null_pointer_dereference.hh"
#include "lang/null.hh"

/*
 * Enable/disable null pointer dereference checking for auto collections.
 */
#ifdef CONFIG__SAFETY__CHECK_DEREFERENCE
   #define COLLECTIONS__POINTERS__AUTO_COLLECTION__CHECK_DEREFERENCE (true)
#else
   #define COLLECTIONS__POINTERS__AUTO_COLLECTION__CHECK_DEREFERENCE (false)
#endif

namespace collections {
namespace pointers {
/*
 * Imports.
 */
using collections::abstract::collection;
using functors::iter_functors;
using lang::exceptions::ex_null_pointer_dereference;

/*
 * Declare auto_collection class.
 */
template <typename T, typename C>
class auto_collections;

/*
 * Auto collection reference.
 *
 * This wrapper class provides auto_collection with reference semantics.
 */
template <typename C>
class auto_collection_ref {
public:
   /*
    * Friend classes.
    */
   template <typename U, typename D> friend class auto_collection;

protected:
   /*
    * Constructor.
    */
   explicit inline auto_collection_ref(C*);

   C* _p;   /* pointer to collection */
};

/*
 * Auto collection reference constructor.
 */
template <typename C>
inline auto_collection_ref<C>::auto_collection_ref(C* p) : _p(p) { }

/*
 * Auto collection.
 * T is the element type and C is the collection type.
 */
template < typename T, typename C = collection<T> >
class auto_collection {
public:
   /*
    * Constructor.
    */
   explicit inline auto_collection(C* = NULL);

   /*
    * Copy constructors.
    *
    * The constructed auto_collection assumes ownership 
    * and the original auto_collection is set to NULL.
    *
    * Note that a const auto_collection cannot be copied, 
    * since it is not allowed to transfer ownership.
    */
   auto_collection(auto_collection<T,C>&);
   template <typename D> auto_collection(auto_collection<T,D>&);

   /*
    * Assignment operators.
    *
    * The assigned auto_collection assumes ownership 
    * and the original auto_collection is set to NULL.
    *
    * Any collection previously owned by the assigned 
    * auto_collection is deleted.
    *
    * Note that a const auto_collection cannot appear on 
    * the right hand side of an assignment, since it is 
    * not allowed to transfer ownership.
    */
   auto_collection<T,C>& operator=(auto_collection<T,C>&);
   template <typename D> auto_collection<T,C>& operator=(auto_collection<T,D>&);

   /*
    * Destructor.
    * Delete the owned collection (if any) as well as the 
    * contents of the collection (if they are non-const items).
    */
   ~auto_collection();

   /*
    * Smart pointer dereference.
    * Return owned collection by reference.
    */
   inline C& operator*() const;

   /*
    * Smart pointer dereference.
    * Return pointer to the owned collection.
    */
   inline C* operator->() const;

   /*
    * Get raw pointer to owned collection.
    * The auto_collection still owns the collection.
    */
   inline C* get() const;

   /*
    * Release the owned collection (set the auto_collection to NULL).
    * Return a pointer to the collection.
    */
   C* release();

   /*
    * Delete the owned collection and set the auto_collection
    * to own the given pointer (defaults to NULL) instead.
    */
   void reset(C* = NULL);
   
   /*
    * Automatic conversions to/from auto_collection_ref.
    *
    * These conversions are necessary in order to allow
    * functions to return auto collections by value (since 
    * copying a const auto_collection is forbidden).
    */
   auto_collection(auto_collection_ref<C>);
   auto_collection<T,C>& operator=(auto_collection_ref<C>);
   template <typename D> operator auto_collection_ref<D>();

   /*
    * Automatic conversion to auto_collection of different type.
    */
   template <typename D> operator auto_collection<T,D>();

protected:
   C* _p;   /* pointer to owned collection */
};

/*
 * Constructor.
 */
template <typename T, typename C>
inline auto_collection<T,C>::auto_collection(C* p) : _p(p) { }

/*
 * Copy constructors.
 */
template <typename T, typename C>
auto_collection<T,C>::auto_collection(auto_collection<T,C>& a)
 : _p(a.release())
{ }

template <typename T, typename C>
template <typename D> 
auto_collection<T,C>::auto_collection(auto_collection<T,D>& a)
 : _p(a.release())
{ }

/*
 * Assignment operators.
 */
template <typename T, typename C>
auto_collection<T,C>& auto_collection<T,C>::operator=(auto_collection<T,C>& a) {
   this->reset(a.release());
   return *this;
}

template <typename T, typename C>
template <typename D> 
auto_collection<T,C>& auto_collection<T,C>::operator=(auto_collection<T,D>& a) {
   this->reset(a.release());
   return *this;
}

/*
 * Destructor.
 * Delete the owned collection (if any) as well as the 
 * contents of the collection (if they are non-const items).
 */
template <typename T, typename C>
auto_collection<T,C>::~auto_collection() {
   if (_p != NULL)
      _p->iter(iter_functors<T>::f_delete());
   delete _p;
}

/*
 * Smart pointer dereference.
 * Return owned collection by reference.
 */
template <typename T, typename C>
inline C& auto_collection<T,C>::operator*() const {
   #ifdef COLLECTIONS__POINTERS__AUTO_COLLECTION__CHECK_DEREFERENCE
      if (_p == NULL)
         throw ex_null_pointer_dereference();
   #endif
   return *_p;
}

/*
 * Smart pointer dereference.
 * Return pointer to the owned collection.
 */
template <typename T, typename C>
inline C* auto_collection<T,C>::operator->() const {
   #ifdef COLLECTIONS__POINTERS__AUTO_COLLECTION__CHECK_DEREFERENCE
      if (_p == NULL)
         throw ex_null_pointer_dereference();
   #endif
   return _p;
}

/*
 * Get raw pointer to owned collection.
 * The auto_collection still owns the collection.
 */
template <typename T, typename C>
inline C* auto_collection<T,C>::get() const {
   return _p;
}

/*
 * Release the owned collection (set the auto_collection to NULL).
 * Return a pointer to the collection.
 */
template <typename T, typename C>
C* auto_collection<T,C>::release() {
   C* p = _p;
   _p = NULL;
   return p;
}

/*
 * Delete the owned collection and set the auto_collection
 * to own the given pointer (defaults to NULL) instead.
 */
template <typename T, typename C>
void auto_collection<T,C>::reset(C* p) {
   if (_p != p) {
      if (_p != NULL)
         _p->iter(iter_functors<T>::f_delete());
      delete _p;
      _p = p;
   }
}

/*
 * Automatic conversions to/from auto_collection_ref.
 */
template <typename T, typename C>
auto_collection<T,C>::auto_collection(auto_collection_ref<C> r) : _p(r._p) { }

template <typename T, typename C>
auto_collection<T,C>& auto_collection<T,C>::operator=(auto_collection_ref<C> r) {
   if (r._p != _p) {
      if (_p != NULL)
         _p->iter(iter_functors<T>::f_delete());
      delete _p;
      _p = r._p;
   }
   return *this;
}
      
template <typename T, typename C>
template <typename D> 
auto_collection<T,C>::operator auto_collection_ref<D>() {
   return auto_collection_ref<D>(this->release());
}

/*
 * Automatic conversion to auto_collection of different type.
 */
template <typename T, typename C>
template <typename D> 
auto_collection<T,C>::operator auto_collection<T,D>() {
   return auto_collection<T,D>(this->release());
}

} /* namespace pointers */
} /* namespace collections */

#endif
