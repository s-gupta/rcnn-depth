/*
 * Counted collection.
 *
 * A counted collection is a counted pointer to a collection that also owns 
 * the items in that collection.  Upon destruction of the last reference to the 
 * collection, both it and the items it contains are deleted (unless ownership 
 * is explicitly released by the user).  Copying a counted_collection increases
 * the reference count of the collection.
 *
 * If one desires to automatically manage the collection, but not its 
 * contents, consider using a counted_ptr instead.
 *
 * Note that collections of const elements may safely be used with a 
 * counted_collection.  In this case, the const elements are not deleted 
 * during destruction.
 *
 * WARNING: A collection containing multiple references to the same element 
 * should not be owned by any counted_collection.  Otherwise, multiple deletes 
 * of the same object will be attempted at destruction time.
 */
#ifndef COLLECTIONS__POINTERS__COUNTED_COLLECTION_HH
#define COLLECTIONS__POINTERS__COUNTED_COLLECTION_HH

#include "collections/abstract/collection.hh"
#include "concurrent/threads/synchronization/counter.hh"
#include "config/safety.hh"
#include "functors/iterable_functors.hh"
#include "lang/exceptions/ex_null_pointer_dereference.hh"
#include "lang/null.hh"

/*
 * Enable/disable null pointer dereference checking for counted collections.
 */
#if CONFIG__SAFETY__CHECK_DEREFERENCE
   #define COLLECTIONS__POINTERS__COUNTED_COLLECTION__CHECK_DEREFERENCE (true)
#else
   #define COLLECTIONS__POINTERS__COUNTED_COLLECTION__CHECK_DEREFERENCE (false)
#endif

namespace collections {
namespace pointers {
/*
 * Imports.
 */
using collections::abstract::collection;
using concurrent::threads::synchronization::counter;
using functors::iter_functors;
using lang::exceptions::ex_null_pointer_dereference;

/*
 * Auto collection.
 * T is the element type and C is the collection type.
 */
template < typename T, typename C = collection<T> >
class counted_collection {
public:
   /*
    * Friend classes.
    */
   template <typename U, typename D> friend class counted_collection;

   /*
    * Constructor.
    */
   explicit counted_collection(C* = NULL);
   
   /*
    * Copy constructors.
    * Set the copy to point to the same collection as the original.
    */
   counted_collection(const counted_collection<T,C>&);
   template <typename D> counted_collection(const counted_collection<T,D>&);
   
   /*
    * Assignment operators.
    * Set the copy to point to the same collection as the original.
    */
   counted_collection<T,C>& operator=(const counted_collection<T,C>&);
   template <typename D> counted_collection<T,C>& operator=(
      const counted_collection<T,D>&
   );

   /*
    * Destructor.
    * Delete the collection and its contents (if they are non-const items)
    * in the case that no other counted collections refer to it.
    */
   ~counted_collection();

   /*
    * Smart pointer dereference.
    * Return collection by reference.
    */
   inline C& operator*() const;

   /*
    * Smart pointer dereference.
    * Return pointer to collection.
    */
   inline C* operator->() const;

   /*
    * Get raw pointer to object.
    */
   inline C* get() const;

   /*
    * Release the collection (set the counted_collection to NULL).
    * Return a pointer to the collection.
    */
   C* release();
    
   /*
    * Delete the collection if no other counted collections refer to it.
    * Set the counted_collection to refer to the new collection (defaults to
    * NULL).
    */
   void reset(C* = NULL);

   /*
    * Automatic conversion to counted_collection of different type.
    */
   template <typename D> operator counted_collection<T,D>() const;

protected:
   C*       _p;      /* pointer to collection */
   counter* _count;  /* reference count for collection */
};

/*
 * Constructor.
 */
template <typename T, typename C>
counted_collection<T,C>::counted_collection(C* p)
 : _p(p), 
   _count(new counter(1))
{ }

/*
 * Copy constructors.
 */
template <typename T, typename C>
counted_collection<T,C>::counted_collection(const counted_collection<T,C>& c)
 : _p(c._p),
   _count(c._count)
{
   _count->increment();
}

template <typename T, typename C>
template <typename D> 
counted_collection<T,C>::counted_collection(const counted_collection<T,D>& c)
 : _p(c._p),
   _count(c._count)
{
   _count->increment();
}

/*
 * Assignment operators.
 */
template <typename T, typename C>
counted_collection<T,C>& counted_collection<T,C>::operator=(
   const counted_collection<T,C>& c)
{
   if (_p != c._p) {
      if (_count->decrement() == 0) {
         if (_p != NULL)
            _p->iter(iter_functors<T>::f_delete());
         delete _p;
         delete _count;
      }
      _p     = c._p;
      _count = c._count;
      _count->increment();
   }
   return *this;
}

template <typename T, typename C>
template <typename D> 
counted_collection<T,C>& counted_collection<T,C>::operator=(
   const counted_collection<T,D>& c)
{
   if (_p != c._p) {
      if (_count->decrement() == 0) {
         if (_p != NULL)
            _p->iter(iter_functors<T>::f_delete());
         delete _p;
         delete _count;
      }
      _p     = c._p;
      _count = c._count;
      _count->increment();
   }
   return *this;
}

/*
 * Destructor.
 * Delete the collection and its contents (if they are non-const items)
 * in the case that no other counted collections refer to it.
 */
template <typename T, typename C>
counted_collection<T,C>::~counted_collection() {
   if (_count->decrement() == 0) {
      if (_p != NULL)
         _p->iter(iter_functors<T>::f_delete());
      delete _p;
      delete _count;
   }
}

/*
 * Smart pointer dereference.
 * Return collection by reference.
 */
template <typename T, typename C>
inline C& counted_collection<T,C>::operator*() const {
   #if COLLECTIONS__POINTERS__COUNTED_COLLECTION__CHECK_DEREFERENCE
      if (_p == NULL)
         throw ex_null_pointer_dereference();
   #endif
   return *_p;
}

/*
 * Smart pointer dereference.
 * Return pointer to collection.
 */
template <typename T, typename C>
inline C* counted_collection<T,C>::operator->() const {
   #if COLLECTIONS__POINTERS__COUNTED_COLLECTION__CHECK_DEREFERENCE
      if (_p == NULL)
         throw ex_null_pointer_dereference();
   #endif
   return *_p;
}

/*
 * Get raw pointer to object.
 */
template <typename T, typename C>
inline C* counted_collection<T,C>::get() const {
   return _p;
}

/*
 * Release the collection (set the counted_collection to NULL).
 * Return a pointer to the collection.
 */
template <typename T, typename C>
C* counted_collection<T,C>::release() {
   T* p = _p;
   _p = NULL;
   if (_count->decrement() == 0)
      delete _count;
   _count = new counter(1);
   return p;
}
 
/*
 * Delete the collection if no other counted collections refer to it.
 * Set the counted_collection to refer to the new collection (defaults to NULL).
 */
template <typename T, typename C>
void counted_collection<T,C>::reset(C* p) {
   if (_p != p) {
      if (_count->decrement() == 0) {
         if (_p != NULL)
            _p->iter(iter_functors<T>::f_delete());
         delete _p;
         delete _count;
      }
      _p     = p;
      _count = new counter(1);
   }
}

/*
 * Automatic conversion to counted_collection of different type.
 */
template <typename T, typename C>
template <typename D> 
counted_collection<T,C>::operator counted_collection<T,D>() const {
   return counted_collection<T,D>(*this);
}

} /* namespace pointers */
} /* namespace collections */

#endif
