/*
 * Counted pointer.
 *
 * The set of counted pointers that point to an object own that object and 
 * delete it upon the destruction of the last counted pointer that refers
 * to it (unless ownership is explicitly released by the user).  Copying a 
 * counted_ptr increments the reference count of the object.
 */
#ifndef LANG__POINTERS__COUNTED_PTR_HH
#define LANG__POINTERS__COUNTED_PTR_HH

#include "concurrent/threads/synchronization/counter.hh"
#include "config/safety.hh"
#include "lang/exceptions/ex_null_pointer_dereference.hh"
#include "lang/null.hh"

/*
 * Enable/disable null pointer dereference checking for counted pointers.
 */
#if CONFIG__SAFETY__CHECK_DEREFERENCE
   #define LANG__POINTERS__COUNTED_PTR__CHECK_DEREFERENCE (true)
#else
   #define LANG__POINTERS__COUNTED_PTR__CHECK_DEREFERENCE (false)
#endif

namespace lang {
namespace pointers {
/*
 * Imports.
 */
using concurrent::threads::synchronization::counter;
using lang::exceptions::ex_null_pointer_dereference;

/*
 * Counted pointer.
 */
template <typename T>
class counted_ptr {
public:
   /*
    * Friend classes.
    */
   template <typename U> friend class counted_ptr;
   
   /*
    * Constructor.
    */
   explicit counted_ptr(T* = NULL);

   /*
    * Copy constructors.
    * Set the copy to point to the same object as the original.
    */
   counted_ptr(const counted_ptr<T>&);
   template <typename U> counted_ptr(const counted_ptr<U>&);
   
   /*
    * Assignment operators.
    * Set the copy to point to the same object as the original.
    */
   counted_ptr<T>& operator=(const counted_ptr<T>&);
   template <typename U> counted_ptr<T>& operator=(const counted_ptr<U>&);

   /*
    * Destructor.
    * Delete the object if no other counted pointers refer to it.
    */
   ~counted_ptr();

   /*
    * Smart pointer dereference.
    * Return object by reference.
    */
   inline T& operator*() const;

   /*
    * Smart pointer dereference.
    * Return pointer to object.
    */
   inline T* operator->() const;

   /*
    * Get raw pointer to object.
    */
   inline T* get() const;

   /*
    * Release the object (set the counted_ptr to NULL).
    * Return a pointer to the object.
    */
   T* release();
    
   /*
    * Delete the object if no other counted pointers refer to it.
    * Set the counted_ptr to refer to the new object (defaults to NULL).
    */
   void reset(T* = NULL);

   /*
    * Automatic conversion to counted_ptr of different type.
    */
   template <typename U> operator counted_ptr<U>() const;

protected:
   T*       _p;      /* pointer to object */
   counter* _count;  /* reference count for object */
};

/*
 * Constructor.
 */
template <typename T>
counted_ptr<T>::counted_ptr(T* p)
 : _p(p),
   _count(new counter(1))
{ }

/*
 * Copy constructors.
 */
template <typename T>
counted_ptr<T>::counted_ptr(const counted_ptr<T>& c)
 : _p(c._p),
   _count(c._count)
{ 
   _count->increment();
}

template <typename T>
template <typename U> 
counted_ptr<T>::counted_ptr(const counted_ptr<U>& c)
 : _p(c._p),
   _count(c._count)
{
   _count->increment();
}

/*
 * Assignment operators.
 */
template <typename T>
counted_ptr<T>& counted_ptr<T>::operator=(const counted_ptr<T>& c) {
   if (_p != c._p) {
      if (_count->decrement() == 0) {
         delete _p;
         delete _count;
      }
      _p     = c._p;
      _count = c._count;
      _count->increment();
   }
   return *this;
}

template <typename T>
template <typename U>
counted_ptr<T>& counted_ptr<T>::operator=(const counted_ptr<U>& c) {
   if (_p != c._p) {
      if (_count->decrement() == 0) {
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
 * Delete the object if no other counted pointers refer to it.
 */
template <typename T>
counted_ptr<T>::~counted_ptr() {
   if (_count->decrement() == 0) {
      delete _p;
      delete _count;
   }
}

/*
 * Smart pointer dereference.
 * Return object by reference.
 */
template <typename T>
inline T& counted_ptr<T>::operator*() const {
   #if LANG__POINTERS__COUNTED_PTR__CHECK_DEREFERENCE
      if (_p == NULL)
         throw ex_null_pointer_dereference();
   #endif
   return *_p;
}

/*
 * Smart pointer dereference.
 * Return pointer to object.
 */
template <typename T>
inline T* counted_ptr<T>::operator->() const {
   #if LANG__POINTERS__COUNTED_PTR__CHECK_DEREFERENCE
      if (_p == NULL)
         throw ex_null_pointer_dereference();
   #endif
   return _p;
}

/*
 * Get raw pointer to object.
 */
template <typename T>
inline T* counted_ptr<T>::get() const {
   return _p;
}

/*
 * Release the object (set the counted_ptr to NULL).
 * Return a pointer to the object.
 */
template <typename T>
T* counted_ptr<T>::release() {
   T* p = _p;
   _p = NULL;
   if (_count->decrement() == 0)
      delete _count;
   _count = new counter(1);
   return p;
}
 
/*
 * Delete the object if no other counted pointers refer to it.
 * Set the counted_ptr to refer to the new object (defaults to NULL).
 */
template <typename T>
void counted_ptr<T>::reset(T* p) {
   if (_p != p) {
      if (_count->decrement() == 0) {
         delete _p;
         delete _count;
      }
      _p     = p;
      _count = new counter(1);
   }
}

/*
 * Automatic conversion to counted_ptr of different type.
 */
template <typename T>
template <typename U>
counted_ptr<T>::operator counted_ptr<U>() const {
   return counted_ptr<U>(*this);
}

} /* namespace pointers */
} /* namespace lang */

#endif
