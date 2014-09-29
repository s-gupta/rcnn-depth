/*
 * Safe pointer.
 *
 * A safe pointer is a pointer that throws an exception on NULL dereference.
 * No ownership semantics are associated with safe pointers.
 */
#ifndef LANG__POINTERS__SAFE_PTR_HH
#define LANG__POINTERS__SAFE_PTR_HH

#include "config/safety.hh"
#include "lang/exceptions/ex_null_pointer_dereference.hh"
#include "lang/null.hh"

/*
 * Enable/disable null pointer dereference checking for safe pointers.
 */
#if CONFIG__SAFETY__CHECK_DEREFERENCE
   #define LANG__POINTERS__SAFE_PTR__CHECK_DEREFERENCE (true)
#else
   #define LANG__POINTERS__SAFE_PTR__CHECK_DEREFERENCE (false)
#endif

namespace lang {
namespace pointers {
/*
 * Imports.
 */
using lang::exceptions::ex_null_pointer_dereference;

/*
 * Safe pointer.
 */
template <typename T>
class safe_ptr {
public:
   /*
    * Friend classes.
    */
   template <typename U> friend class safe_ptr;
   
   /*
    * Constructor.
    */
   explicit safe_ptr(T* = NULL);

   /*
    * Copy constructors.
    * Set the copy to point to the same object as the original.
    */
   safe_ptr(const safe_ptr<T>&);
   template <typename U> safe_ptr(const safe_ptr<U>&);
   
   /*
    * Assignment operators.
    * Set the copy to point to the same object as the original.
    */
   safe_ptr<T>& operator=(const safe_ptr<T>&);
   template <typename U> safe_ptr<T>& operator=(const safe_ptr<U>&);

   /*
    * Destructor.
    */
   ~safe_ptr();

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
    * Set the safe_ptr to refer to the new object (defaults to NULL).
    */
   void reset(T* = NULL);

   /*
    * Automatic conversion to safe_ptr of different type.
    */
   template <typename U> operator safe_ptr<U>() const;

protected:
   T* _p;   /* pointer to object */
};

/*
 * Constructor.
 */
template <typename T>
safe_ptr<T>::safe_ptr(T* p)
 : _p(p)
{ }

/*
 * Copy constructors.
 */
template <typename T>
safe_ptr<T>::safe_ptr(const safe_ptr<T>& s)
 : _p(s._p)
{ }

template <typename T>
template <typename U> 
safe_ptr<T>::safe_ptr(const safe_ptr<U>& s)
 : _p(s._p)
{ }

/*
 * Assignment operators.
 */
template <typename T>
safe_ptr<T>& safe_ptr<T>::operator=(const safe_ptr<T>& s) {
   _p = s._p;
   return *this;
}

template <typename T>
template <typename U>
safe_ptr<T>& safe_ptr<T>::operator=(const safe_ptr<U>& s) {
   _p = s._p;
   return *this;
}

/*
 * Destructor.
 */
template <typename T>
safe_ptr<T>::~safe_ptr() {
   /* do nothing */
}

/*
 * Smart pointer dereference.
 * Return object by reference.
 */
template <typename T>
inline T& safe_ptr<T>::operator*() const {
   #if LANG__POINTERS__SAFE_PTR__CHECK_DEREFERENCE
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
inline T* safe_ptr<T>::operator->() const {
   #if LANG__POINTERS__SAFE_PTR__CHECK_DEREFERENCE
      if (_p == NULL)
         throw ex_null_pointer_dereference();
   #endif
   return _p;
}

/*
 * Get raw pointer to object.
 */
template <typename T>
inline T* safe_ptr<T>::get() const {
   return _p;
}

/*
 * Set the safe_ptr to refer to the new object (defaults to NULL).
 */
template <typename T>
void safe_ptr<T>::reset(T* p) {
   _p = p;
}

/*
 * Automatic conversion to safe_ptr of different type.
 */
template <typename T>
template <typename U>
safe_ptr<T>::operator safe_ptr<U>() const {
   return safe_ptr<U>(*this);
}

} /* namespace pointers */
} /* namespace lang */

#endif
