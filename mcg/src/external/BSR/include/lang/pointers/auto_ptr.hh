/*
 * Auto pointer.
 *
 * An auto pointer owns the object it holds a pointer to and deletes this 
 * object upon its own destruction.  Copying an auto_ptr copies the pointer
 * and transfers ownership.  For any given object, at most one auto_ptr should
 * own it.
 */
#ifndef LANG__POINTERS__AUTO_PTR_HH
#define LANG__POINTERS__AUTO_PTR_HH

#include "config/safety.hh"
#include "lang/exceptions/ex_null_pointer_dereference.hh"
#include "lang/null.hh"

/*
 * Enable/disable null pointer dereference checking for auto pointers.
 */
#ifdef CONFIG__SAFETY__CHECK_DEREFERENCE
   #define LANG__POINTERS__AUTO_PTR__CHECK_DEREFERENCE (true)
#else
   #define LANG__POINTERS__AUTO_PTR__CHECK_DEREFERENCE (false)
#endif

namespace lang {
namespace pointers {
/*
 * Imports.
 */
using lang::exceptions::ex_null_pointer_dereference;

/*
 * Declare auto_ptr class.
 */
template <typename T>
class auto_ptr;

/*
 * Auto pointer reference.
 *
 * This wrapper class provides auto_ptr with reference semantics.
 */
template <typename T>
class auto_ptr_ref {
public:
   /*
    * Friend classes.
    */
   template <typename U> friend class auto_ptr;
   
protected:
   /* 
    * Constructor.
    */
   explicit inline auto_ptr_ref(T*);
   
   T* _p;   /* pointer */
};

/*
 * Auto pointer reference constructor.
 */
template <typename T>
inline auto_ptr_ref<T>::auto_ptr_ref(T* p) : _p(p) { }

/*
 * Auto pointer.
 */
template <typename T>
class auto_ptr {
public:
   /*
    * Constructor.
    */
   explicit inline auto_ptr(T* = NULL);

   /*
    * Copy constructors.
    *
    * The constructed auto_ptr assumes ownership
    * and the original auto_ptr is set to NULL.
    *
    * Note that a const auto_ptr cannot be copied, 
    * since it is not allowed to transfer ownership.
    */
   auto_ptr(auto_ptr<T>&);
   template <typename U> auto_ptr(auto_ptr<U>&);

   /*
    * Assignment operators.
    *
    * The assigned auto_ptr assumes ownership
    * and the original auto_ptr is set to NULL.
    *
    * Any object previously owned by the assigned
    * auto_ptr is deleted.
    * 
    * Note that a const auto_ptr cannot appear on  
    * the right hand side of an assignment, since 
    * it is not allowed to transfer ownership.
    */
   auto_ptr<T>& operator=(auto_ptr<T>&);
   template <typename U> auto_ptr<T>& operator=(auto_ptr<U>&);

   /*
    * Destructor.
    * Delete the owned object (if any).
    */
   inline ~auto_ptr();
   
   /*
    * Smart pointer dereference.
    * Return owned object by reference.
    */
   inline T& operator*() const;

   /*
    * Smart pointer dereference.
    * Return pointer to the owned object.
    */
   inline T* operator->() const;

   /*
    * Get raw pointer to owned object.
    * The auto_ptr still owns the object.
    */
   inline T* get() const;

   /*
    * Release the owned object (set the auto_ptr to NULL).
    * Return a pointer to the object.
    */
   T* release();

   /*
    * Delete the owned object and set the auto_ptr to 
    * own the given pointer (defaults to NULL) instead.
    */
   void reset(T* = NULL);
   
   /*
    * Automatic conversions to/from auto_ptr_ref.
    *
    * These conversions are necessary in order to allow
    * functions to return auto pointers by value (since 
    * copying a const auto_ptr is forbidden).
    */
   auto_ptr(auto_ptr_ref<T>);
   auto_ptr<T>& operator=(auto_ptr_ref<T>);
   template <typename U> operator auto_ptr_ref<U>();

   /*
    * Automatic conversion to auto_ptr of different type.
    */
   template <typename U> operator auto_ptr<U>();
         
protected:
   T* _p;   /* pointer to owned object */
};

/*
 * Constructor.
 */
template <typename T>
inline auto_ptr<T>::auto_ptr(T* p) : _p(p) { }

/*
 * Copy constructors.
 */
template <typename T>
auto_ptr<T>::auto_ptr(auto_ptr<T>& a) : _p(a.release()) { }

template <typename T>
template <typename U>
auto_ptr<T>::auto_ptr(auto_ptr<U>& a) : _p(a.release()) { }

/*
 * Assignment operators.
 */
template <typename T>
auto_ptr<T>& auto_ptr<T>::operator=(auto_ptr<T>& a) {
   this->reset(a.release());
   return *this;
}

template <typename T>
template <typename U>
auto_ptr<T>& auto_ptr<T>::operator=(auto_ptr<U>& a) {
   this->reset(a.release());
   return *this;
}

/*
 * Destructor.
 * Delete the owned object (if any).
 */
template <typename T>
inline auto_ptr<T>::~auto_ptr() {
   delete _p;
}

/*
 * Smart pointer dereference.
 * Return owned object by reference.
 */
template <typename T>
inline T& auto_ptr<T>::operator*() const {
   #ifdef LANG__POINTERS__AUTO_PTR__CHECK_DEREFERENCE
      if (_p == NULL)
         throw ex_null_pointer_dereference();
   #endif
   return *_p;
}

/*
 * Smart pointer dereference.
 * Return pointer to the owned object.
 */
template <typename T>
inline T* auto_ptr<T>::operator->() const {
   #ifdef LANG__POINTERS__AUTO_PTR__CHECK_DEREFERENCE
      if (_p == NULL)
         throw ex_null_pointer_dereference();
   #endif
   return _p;
}

/*
 * Get raw pointer to owned object.
 * The auto_ptr still owns the object.
 */
template <typename T>
inline T* auto_ptr<T>::get() const {
   return _p;
}

/*
 * Release the owned object (set the auto_ptr to NULL).
 * Return a pointer to the object.
 */
template <typename T>
T* auto_ptr<T>::release() {
   T* p = _p;
   _p = NULL;
   return p;
}

/*
 * Delete the owned object and set the auto_ptr to 
 * own the given pointer (defaults to NULL) instead.
 */
template <typename T>
void auto_ptr<T>::reset(T* p) {
   if (_p != p) {
      delete _p;
      _p = p;
   }
}

/*
 * Automatic conversions to/from auto_ptr_ref.
 */
template <typename T>
auto_ptr<T>::auto_ptr(auto_ptr_ref<T> r) : _p(r._p) { }

template <typename T>
auto_ptr<T>& auto_ptr<T>::operator=(auto_ptr_ref<T> r) {
   if (r._p != _p) {
      delete _p;
      _p = r._p;
   }
   return *this;
}

template <typename T>
template <typename U>
auto_ptr<T>::operator auto_ptr_ref<U>() {
   return auto_ptr_ref<U>(this->release());
}

/*
 * Automatic conversion to auto_ptr of different type.
 */
template <typename T>
template <typename U>
auto_ptr<T>::operator auto_ptr<U>() {
   return auto_ptr<U>(this->release());
}

} /* namespace pointers */
} /* namespace lang */

#endif
