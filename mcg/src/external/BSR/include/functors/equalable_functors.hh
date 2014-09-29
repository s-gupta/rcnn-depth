/*
 * Equalable functors.
 */
#ifndef FUNCTORS__EQUALABLE_FUNCTORS_HH
#define FUNCTORS__EQUALABLE_FUNCTORS_HH

#include "interfaces/equalable.hh"

namespace functors {
/*
 * Imports.
 */
using interfaces::equalable;

/*
 * Abstract functor for testing equality.
 * Note that the functor is the same regardless of the type const qualifier.
 */
template <typename T>
class equalable_functor {
public:
   virtual ~equalable_functor() { }
   virtual bool operator()(const T&, const T&) const = 0;
};

template <typename T>
class equalable_functor<const T> : public equalable_functor<T> { };

/*
 * Functor for testing equality between equalable items.
 * Call the underlying equality test function.
 */
template <typename T>
class equal_functor : public equalable_functor<const T> {
public:
   bool operator()(const T& t0, const T& t1) const {
      return equalable<T>::is_equal(t0, t1);
   }
};

template <typename T>
class equal_functor<const T> : public equal_functor<T> { };

/*
 * Functor that returns inverse of equality test (not equal test).
 */
template <typename T>
class equal_functor_inverse : public equalable_functor<const T> {
public:
   equal_functor_inverse() : _f(equal_functor<T>()) { }
   bool operator()(const T& t0, const T& t1) const {
      return (!(_f(t0, t1)));
   }
private:
   const equal_functor<T> _f;
};

template <typename T>
class equal_functor_inverse<const T> : public equal_functor_inverse<T> { };

/*
 * Functor that always returns true for equality test.
 */
template <typename T>
class equal_functor_true : public equalable_functor<const T> {
public:
   bool operator()(const T& t0, const T& t1) const {
      return true;
   }
};

template <typename T>
class equal_functor_true<const T> : public equal_functor_true<T> { };

/*
 * Functor that always returns false for equality test.
 */
template <typename T>
class equal_functor_false : public equalable_functor<const T> {
public:
   bool operator()(const T& t0, const T& t1) const {
      return false;
   }
};

template <typename T>
class equal_functor_false<const T> : public equal_functor_false<T> { };

/*
 * Functor for testing if items are both references to the same object in 
 * memory.
 */
template <typename T>
class equal_functor_memloc : public equalable_functor<const T> {
public:
   bool operator()(const T& t0, const T& t1) const {
      const T* p0 = &t0;
      const T* p1 = &t1;
      return (p0 == p1);
   }
};

template <typename T>
class equal_functor_memloc<const T> : public equal_functor_memloc<T> { };

/*
 * Functor for invoking equality test through pointers to the items.
 */
template <typename T, typename Ptr>
class ptr_equal_functor : public equalable_functor<const Ptr> {
public:
   ptr_equal_functor() : _f() { }
   bool operator()(const Ptr& p0, const Ptr& p1) const {
      return _f(*p0, *p1);
   }
private:
   const equal_functor<T> _f;
};

template <typename T, typename Ptr>
class ptr_equal_functor<T, const Ptr> : public ptr_equal_functor<T,Ptr> { };

/*
 * Functor for invoking inverse of equality test through pointers to the items.
 */
template <typename T, typename Ptr>
class ptr_equal_functor_inverse : public equalable_functor<const Ptr> {
public:
   ptr_equal_functor_inverse() : _f() { }
   bool operator()(const Ptr& p0, const Ptr& p1) const {
      return _f(*p0, *p1);
   }
private:
   const equal_functor_inverse<T> _f;
};

template <typename T, typename Ptr>
class ptr_equal_functor_inverse<T, const Ptr>
 : public ptr_equal_functor_inverse<T,Ptr> { };

/*
 * Functor that always returns true for equality test through pointers to the
 * items.
 */
template <typename T, typename Ptr>
class ptr_equal_functor_true : public equalable_functor<const Ptr> {
public:
   bool operator()(const Ptr& p0, const Ptr& p1) const {
      return true;
   }
};

template <typename T, typename Ptr>
class ptr_equal_functor_true<T, const Ptr>
 : public ptr_equal_functor_true<T,Ptr> { };

/*
 * Functor that always returns false for equality test through pointers to the
 * items.
 */
template <typename T, typename Ptr>
class ptr_equal_functor_false : public equalable_functor<const Ptr> {
public:
   bool operator()(const Ptr& p0, const Ptr& p1) const {
      return false;
   }
};

template <typename T, typename Ptr>
class ptr_equal_functor_false<T, const Ptr>
 : public ptr_equal_functor_false<T,Ptr> { };

/*
 * Functor for testing if items are both references to the same object in 
 * memory through pointers to the items.
 */
template <typename T, typename Ptr>
class ptr_equal_functor_memloc : public equalable_functor<const Ptr> {
public:
   ptr_equal_functor_memloc() : _f() { }
   bool operator()(const Ptr& p0, const Ptr& p1) const {
      return _f(*p0, *p1);
   }
private:
   const equal_functor_memloc<T> _f;
};

template <typename T, typename Ptr>
class ptr_equal_functor_memloc<T, const Ptr>
 : public ptr_equal_functor_memloc<T,Ptr> { };

/*
 * Globally accessible set of default equality test functors.
 */
template <typename T>
class equal_functors {
public:
   static const equal_functor<T>&         f_equal();
   static const equal_functor_inverse<T>& f_not_equal();
   static const equal_functor_true<T>&    f_true();
   static const equal_functor_false<T>&   f_false();
   static const equal_functor_memloc<T>&  f_equal_memloc();
};

template <typename T>
const equal_functor<T>& equal_functors<T>::f_equal() {
   static const equal_functor<T>* f = new equal_functor<T>();
   return *f;
}

template <typename T>
const equal_functor_inverse<T>& equal_functors<T>::f_not_equal() {
   static const equal_functor_inverse<T>* f = new equal_functor_inverse<T>();
   return *f;
}

template <typename T>
const equal_functor_true<T>& equal_functors<T>::f_true() {
   static const equal_functor_true<T>* f = new equal_functor_true<T>();
   return *f;
}

template <typename T>
const equal_functor_false<T>& equal_functors<T>::f_false() {
   static const equal_functor_false<T>* f = new equal_functor_false<T>();
   return *f;
}

template <typename T>
const equal_functor_memloc<T>& equal_functors<T>::f_equal_memloc() {
   static const equal_functor_memloc<T>* f = new equal_functor_memloc<T>();
   return *f;
}

/*
 * Globally accessible set of default pointer equality test functors.
 */
template <typename T, typename Ptr>
class ptr_equal_functors {
public:
   static const ptr_equal_functor<T,Ptr>&         f_equal();
   static const ptr_equal_functor_inverse<T,Ptr>& f_not_equal();
   static const ptr_equal_functor_true<T,Ptr>&    f_true();
   static const ptr_equal_functor_false<T,Ptr>&   f_false();
   static const ptr_equal_functor_memloc<T,Ptr>&  f_equal_memloc();
};

template <typename T, typename Ptr>
const ptr_equal_functor<T,Ptr>& ptr_equal_functors<T,Ptr>::f_equal() {
   static const ptr_equal_functor<T,Ptr>* f
      = new ptr_equal_functor<T,Ptr>();
   return *f;
}

template <typename T, typename Ptr>
const ptr_equal_functor_inverse<T,Ptr>& ptr_equal_functors<T,Ptr>::f_not_equal() {
   static const ptr_equal_functor_inverse<T,Ptr>* f
      = new ptr_equal_functor_inverse<T,Ptr>();
   return *f;
}

template <typename T, typename Ptr>
const ptr_equal_functor_true<T,Ptr>& ptr_equal_functors<T,Ptr>::f_true() {
   static const ptr_equal_functor_true<T,Ptr>* f
      = new ptr_equal_functor_true<T,Ptr>();
   return *f;
}

template <typename T, typename Ptr>
const ptr_equal_functor_false<T,Ptr>& ptr_equal_functors<T,Ptr>::f_false() {
   static const ptr_equal_functor_false<T,Ptr>* f
      = new ptr_equal_functor_false<T,Ptr>();
   return *f;
}

template <typename T, typename Ptr>
const ptr_equal_functor_memloc<T,Ptr>& ptr_equal_functors<T,Ptr>::f_equal_memloc() {
   static const ptr_equal_functor_memloc<T,Ptr>* f
      = new ptr_equal_functor_memloc<T,Ptr>();
   return *f;
}

/*
 * Functor for testing pointer equality.
 */
template <typename T>
class equal_functor<T*>
 : public equalable_functor<T* const> {
public:
   bool operator()(T* const& t0, T* const& t1) const { 
      return (t0 == t1);
   }
};
 
/*
 * Functors for testing equality of built-in data types.
 */
template <>
class equal_functor<bool>
 : public equalable_functor<const bool> {
public:
   bool operator()(const bool& b0, const bool& b1) const { 
      return (b0 == b1);
   }
};

template <>
class equal_functor<char>
 : public equalable_functor<const char> {
public:
   bool operator()(const char& c0, const char& c1) const { 
      return (c0 == c1);
   }
};

template <>
class equal_functor<unsigned char>
 : public equalable_functor<const unsigned char> {
public:
   bool operator()(const unsigned char& c0, const unsigned char& c1) const { 
      return (c0 == c1);
   }
};

template <>
class equal_functor<short>
 : public equalable_functor<const short> {
public:
   bool operator()(const short& i0, const short& i1) const { 
      return (i0 == i1);
   }
};

template <>
class equal_functor<unsigned short>
 : public equalable_functor<const unsigned short> {
public:
   bool operator()(const unsigned short& i0, const unsigned short& i1) const { 
      return (i0 == i1);
   }
};

template <>
class equal_functor<int>
 : public equalable_functor<const int> {
public:
   bool operator()(const int& i0, const int& i1) const { 
      return (i0 == i1);
   }
};

template <>
class equal_functor<unsigned int>
 : public equalable_functor<const unsigned int> {
public:
   bool operator()(const unsigned int& i0, const unsigned int& i1) const { 
      return (i0 == i1);
   }
};

template <>
class equal_functor<long>
 : public equalable_functor<const long> {
public:
   bool operator()(const long& i0, const long& i1) const { 
      return (i0 == i1);
   } 
};

template <>
class equal_functor<unsigned long>
 : public equalable_functor<const unsigned long> {
public:
   bool operator()(const unsigned long& i0, const unsigned long& i1) const { 
      return (i0 == i1);
   }
};

template <>
class equal_functor<long long>
 : public equalable_functor<const long long> {
public:
   bool operator()(const long long& i0, const long long& i1) const { 
      return (i0 == i1);
   } 
};

template <>
class equal_functor<unsigned long long>
 : public equalable_functor<const unsigned long long> {
public:
   bool operator()(
      const unsigned long long& i0,
      const unsigned long long& i1) const 
   { 
      return (i0 == i1);
   }
};

template <>
class equal_functor<float>
 : public equalable_functor<const float> {
public:
   bool operator()(const float& f0, const float& f1) const { 
      return (f0 == f1);
   } 
};

template <>
class equal_functor<double>
 : public equalable_functor<const double> {
public:
   bool operator()(const double& d0, const double& d1) const { 
      return (d0 == d1);
   }
};

template <>
class equal_functor<long double>
 : public equalable_functor<const long double> {
public:
   bool operator()(const long double& d0, const long double& d1) const { 
      return (d0 == d1);
   }
};

} /* namespace functors */

#endif
