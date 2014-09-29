/*
 * Comparable functors.
 */
#ifndef FUNCTORS__COMPARABLE_FUNCTORS_HH
#define FUNCTORS__COMPARABLE_FUNCTORS_HH

#include "interfaces/comparable.hh"

namespace functors {
/*
 * Imports.
 */
using interfaces::comparable;

/*
 * Abstract functor for comparing items.
 * Note that the functor is the same regardless of the type const qualifier.
 */
template <typename T>
class comparable_functor {
public:
   virtual ~comparable_functor() { }
   virtual int operator()(const T&, const T&) const = 0;
};

template <typename T>
class comparable_functor<const T> : public comparable_functor<T> { };

/*
 * Functor for comparing comparable items.
 * Call the underlying compare function.
 */
template <typename T>
class compare_functor : public comparable_functor<const T> {
public:
   int operator()(const T& t0, const T& t1) const {
      return comparable<T>::compare(t0, t1);
   }
};

template <typename T>
class compare_functor<const T> : public compare_functor<T> { };

/*
 * Functor that reverses normal comparison order.
 */
template <typename T>
class compare_functor_reverse : public comparable_functor<const T> {
public:
   compare_functor_reverse() : _f() { }
   int operator()(const T& t0, const T& t1) const {
      return (-(_f(t0, t1)));
   }
private:
   const compare_functor<T> _f;
};

template <typename T>
class compare_functor_reverse<const T> : public compare_functor_reverse<T> { };

/*
 * Functor that always returns zero (equal) for comparison test.
 */
template <typename T>
class compare_functor_zero : public comparable_functor<const T> {
public:
   int operator()(const T& t0, const T& t1) const {
      return 0;
   }
};

template <typename T>
class compare_functor_zero<const T> : public compare_functor_zero<T> { };

/*
 * Functor that compares location of items in memory.
 * Items are equal iff they are both references to the same object in memory.
 */
template <typename T>
class compare_functor_memloc : public comparable_functor<const T> {
public:
   int operator()(const T& t0, const T& t1) const {
      const T* p0 = &t0;
      const T* p1 = &t1;
      return (p0 < p1) ? (-1) : ((p0 > p1) ? 1 : 0);
   }
};

template <typename T>
class compare_functor_memloc<const T> : public compare_functor_memloc<T> { };

/*
 * Functor for invoking item comparison through pointers to the items. 
 */
template <typename T, typename Ptr>
class ptr_compare_functor : public comparable_functor<const Ptr> {
public:
   ptr_compare_functor() : _f() { }
   int operator()(const Ptr& p0, const Ptr& p1) const {
      return _f(*p0, *p1);
   }
private:
   const compare_functor<T> _f;
};

template <typename T, typename Ptr>
class ptr_compare_functor<T, const Ptr>
 : public ptr_compare_functor<T,Ptr> { };

/*
 * Functor for invoking reverse item comparison through pointers to the items.
 */
template <typename T, typename Ptr>
class ptr_compare_functor_reverse : public comparable_functor<const Ptr> {
public:
   ptr_compare_functor_reverse() : _f() { }
   int operator()(const Ptr& p0, const Ptr& p1) const {
      return _f(*p0, *p1);
   }
private:
   const compare_functor_reverse<T> _f;
};

template <typename T, typename Ptr>
class ptr_compare_functor_reverse<T, const Ptr>
 : public ptr_compare_functor_reverse<T,Ptr> { };

/*
 * Functor that always returns zero (equal) for comparison of pointers to items.
 */
template <typename T, typename Ptr>
class ptr_compare_functor_zero : public comparable_functor<const Ptr> {
public:
   int operator()(const Ptr& p0, const Ptr& p1) const {
      return 0;
   }
};

template <typename T, typename Ptr>
class ptr_compare_functor_zero<T, const Ptr>
 : public ptr_compare_functor_zero<T,Ptr> { };

/*
 * Functor that compares location of items in memory through pointers to the
 * items.  Items are equal iff they are both references to the same object in
 * memory.
 */
template <typename T, typename Ptr>
class ptr_compare_functor_memloc : public comparable_functor<const Ptr> {
public:
   ptr_compare_functor_memloc() : _f() { }
   int operator()(const Ptr& p0, const Ptr& p1) const {
      return _f(*p0, *p1);
   }
private:
   const compare_functor_memloc<T> _f;
};

template <typename T, typename Ptr>
class ptr_compare_functor_memloc<T, const Ptr>
 : public ptr_compare_functor_memloc<T,Ptr> { };

/*
 * Globally accessible set of default comparison functors.
 */
template <typename T>
class compare_functors {
public:
   static const compare_functor<T>&         f_compare();
   static const compare_functor_reverse<T>& f_compare_reverse();
   static const compare_functor_zero<T>&    f_zero();
   static const compare_functor_memloc<T>&  f_compare_memloc();
};

template <typename T>
const compare_functor<T>& compare_functors<T>::f_compare() {
   static const compare_functor<T>* f = new compare_functor<T>();
   return *f;
}

template <typename T>
const compare_functor_reverse<T>& compare_functors<T>::f_compare_reverse() {
   static const compare_functor_reverse<T>* f = new compare_functor_reverse<T>();
   return *f;
}

template <typename T>
const compare_functor_zero<T>& compare_functors<T>::f_zero() {
   static const compare_functor_zero<T>* f = new compare_functor_zero<T>();
   return *f;
}

template <typename T>
const compare_functor_memloc<T>& compare_functors<T>::f_compare_memloc() {
   static const compare_functor_memloc<T>* f = new compare_functor_memloc<T>();
   return *f;
}

/*
 * Globally accessible set of default pointer comparison functors.
 */
template <typename T, typename Ptr>
class ptr_compare_functors {
public:
   static const ptr_compare_functor<T,Ptr>&         f_compare();
   static const ptr_compare_functor_reverse<T,Ptr>& f_compare_reverse();
   static const ptr_compare_functor_zero<T,Ptr>&    f_zero();
   static const ptr_compare_functor_memloc<T,Ptr>&  f_compare_memloc();
};

template <typename T, typename Ptr>
const ptr_compare_functor<T,Ptr>& ptr_compare_functors<T,Ptr>::f_compare() {
   static const ptr_compare_functor<T,Ptr>* f
      = new ptr_compare_functor<T,Ptr>();
   return *f;
}

template <typename T, typename Ptr>
const ptr_compare_functor_reverse<T,Ptr>& ptr_compare_functors<T,Ptr>::f_compare_reverse() {
   static const ptr_compare_functor_reverse<T,Ptr>* f
      = new ptr_compare_functor_reverse<T,Ptr>();
   return *f;
}

template <typename T, typename Ptr>
const ptr_compare_functor_zero<T,Ptr>& ptr_compare_functors<T,Ptr>::f_zero() {
   static const ptr_compare_functor_zero<T,Ptr>* f
      = new ptr_compare_functor_zero<T,Ptr>();
   return *f;
}

template <typename T, typename Ptr>
const ptr_compare_functor_memloc<T,Ptr>& ptr_compare_functors<T,Ptr>::f_compare_memloc() {
   static const ptr_compare_functor_memloc<T,Ptr>* f
      = new ptr_compare_functor_memloc<T,Ptr>();
   return *f;
}

/*
 * Functor for comparing pointers.
 */
template <typename T>
class compare_functor<T*> : public comparable_functor<T* const> {
public:
   int operator()(T* const& t0, T* const& t1) const { 
      return (t0 < t1) ? (-1) : ((t0 > t1) ? 1 : 0);
   }
};

/*
 * Functors for comparing built-in data types.
 */
template <>
class compare_functor<bool>
 : public comparable_functor<const bool> {
public:
   int operator()(const bool& b0, const bool& b1) const { 
      return (static_cast<int>(b0)) - (static_cast<int>(b1));
   }
};

template <>
class compare_functor<char>
 : public comparable_functor<const char> {
public:
   int operator()(const char& c0, const char& c1) const { 
      return (static_cast<int>(c0)) - (static_cast<int>(c1));
   }
};

template <>
class compare_functor<unsigned char>
 : public comparable_functor<const unsigned char> {
public:
   int operator()(const unsigned char& c0, const unsigned char& c1) const { 
      return (static_cast<int>(c0)) - (static_cast<int>(c1));
   }
};

template <>
class compare_functor<short>
 : public comparable_functor<const short> {
public:
   int operator()(const short& i0, const short& i1) const { 
      return (static_cast<int>(i0)) - (static_cast<int>(i1));
   }
};

template <>
class compare_functor<unsigned short>
 : public comparable_functor<const unsigned short> {
public:
   int operator()(const unsigned short& i0, const unsigned short& i1) const { 
      return (static_cast<int>(i0)) - (static_cast<int>(i1));
   }
};

template <>
class compare_functor<int>
 : public comparable_functor<const int> {
public:
   int operator()(const int& i0, const int& i1) const { 
      return (i0 < i1) ? (-1) : ((i0 > i1) ? 1 : 0); 
   }
};

template <>
class compare_functor<unsigned int>
 : public comparable_functor<const unsigned int> {
public:
   int operator()(const unsigned int& i0, const unsigned int& i1) const { 
      return (i0 < i1) ? (-1) : ((i0 > i1) ? 1 : 0); 
   }
};

template <>
class compare_functor<long>
 : public comparable_functor<const long> {
public:
   int operator()(const long& i0, const long& i1) const { 
      return (i0 < i1) ? (-1) : ((i0 > i1) ? 1 : 0); 
   } 
};

template <>
class compare_functor<unsigned long>
 : public comparable_functor<const unsigned long> {
public:
   int operator()(const unsigned long& i0, const unsigned long& i1) const { 
      return (i0 < i1) ? (-1) : ((i0 > i1) ? 1 : 0); 
   }
};

template <>
class compare_functor<long long>
 : public comparable_functor<const long long> {
public:
   int operator()(const long long& i0, const long long& i1) const { 
      return (i0 < i1) ? (-1) : ((i0 > i1) ? 1 : 0); 
   } 
};

template <>
class compare_functor<unsigned long long>
 : public comparable_functor<const unsigned long long> {
public:
   int operator()(
      const unsigned long long& i0,
      const unsigned long long& i1) const 
   { 
      return (i0 < i1) ? (-1) : ((i0 > i1) ? 1 : 0); 
   }
};

template <>
class compare_functor<float>
 : public comparable_functor<const float> {
public:
   int operator()(const float& f0, const float& f1) const { 
      return (f0 < f1) ? (-1) : ((f0 > f1) ? 1 : 0); 
   } 
};

template <>
class compare_functor<double>
 : public comparable_functor<const double> {
public:
   int operator()(const double& d0, const double& d1) const { 
      return (d0 < d1) ? (-1) : ((d0 > d1) ? 1 : 0); 
   }
};

template <>
class compare_functor<long double>
 : public comparable_functor<const long double> {
public:
   int operator()(const long double& d0, const long double& d1) const { 
      return (d0 < d1) ? (-1) : ((d0 > d1) ? 1 : 0); 
   }
};

} /* namespace functors */

#endif
