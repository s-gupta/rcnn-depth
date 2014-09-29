/*
 * Distanceable functors.
 */
#ifndef FUNCTORS__DISTANCEABLE_FUNCTORS_HH
#define FUNCTORS__DISTANCEABLE_FUNCTORS_HH

#include "interfaces/distanceable.hh"

namespace functors {
/*
 * Imports.
 */
using interfaces::distanceable;

/*
 * Abstract functor for computing distance.
 * Note that the functor is the same regardless of the type const qualifier.
 */
template <typename T, typename V>
class distanceable_functor {
public:
   virtual ~distanceable_functor() { }
   virtual V operator()(const T&, const T&) const = 0;
};

template <typename T, typename V>
class distanceable_functor<const T, V> : public distanceable_functor<T,V> { };

/*
 * Functor for computing distance between distanceable items.
 * Call the underlying distance computation function.
 */
template <typename T, typename V>
class distance_functor : public distanceable_functor<const T, V> {
public:
   V operator()(const T& t0, const T& t1) const {
      return distanceable<T,V>::distance(t0, t1);
   }
};

template <typename T, typename V>
class distance_functor<const T, V> : public distance_functor<T,V> { };

/*
 * Functor for invoking distance computation through pointers to items.
 */
template <typename T, typename V, typename Ptr>
class ptr_distance_functor : public distanceable_functor<const Ptr, V> {
public:
   ptr_distance_functor() : _f() { }
   V operator()(const Ptr& p0, const Ptr& p1) const {
      return _f(*p0, *p1);
   }
private:
   const distance_functor<T,V> _f;
};

template <typename T, typename V, typename Ptr>
class ptr_distance_functor<T, V, const Ptr>
 : public ptr_distance_functor<T,V,Ptr> { };

/*
 * Globally accessible set of default distance functors.
 */
template <typename T, typename V>
class distance_functors {
public:
   static const distance_functor<T,V>& f_distance();
};

template <typename T, typename V>
const distance_functor<T,V>& distance_functors<T,V>::f_distance() {
   static const distance_functor<T,V>* f = new distance_functor<T,V>();
   return *f;
}

/*
 * Globally accessible set of default pointer distance functors.
 */
template <typename T, typename V, typename Ptr>
class ptr_distance_functors {
public:
   static const ptr_distance_functor<T,V,Ptr>& f_distance();
};

template <typename T, typename V, typename Ptr>
const ptr_distance_functor<T,V,Ptr>& ptr_distance_functors<T,V,Ptr>::f_distance() {
   static const ptr_distance_functor<T,V,Ptr>* f
      = new ptr_distance_functor<T,V,Ptr>();
   return *f;
}

/*
 * Functors for computing distances for built-in data types.
 */
template <>
class distance_functor<bool, bool>
 : public distanceable_functor<const bool, bool> {
public:
   bool operator()(const bool& b0, const bool& b1) const {
      return (b0 != b1);
   }
};

template <>
class distance_functor<char, char>
 : public distanceable_functor<const char, char> {
public:
   char operator()(const char& c0, const char& c1) const {
      return (c0 < c1) ? (c1 - c0) : (c0 - c1);
   }
};

template <>
class distance_functor<unsigned char, unsigned char>
 : public distanceable_functor<const unsigned char, unsigned char> {
public:
   unsigned char operator()(
      const unsigned char& c0,
      const unsigned char& c1) const 
   {
      return (c0 < c1) ? (c1 - c0) : (c0 - c1);
   }
};

template <>
class distance_functor<short, short>
 : public distanceable_functor<const short, short> {
public:
   short operator()(const short& i0, const short& i1) const {
      return (i0 < i1) ? (i1 - i0) : (i0 - i1);
   }
};

template <>
class distance_functor<unsigned short, unsigned short>
 : public distanceable_functor<const unsigned short, unsigned short> {
public:
   unsigned short operator()(
      const unsigned short& i0,
      const unsigned short& i1) const
   {
      return (i0 < i1) ? (i1 - i0) : (i0 - i1);
   }
};

template <>
class distance_functor<int, int>
 : public distanceable_functor<const int, int> {
public:
   int operator()(const int& i0, const int& i1) const {
      return (i0 < i1) ? (i1 - i0) : (i0 - i1);
   }
};

template <>
class distance_functor<unsigned int, unsigned int>
 : public distanceable_functor<const unsigned int, unsigned int> {
public:
   unsigned int operator()(
      const unsigned int& i0,
      const unsigned int& i1) const
   {
      return (i0 < i1) ? (i1 - i0) : (i0 - i1);
   }
};

template <>
class distance_functor<long, long>
 : public distanceable_functor<const long, long> {
public:
   long operator()(const long& i0, const long& i1) const {
      return (i0 < i1) ? (i1 - i0) : (i0 - i1);
   }
};

template <>
class distance_functor<unsigned long, unsigned long>
 : public distanceable_functor<const unsigned long, unsigned long> {
public:
   unsigned long operator()(
      const unsigned long& i0,
      const unsigned long& i1) const 
   {
      return (i0 < i1) ? (i1 - i0) : (i0 - i1);
   }
};

template <>
class distance_functor<long long, long long>
 : public distanceable_functor<const long long, long long> {
public:
   long long operator()(const long long& i0, const long long& i1) const {
      return (i0 < i1) ? (i1 - i0) : (i0 - i1);
   }
};

template <>
class distance_functor<unsigned long long, unsigned long long>
 : public distanceable_functor<const unsigned long long, unsigned long long> {
public:
   unsigned long long operator()(
      const unsigned long long& i0,
      const unsigned long long& i1) const
   {
      return (i0 < i1) ? (i1 - i0) : (i0 - i1);
   }
};

template <>
class distance_functor<float, float>
 : public distanceable_functor<const float, float> {
public:
   float operator()(const float& f0, const float& f1) const {
      return (f0 < f1) ? (f1 - f0) : (f0 - f1);
   }
};

template <>
class distance_functor<double, double>
 : public distanceable_functor<const double, double> {
public:
   double operator()(const double& d0, const double& d1) const {
      return (d0 < d1) ? (d1 - d0) : (d0 - d1);
   }
};

template <>
class distance_functor<long double, long double>
 : public distanceable_functor<const long double, long double> {
public:
   long double operator()(const long double& d0, const long double& d1) const {
      return (d0 < d1) ? (d1 - d0) : (d0 - d1);
   }
};

} /* namespace functors */

#endif
