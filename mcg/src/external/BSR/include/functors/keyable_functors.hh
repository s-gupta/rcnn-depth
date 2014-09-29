/*
 * Keyable functors.
 */
#ifndef FUNCTORS__KEYABLE_FUNCTORS_HH
#define FUNCTORS__KEYABLE_FUNCTORS_HH

#include "collections/abstract/collection.hh"
#include "interfaces/keyable.hh"

namespace functors {
/*
 * Imports.
 */
using collections::abstract::collection;
using interfaces::keyable;

/*
 * Abstract functor for compare operations.
 * Note that the functor is the same regardless of the type const qualifier.
 */
template <typename T, typename K>
class keyable_compare_functor {
public:
   virtual ~keyable_compare_functor() { }
   virtual int operator()(const T&, const K&) const = 0;
};

template <typename T, typename K>
class keyable_compare_functor<const T,K>
 : public keyable_compare_functor<T,K> { };

template <typename T, typename K>
class keyable_compare_functor<T,const K>
 : public keyable_compare_functor<T,K> { };

template <typename T, typename K>
class keyable_compare_functor<const T,const K>
 : public keyable_compare_functor<T,K> { };

/*
 * Abstract functor for split operations.
 */
template <typename T, typename K>
class keyable_split_functor {
public:
   virtual ~keyable_split_functor() { }
   virtual K operator()(const collection<T>&) const = 0;
};

/*
 * Functor for comparison of a keyable object to a key.
 * Call the underlying compare function.
 */
template <typename T, typename K>
class key_compare_functor : public keyable_compare_functor<const T,K> {
public:
   int operator()(const T& t, const K& k) const {
      return keyable<T,K>::compare(t, k);
   }
};

template <typename T, typename K>
class key_compare_functor<const T,K> : public key_compare_functor<T,K> { };

/*
 * Functor that reverses normal key comparison order
 */
template <typename T, typename K>
class key_compare_functor_reverse : public keyable_compare_functor<const T,K> {
public:
   key_compare_functor_reverse() : _f() { }
   int operator()(const T& t, const K& k) const {
      return (-(_f(t, k)));
   }
private:
   const key_compare_functor<T,K> _f;
};

template <typename T, typename K>
class key_compare_functor_reverse<const T,K>
 : public key_compare_functor_reverse<T,K> { };

/*
 * Functor for splitting of collections of keyable objects.
 * Call the underlying split function.
 */
template <typename T, typename K>
class key_split_functor : public keyable_split_functor<T,K> {
public:
   K operator()(const collection<T>& c) const {
      return keyable<T,K>::key_split(c);
   }
};

/*
 * Globally accessible set of default key comparison functors.
 */
template <typename T, typename K>
class key_compare_functors {
public:
   static const key_compare_functor<T,K>&         f_key_compare();
   static const key_compare_functor_reverse<T,K>& f_key_compare_reverse();
};

template <typename T, typename K>
const key_compare_functor<T,K>& key_compare_functors<T,K>::f_key_compare() {
   static const key_compare_functor<T,K>* f
      = new key_compare_functor<T,K>();
   return *f;
}

template <typename T, typename K>
const key_compare_functor_reverse<T,K>& key_compare_functors<T,K>::f_key_compare_reverse() {
   static const key_compare_functor_reverse<T,K>* f
      = new key_compare_functor_reverse<T,K>();
   return *f;
}

/*
 * Globally accessible set of default key split functors.
 */
template <typename T, typename K>
class key_split_functors {
public:
   static const key_split_functor<T,K>& f_key_split();
};

template <typename T, typename K>
const key_split_functor<T,K>& key_split_functors<T,K>::f_key_split() {
   static const key_split_functor<T,K>* f = new key_split_functor<T,K>();
   return *f;
}

} /* namespace functors */

#endif
