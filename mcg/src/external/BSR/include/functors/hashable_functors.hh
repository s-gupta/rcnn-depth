/*
 * Hashable functors.
 */
#ifndef FUNCTORS__HASHABLE_FUNCTORS_HH
#define FUNCTORS__HASHABLE_FUNCTORS_HH

#include "functors/mappable_functors.hh"
#include "interfaces/hashable.hh"

namespace functors {
/*
 * Imports.
 */
using interfaces::hashable;

/*
 * Functor for hashing items.
 * Call the underlying hash method.
 */
template <typename T, typename H>
class hash_functor : public mappable_functor<const T, H> {
public:
   using mappable_functor<const T, H>::operator();
   H& operator()(const T&) const;
};

template <typename T, typename H>
H& hash_functor<T,H>::operator()(const T& t) const {
   return hashable<T,H>::hash(t);
}

/*
 * Functor for hashing items with const qualifier.
 * The same functor is used as for types without the qualifier.
 */
template <typename T, typename H>
class hash_functor<const T, H> : public hash_functor<T,H> { };

/*
 * Globally accessible set of default hash functors.
 */
template <typename T, typename H>
class hash_functors {
public:
   static const hash_functor<T,H>& f_hash();
};

template <typename T, typename H>
const hash_functor<T,H>& hash_functors<T,H>::f_hash() {
   static const hash_functor<T,H>* f = new hash_functor<T,H>();
   return *f;
}

} /* namespace functors */

#endif
