/*
 * Key-hashable functors.
 */
#ifndef FUNCTORS__KEY_HASHABLE_HH
#define FUNCTORS__KEY_HASHABLE_HH

#include "functors/mappable_functors.hh"
#include "interfaces/key_hashable.hh"

namespace functors {
/*
 * Imports.
 */
using interfaces::key_hashable;

/*
 * Functor for hashing items given a key.
 * Call the underlying hash method.
 */
template <typename T, typename K, typename H>
class key_hash_functor : public mappable_functor<const T, H> {
public:
   /*
    * Constructors.
    * Set the key to use for the hash.
    */
   key_hash_functor(const K&);
   key_hash_functor(const key_hash_functor<T,K,H>&);

   /*
    * Hash method.
    */
   using mappable_functor<const T, H>::operator();
   H& operator()(const T&) const;

private:
   const K& _key;
};

template <typename T, typename K, typename H>
key_hash_functor<T,K,H>::key_hash_functor(const K& k)
 : _key(k)
{ }

template <typename T, typename K, typename H>
key_hash_functor<T,K,H>::key_hash_functor(const key_hash_functor<T,K,H>& f)
 : _key(f._key)
{ }

template <typename T, typename K, typename H>
H& key_hash_functor<T,K,H>::operator()(const T& t) const {
   return key_hashable<T,K,H>::hash(t, _key);
}

/*
 * Functor for hashing items with const qualifier (given a key).
 * The same functor is used as for types without the qualifier.
 */
template <typename T, typename K, typename H>
class key_hash_functor<const T, K, H> : public key_hash_functor<T,K,H> { };

} /* namespace functors */

#endif
