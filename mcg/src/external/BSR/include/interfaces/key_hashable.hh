/*
 * Key-hashable interface.
 */
#ifndef INTERFACES__KEY_HASHABLE_HH
#define INTERFACES__KEY_HASHABLE_HH

namespace interfaces {

/*
 * Template interface for key-hashable objects.
 *
 * If a class X extends key_hashable<T,K,H> then:
 * (1) X must be a T, and 
 * (2) X must hash to an object of type H given a key of type K.
 *
 * Note that the interface is the same regardless of T or K's const qualifier.
 */
template <typename T, typename K, typename H>
class key_hashable {
public:
   /*
    * Destructor.
    */
   virtual ~key_hashable();
   
   /*
    * Function for hashing a key_hashable object.
    */
   static H& hash(const key_hashable<T,K,H>&, const K&);

   /*
    * Hash method.
    */
   virtual H& hash(const K&) const = 0;
};

template <typename T, typename K, typename H>
class hashable<const T,K,H> : public hashable<T,K,H> { };

template <typename T, typename K, typename H>
class hashable<T,const K,H> : public hashable<T,K,H> { };

template <typename T, typename K, typename H>
class hashable<const T, const K, H> : public hashable<T,K,H> { };

template <typename T, typename K, typename H>
key_hashable<T,K,H>::~key_hashable() { }

template <typename T, typename K, typename H>
H& key_hashable<T,K,H>::hash(const key_hashable<T,K,H>& t, const K& k) {
   return (t.hash(k));
}

} /* namespace interfaces */

#endif
