/*
 * Keyable interface.
 */
#ifndef INTERFACES__KEYABLE_HH
#define INTERFACES__KEYABLE_HH

#include "collections/abstract/collection.hh"

namespace interfaces {
/*
 * Imports.
 */
using collections::abstract::collection;

/*
 * Template interface for keyable objects.
 *
 * If a class X extends keyable<T,K> then:
 * (1) X must be a T,
 * (2) X must implement comparisons between itself and keys of type K, and 
 * (3) It must be possible to determine a key of type K that splits any 
 *     collection of X's into two disjoint subcollections that are 
 *     approximately equal in size.
 *
 * The goal of the keyable interface is to enable the use of generic search
 * algorithms.  For example, a balanced search tree whose internal nodes 
 * are of type K and whos leaf nodes are of type T can be built (using a 
 * generic algorithm) for any class that implements keyable<T,K>.
 */
template <typename T, typename K>
class keyable {
public:
   /*
    * Destructor.
    */
   virtual ~keyable();
   
   /*
    * Function for comparing a keyable object to a key.
    */
   static int compare(const keyable<T,K>&, const K&);
   
   /*
    * Virtual key comparison method.
    */
   virtual int compare_to(const K&) const = 0;

   /*
    * Compute a key that can be used to split the collection into two disjoint
    * subcollections.  The implementation should be in T::key_split(...).
    */
   static K key_split(const collection<T>&);
};

/*
 * Destructor.
 */
template <typename T, typename K>
keyable<T,K>::~keyable() { }

/*
 * Function for comparing a keyable object to a key.
 */
template <typename T, typename K>
int keyable<T,K>::compare(const keyable<T,K>& t, const K& k) {
   return t.compare_to(k);
}

/*
 * Compute a key that can be used to split the collection into two disjoint
 * subcollections.  The implementation should be in T::key_split(...).
 */
template <typename T, typename K>
K keyable<T,K>::key_split(const collection<T>& c) {
   return T::key_split(c);
}

} /* namespace interfaces */

#endif
