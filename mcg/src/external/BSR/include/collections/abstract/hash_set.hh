/*
 * Abstract hash set.
 *
 * A hash set is a set which stores items in bins according to their hash.
 *
 * Hash sets do not allow duplicate items.
 * Adding an item already in the set has no effect.
 */
#ifndef COLLECTIONS__ABSTRACT__HASH_SET_HH
#define COLLECTIONS__ABSTRACT__HASH_SET_HH

#include "collections/abstract/collection.hh"
#include "collections/abstract/set.hh"
#include "collections/list.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"

namespace collections {
namespace abstract {
/*
 * Imports.
 */
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

/*
 * Abstract base class for hash sets.
 * T is the element type.
 * H is the hash type.
 */
template <typename T, typename H>
class hash_set : virtual public set<T> {
public:
   /*
    * Destructor.
    */
   virtual ~hash_set() = 0;

   /*
    * Add element(s) to the hash set.
    * Return a reference to the set.
    */
   virtual hash_set<T,H>& add(T&) = 0;
   virtual hash_set<T,H>& add(const collection<T>&) = 0;
   
   /*
    * Remove element(s) from the hash set.
    * Return a reference to the set.
    */
   virtual hash_set<T,H>& remove(const T&) = 0;
   virtual hash_set<T,H>& remove(const collection<T>&) = 0;

   /*
    * Remove all element(s) with the given hash from the set.
    * Return a reference to the set.
    */
   virtual hash_set<T,H>& remove_hash(const H&) = 0;
    
   /*
    * Remove all element(s) from the hash set.
    * Return a reference to the set.
    */
   virtual hash_set<T,H>& clear() = 0;
   
   /*
    * Search.
    * Return the element in the set matching the given element.
    * Throw an exception (ex_not_found) when attempting to find an element not
    * contained in the set.
    */
   virtual bool contains(const T&) const = 0;
   virtual T& find(const T&) const = 0;
  
   /*
    * Search.
    * Find an element with the given hash.
    * Throw an exception (ex_not_found) when attempting to find an element not
    * contained in the set.
    */
   virtual bool contains_hash(const H&) const = 0;
   virtual T& find_hash(const H&) const = 0;
   
   /*
    * Search.
    * Find all elements in the set with the same hash as the given element.
    */
   virtual collections::list<T> find_bin(const T&) const = 0;

   /*
    * Search.
    * Find all elements in the set with the given hash.
    */
   virtual collections::list<T> find_bin_hash(const H&) const = 0;

   /*
    * Size.
    */
   virtual unsigned long size() const = 0;

   /*
    * Return iterator over elements.
    */
   virtual auto_ptr< iterator<T> > iter_create() const = 0;
};

/*
 * Pure virtual destructor.
 */
template <typename T, typename H>
hash_set<T,H>::~hash_set() { }

} /* namespace abstract */
} /* namespace collections */

#endif
