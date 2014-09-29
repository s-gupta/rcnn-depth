/*
 * Abstract hash map.
 *
 * A hash map is a map which stores t -> u mappings in bins according to the 
 * hash of t.
 *
 * Hash maps do not allow duplicate items in the domain.
 * Note that maps are not derived from the abstract collection class.
 */
#ifndef COLLECTIONS__ABSTRACT__HASH_MAP_HH
#define COLLECTIONS__ABSTRACT__HASH_MAP_HH

#include "collections/abstract/collection.hh"
#include "collections/abstract/map.hh"
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
 * Abstract base class for hash maps.
 * T is the domain element type.
 * U is the range element type.
 * H is the hash type for elements of the domain.
 */
template <typename T, typename U, typename H>
class hash_map : virtual public map<T,U> {
public:
   /*
    * Destructor.
    */
   virtual ~hash_map() = 0;
   
   /*
    * Add the mapping t -> u, overwriting any existing mapping for t.
    * Return a reference to the map.
    */
   virtual hash_map<T,U,H>& add(T& /* t */, U& /* u */) = 0;

   /*
    * Add the mapping t -> u for corresponding elements of the collections.
    * Overwriting previous mappings if they exist.
    * The collections must have the same size.
    * Return a reference to the map.
    */
   virtual hash_map<T,U,H>& add(
      const collection<T>&,
      const collection<U>&
   ) = 0;

   /*
    * Add all mappings in the given map to the current map.
    * Overwriting previous mappings if they exist.
    * Return a reference to the map.
    */
   virtual hash_map<T,U,H>& add(const map<T,U>&) = 0;
   
   /*
    * Remove element(s) and their corresponding images from the map.
    * Return a reference to the map.
    */
   virtual hash_map<T,U,H>& remove(const T&) = 0;
   virtual hash_map<T,U,H>& remove(const collection<T>&) = 0;
   
   /*
    * Remove all element(s) with the given hash and their images from the map.
    * Return a reference to the map.
    */
   virtual hash_map<T,U,H>& remove_hash(const H&) = 0;

   /*
    * Clear the map, resetting it to the empty map.
    * Return a reference to the map.
    */
   virtual hash_map<T,U,H>& clear() = 0;

   /*
    * Search.
    * Return the element in the domain matching the given element.
    * Throw an exception (ex_not_found) when attempting to find an element not
    * contained in the map.
    */
   virtual bool contains(const T&) const = 0;
   virtual T& find(const T&) const = 0;
   
   /*
    * Search.
    * Find an element in the domain with the given hash.
    * Throw an exception (ex_not_found) when attempting to find an element not
    * contained in the set.
    */
   virtual bool contains_hash(const H&) const = 0;
   virtual T& find_hash(const H&) const = 0;

   /*
    * Search.
    * Find all elements in the domain of the map with the same hash as the 
    * given element.
    */
   virtual collections::list<T> find_bin(const T&) const = 0;

   /*
    * Search.
    * Find all elements in the domain of the map with the given hash.
    */
   virtual collections::list<T> find_bin_hash(const H&) const = 0;

   /*
    * Search.
    * Return the image of the given domain element under the map.
    * Throw an exception (ex_not_found) when attempting to find the image of
    * an element not contained in the map.
    */
   virtual U& find_image(const T&) const = 0;

   /*
    * Size.
    * Return the number of elements in the domain.
    */
   virtual unsigned long size() const = 0;
   
   /*
    * Return a list of elements in the domain of the map and a list of their 
    * corresponding images in the range.  Hence, range() does not necessarily 
    * return a unique list of elements as two or more elements of the domain 
    * might have the same image.
    */
   virtual collections::list<T> domain() const = 0;
   virtual collections::list<U> range() const = 0;

   /*
    * Return lists of elements in the domain and range.
    * This method allows one to atomically capture the contents of a hash map
    * if there is a potential for the hash map to change between calls of the
    * above domain() and range() methods.
    */
   virtual void contents(
      auto_ptr< collections::list<T> >&,  /* domain */
      auto_ptr< collections::list<U> >&   /* range */
   ) const = 0;

   /*
    * Return iterator over elements in the domain.
    */
   virtual auto_ptr< iterator<T> > iter_create() const = 0;
};

/*
 * Pure virtual destructor.
 */
template <typename T, typename U, typename H>
hash_map<T,U,H>::~hash_map() { }

} /* namespace abstract */
} /* namespace collections */

#endif
