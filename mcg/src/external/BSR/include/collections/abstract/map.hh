/*
 * Abstract map.
 *
 * Maps associate each element of their domain with exactly one element of 
 * their range.  The domain of a map is a unique collection of elements.  The 
 * range may contain multiple copies of the same element, each with a different
 * preimage in the domain.  Equivalently, the mapping stored in a map is onto.
 *
 * Maps are not derived from the abstract collection class because elements 
 * must be added to a map as a pair (t -> u mapping).  A collection view of 
 * the domain and range of a map is accessible via the domain() and range() 
 * methods.
 */
#ifndef COLLECTIONS__ABSTRACT__MAP_HH
#define COLLECTIONS__ABSTRACT__MAP_HH

#include "collections/abstract/collection.hh"
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
 * Abstract base class for maps.
 * T is the type of elements in the domain.
 * U is the type of elements in the range.
 */
template <typename T, typename U>
class map {
public:
   /*
    * Destructor.
    */
   virtual ~map() = 0;
   
   /*
    * Add the mapping t -> u, overwriting any existing mapping for t.
    * Return a reference to the map.
    */
   virtual map<T,U>& add(T& /* t */, U& /* u */) = 0;

   /*
    * Add the mapping t -> u for corresponding elements of the collections.
    * Overwriting previous mappings if they exist.
    * The collections must have the same size.
    * Return a reference to the map.
    */
   virtual map<T,U>& add(const collection<T>&, const collection<U>&) = 0;

   /*
    * Add all mappings in the given map to the current map.
    * Overwriting previous mappings if they exist.
    * Return a reference to the map.
    */
   virtual map<T,U>& add(const map<T,U>&) = 0;
   
   /*
    * Remove element(s) and their corresponding images from the map.
    * Return a reference to the map.
    */
   virtual map<T,U>& remove(const T&) = 0;
   virtual map<T,U>& remove(const collection<T>&) = 0;

   /*
    * Clear the map, resetting it to the empty map.
    * Return a reference to the map.
    */
   virtual map<T,U>& clear() = 0;

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
    * This method allows one to atomically capture the contents of a map if 
    * there is a potential for the map to change between calls of the above 
    * domain() and range() methods.
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
template <typename T, typename U>
map<T,U>::~map() { }

} /* namespace abstract */
} /* namespace collections */

#endif
