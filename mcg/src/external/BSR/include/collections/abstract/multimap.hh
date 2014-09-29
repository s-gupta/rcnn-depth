/*
 * Abstract multimap.
 *
 * Multimaps associate each element of their domain with one or more elements 
 * of their range.  Both the domain and range may contain multiple copies of 
 * the same element.
 *
 * Multimaps are not derived from the abstract collection class because 
 * elements must be added to a multimap as a pair (t -> u mapping).  A 
 * collection view of the domain and range of a multimap is accessible via 
 * the domain() and range() methods.
 */
#ifndef COLLECTIONS__ABSTRACT__MULTIMAP_HH
#define COLLECTIONS__ABSTRACT__MULTIMAP_HH

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
 * Abstract base class for multimaps.
 * T is the type of elements in the domain.
 * U is the type of elements in the range.
 */
template <typename T, typename U>
class multimap {
public:
   /*
    * Destructor.
    */
   virtual ~multimap() = 0;

   /*
    * Add the mapping t -> u.
    * Return a reference to the multimap.
    */
   virtual multimap<T,U>& add(T& /* t */, U& /* u */) = 0;

   /*
    * Add the mapping t -> u for corresponding elements of the collections.
    * The collections must have the same size.
    * Return a reference to the multimap.
    */
   virtual multimap<T,U>& add(const collection<T>&, const collection<U>&) = 0;

   /*
    * Add all mappings in the given map/multimap to the current multimap.
    * Return a reference to the multimap.
    */
   virtual multimap<T,U>& add(const map<T,U>&) = 0;
   virtual multimap<T,U>& add(const multimap<T,U>&) = 0;

   /*
    * Remove element(s) and their corresponding images from the multimap.
    * Return a reference to the multimap.
    */
   virtual multimap<T,U>& remove(const T&) = 0;
   virtual multimap<T,U>& remove(const collection<T>&) = 0;

   /*
    * Remove all instances of the element(s) from the multimap.
    * Return a reference to the multimap.
    */
   virtual multimap<T,U>& remove_all(const T&) = 0;
   virtual multimap<T,U>& remove_all(const collection<T>&) = 0;

   /*
    * Clear the multimap, resetting it to the empty multimap.
    * Return a reference to the multimap.
    */
   virtual multimap<T,U>& clear() = 0;

   /*
    * Search.
    * Return an element in the domain matching the given element.
    * Throw an exception (ex_not_found) when attempting to find an element not 
    * contained in the map.
    */
   virtual bool contains(const T&) const = 0;
   virtual T& find(const T&) const = 0;

   /*
    * Search.
    * Return a list of all element(s) in the domain matching the given element.
    */
   virtual collections::list<T> find_all(const T&) const = 0;

   /*
    * Search.
    * Return an image of the given domain element under the multimap.
    * Throw an exception (ex_not_found) when attempting to find the image of
    * an element not contained in the multimap.
    */
   virtual U& find_image(const T&) const = 0;

   /*
    * Search.
    * Return a list of all element(s) in the range which are images of the 
    * given domain element.
    */
   virtual collections::list<U> find_images(const T&) const = 0;

   /*
    * Size.
    * Return the number of mappings.
    */
   virtual unsigned long size() const = 0;
   
   /*
    * Return a list of elements in the domain of the multimap and a list of 
    * their corresponding images in the range.  All mappings t -> u in the 
    * multimap appear as a pair of elements (t, u) at corresponding positions
    * in the lists.  Hence, both the domain and range list may have repeated 
    * elements.
    */
   virtual collections::list<T> domain() const = 0;
   virtual collections::list<U> range() const = 0;
   
   /*
    * Return lists of elements in the domain and range.
    * This method allows one to atomically capture the contents of a multimap
    * if there is a potential for the multimap to change between calls of the
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
template <typename T, typename U>
multimap<T,U>::~multimap() { }

} /* namespace abstract */
} /* namespace collections */

#endif
