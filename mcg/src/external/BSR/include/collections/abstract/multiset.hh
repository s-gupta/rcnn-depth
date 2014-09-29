/*
 * Abstract multiset.
 *
 * Multisets allow duplicate items.  Adding an element to a multiset multiple
 * times results in a multiset which contains multiple copies of that element.
 */
#ifndef COLLECTIONS__ABSTRACT__MULTISET_HH
#define COLLECTIONS__ABSTRACT__MULTISET_HH

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
 * Abstract base class for multisets.
 */
template <typename T>
class multiset : virtual public collection<T> {
public:
   /*
    * Destructor.
    */
   virtual ~multiset() = 0;

   /*
    * Add element(s) to the multiset.
    * Return a reference to the multiset.
    */
   virtual multiset<T>& add(T&) = 0;
   virtual multiset<T>& add(const collection<T>&) = 0;

   /*
    * Remove element(s) from the multiset.
    * Return a reference to the multiset.
    */
   virtual multiset<T>& remove(const T&) = 0;
   virtual multiset<T>& remove(const collection<T>&) = 0;

   /*
    * Remove all instances of the element(s) from the multiset.
    * Return a reference to the multiset.
    */
   virtual multiset<T>& remove_all(const T&) = 0;
   virtual multiset<T>& remove_all(const collection<T>&) = 0;

   /*
    * Remove all element(s) from the multiset.
    * Return a reference to the multiset.
    */
   virtual multiset<T>& clear() = 0;
    
   /*
    * Search.
    * Throw an exception (ex_not_found) when attempting to find an element not
    * contained in the multiset.
    */
   virtual bool contains(const T&) const = 0;
   virtual T& find(const T&) const = 0;

   /*
    * Search.
    * Return a list of all element(s) matching the given element.
    */
   virtual collections::list<T> find_all(const T&) const = 0;

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
template <typename T>
multiset<T>::~multiset() { }

} /* namespace abstract */
} /* namespace collections */

#endif
