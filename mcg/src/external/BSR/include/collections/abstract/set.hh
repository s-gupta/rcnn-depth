/*
 * Abstract set.
 *
 * Sets do not allow duplicate items.
 * Adding an element to a set that already contains it has no effect.
 */
#ifndef COLLECTIONS__ABSTRACT__SET_HH
#define COLLECTIONS__ABSTRACT__SET_HH

#include "collections/abstract/collection.hh"
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
 * Abstract base class for sets.
 */
template <typename T>
class set : virtual public collection<T> {
public:
   /*
    * Destructor.
    */
   virtual ~set() = 0;

   /*
    * Add element(s) to the set.
    * Return a reference to the set.
    */
   virtual set<T>& add(T&) = 0;
   virtual set<T>& add(const collection<T>&) = 0;

   /*
    * Remove element(s) from the set.
    * Return a reference to the set.
    */
   virtual set<T>& remove(const T&) = 0;
   virtual set<T>& remove(const collection<T>&) = 0;

   /*
    * Remove all element(s) from the set.
    * Return a reference to the set.
    */
   virtual set<T>& clear() = 0;
    
   /*
    * Search.
    * Throw an exception (ex_not_found) when attempting to find an element not
    * contained in the set.
    */
   virtual bool contains(const T&) const = 0;
   virtual T& find(const T&) const = 0;

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
set<T>::~set() { }

} /* namespace abstract */
} /* namespace collections */

#endif
