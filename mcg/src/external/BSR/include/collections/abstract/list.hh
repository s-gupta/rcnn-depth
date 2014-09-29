/*
 * Abstract list.
 */
#ifndef COLLECTIONS__ABSTRACT__LIST_HH
#define COLLECTIONS__ABSTRACT__LIST_HH

#include "collections/abstract/collection.hh"
#include "functors/comparable_functors.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"

namespace collections {
namespace abstract {
/*
 * Imports.
 */
using functors::comparable_functor;
using functors::compare_functors;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

/*
 * Abstract base class for lists.
 */
template <typename T>
class list : virtual public collection<T> {
public:
   /*
    * Destructor.
    */
   virtual ~list() = 0;

   /*
    * Return head element(s).
    * Throw an exception (ex_not_found) if the list doesn't contain enough
    * elements.
    */
   virtual T& head() const = 0;
   virtual void head(unsigned long, collection<T>&) const = 0;

   /*
    * Return tail element(s). 
    * Throw an exception (ex_not_found) if the list doesn't contain enough
    * elements.
    */
   virtual T& tail() const = 0;
   virtual void tail(unsigned long, collection<T>&) const = 0;

   /*
    * Collection interface for adding element(s) to tail of list.
    * Append element(s) to the tail of the list.
    * Return a reference to the list.
    */
   virtual list<T>& add(T&);
   virtual list<T>& add(const collection<T>&);
    
   /*
    * Addition of element(s) to head of list.
    * Return a reference to the list.
    */
   virtual list<T>& prepend(T&) = 0;
   virtual list<T>& prepend(const collection<T>&) = 0;
   
   /*
    * Addition of element(s) to tail of list.
    * Return a reference to the list.
    */
   virtual list<T>& append(T&) = 0;
   virtual list<T>& append(const collection<T>&) = 0;

   /*
    * Remove and return head element(s).
    * Throw an exception (ex_not_found) if the list doesn't contain enough
    * elements.
    */
   virtual T& remove_head() = 0;
   virtual void remove_head(unsigned long, collection<T>&) = 0;

   /*
    * Remove and return tail element(s).
    * Throw an exception (ex_not_found) if the list doesn't contain enough
    * elements.
    */
   virtual T& remove_tail() = 0;
   virtual void remove_tail(unsigned long, collection<T>&) = 0;

   /*
    * Remove all element(s) from the list.
    * Return a reference to the list.
    */
   virtual list<T>& clear() = 0;

   /*
    * Reverse list.
    * Return a reference to the list.
    */
   virtual list<T>& reverse() = 0;
   
   /*
    * Size.
    */
   virtual unsigned long size() const = 0;

   /*
    * Return iterators over elements.
    */
   virtual auto_ptr< iterator<T> > iter_create() const = 0;
   virtual auto_ptr< iterator<T> > iter_reverse_create() const = 0;
 
   /*
    * Sort elements in ascending order using the given comparison functor.
    */
   virtual void sort(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   ) = 0;

   /*
    * Remove duplicate elements (as defined by the given comparison functor)
    * and place the remaining unique elements in sorted order.
    */
   virtual void unique(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   ) = 0;
};

/*
 * Pure virtual destructor.
 */
template <typename T>
list<T>::~list() { }

/*
 * Collection interface to add an element to the tail of the list.
 * Return a reference to the list.
 */
template <typename T>
list<T>& list<T>::add(T& t) {
   return this->append(t);
}

/*
 * Add all elements in a collection to the tail of the list.
 * Return a reference to the list.
 */
template <typename T>
list<T>& list<T>::add(const collection<T>& c) {
   return this->append(c);
}

} /* namespace abstract */
} /* namespace collections */

#endif
