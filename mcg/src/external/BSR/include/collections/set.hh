/*
 * Sets (thread-safe).
 *
 * Sets do not allow duplicate items.
 * Adding an element to a set that already contains it has no effect.
 *
 * This is a wrapper class for the default set implementation.
 */
#ifndef COLLECTIONS__SET_HH
#define COLLECTIONS__SET_HH

#include "collections/abstract/collection.hh"
#include "collections/abstract/set.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "functors/comparable_functors.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"

/*
 * Include the default set implementation.
 */
#include "collections/splay_set.hh"

namespace collections {
/*
 * Imports.
 */
using collections::abstract::collection;
using concurrent::threads::synchronization::synchronizables::unsynchronized;
using functors::comparable_functor;
using functors::compare_functors;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

/*
 * Declare classes for iterators over sets.
 */
template <typename T, typename Syn>
class set_iterator;

template <typename T, typename Syn>
class set_iterator_reverse;

/*
 * Sets.
 */
template <typename T, typename Syn = unsynchronized>
class set : public abstract::set<T> {
public:
   /*
    * Friend classes.
    */
   friend class set_iterator<T,Syn>;
   friend class set_iterator_reverse<T,Syn>;

   /*
    * Define the iterator types.
    */
   typedef set_iterator<T,Syn> iterator_t;
   typedef set_iterator_reverse<T,Syn> iterator_reverse_t;

   /*
    * Constructor.
    */
   set(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Constructor.
    * Create a set from a collection.
    */
   explicit set(
      const collection<T>&, 
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );
   
   /*
    * Copy constructor.
    */
   set(const set<T,Syn>&);
      
   /*
    * Destructor.
    */
   virtual ~set();

   /*
    * Add element(s) to the set.
    * Return a reference to the set.
    */
   set<T,Syn>& add(T&);
   set<T,Syn>& add(const collection<T>&);

   /*
    * Remove element(s) from the set.
    * Return a reference to the set.
    */
   set<T,Syn>& remove(const T&);
   set<T,Syn>& remove(const collection<T>&);

   /*
    * Remove all element(s) from the set.
    * Return a reference to the set.
    */
   set<T,Syn>& clear();

   /*
    * Search.
    * Throw an exception (ex_not_found) when attempting to find an element not 
    * contained in the set.
    */
   bool contains(const T&) const;
   T& find(const T&) const;

   /*
    * Size.
    */
   unsigned long size() const;

   /*
    * Return iterators over elements.
    */
   auto_ptr< iterator<T> > iter_create() const;
   auto_ptr< iterator<T> > iter_reverse_create() const;

protected:
   typedef splay_set<T,Syn> set_t;  /* type of underlying set */
   set_t _s;                        /* wrapped set object */
};

/*
 * Set iterator.
 */
template <typename T, typename Syn = unsynchronized>
class set_iterator : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit set_iterator(const set<T,Syn>&);

   /*
    * Copy constructor.
    */
   set_iterator(const set_iterator<T,Syn>&);

   /*
    * Destructor.
    */
   virtual ~set_iterator();

   /*
    * Check if there is another item available.
    */
   bool has_next() const;

   /*
    * Return the next item.
    * Throw an exception (ex_not_found) if there are no more items.
    */
   T& next();

protected:
   typename set<T,Syn>::set_t::iterator_t _s_iter;
};

/*
 * Set reverse iterator.
 */
template <typename T, typename Syn = unsynchronized>
class set_iterator_reverse : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit set_iterator_reverse(const set<T,Syn>&);

   /*
    * Copy constructor.
    */
   set_iterator_reverse(const set_iterator_reverse<T,Syn>&);

   /*
    * Destructor.
    */
   virtual ~set_iterator_reverse();

   /*
    * Check if there is another item available.
    */
   bool has_next() const;

   /*
    * Return the next item.
    * Throw an exception (ex_not_found) if there are no more items.
    */
   T& next();

protected:
   typename set<T,Syn>::set_t::iterator_reverse_t _s_iter_rev;
};

/***************************************************************************
 * Set iterator implementation.
 ***************************************************************************/

/*
 * Constructors.
 */
template <typename T, typename Syn>
set_iterator<T,Syn>::set_iterator(const set<T,Syn>& s)
 : _s_iter(s._s)
{ }

template <typename T, typename Syn>
set_iterator_reverse<T,Syn>::set_iterator_reverse(const set<T,Syn>& s)
 : _s_iter_rev(s._s)
{ }

/*
 * Copy constructors.
 */
template <typename T, typename Syn>
set_iterator<T,Syn>::set_iterator(
   const set_iterator<T,Syn>& i)
 : _s_iter(i._s_iter)
{ }

template <typename T, typename Syn>
set_iterator_reverse<T,Syn>::set_iterator_reverse(
   const set_iterator_reverse<T,Syn>& i)
 : _s_iter_rev(i._s_iter_rev)
{ }

/*
 * Destructors.
 * Do nothing as the underlying iterator's destructor handles destruction.
 */
template <typename T, typename Syn>
set_iterator<T,Syn>::~set_iterator() {
   /* do nothing */
}

template <typename T, typename Syn>
set_iterator_reverse<T,Syn>::~set_iterator_reverse() {
   /* do nothing */
}

/*
 * Check if there is another item available.
 */
template <typename T, typename Syn>
bool set_iterator<T,Syn>::has_next() const {
   return _s_iter.has_next();
}

template <typename T, typename Syn>
bool set_iterator_reverse<T,Syn>::has_next() const {
   return _s_iter_rev.has_next();
}

/*
 * Return the next item.
 * Throw an exception (ex_not_found) if there are no more items.
 */
template <typename T, typename Syn>
T& set_iterator<T,Syn>::next() {
   return _s_iter.next();
}

template <typename T, typename Syn>
T& set_iterator_reverse<T,Syn>::next() {
   return _s_iter_rev.next();
}

/***************************************************************************
 * Set implementation.
 ***************************************************************************/

/*
 * Constructor.
 */
template <typename T, typename Syn>
set<T,Syn>::set(const comparable_functor<T>& f)
 : _s(f)
{ }

/*
 * Constructor.
 * Create a set from a collection.
 */
template <typename T, typename Syn>
set<T,Syn>::set(
   const collection<T>& c,
   const comparable_functor<T>& f)
 : _s(c, f)
{ }

/*
 * Copy constructor.
 */
template <typename T, typename Syn>
set<T,Syn>::set(const set<T,Syn>& s)
 : _s(s._s)
{ }
   
/*
 * Destructor.
 * Do nothing as the member object's desctructor handles destruction.
 */
template <typename T, typename Syn>
set<T,Syn>::~set() {
   /* do nothing */
}

/*
 * Add element(s) to the set.
 * Return a reference to the set.
 */
template <typename T, typename Syn>
set<T,Syn>& set<T,Syn>::add(T& t) {
   _s.add(t);
   return *this;
}

template <typename T, typename Syn>
set<T,Syn>& set<T,Syn>::add(const collection<T>& c) {
   _s.add(c);
   return *this;
}

/*
 * Remove element(s) from the set.
 * Return a reference to the set.
 */
template <typename T, typename Syn>
set<T,Syn>& set<T,Syn>::remove(const T& t) {
   _s.remove(t);
   return *this;
}

template <typename T, typename Syn>
set<T,Syn>& set<T,Syn>::remove(const collection<T>& c) {
   _s.remove(c);
   return *this;
}

/*
 * Remove all element(s) from the set.
 * Return a reference to the set.
 */
template <typename T, typename Syn>
set<T,Syn>& set<T,Syn>::clear() {
   _s.clear();
   return *this;
}

/*
 * Search.
 */
template <typename T, typename Syn>
bool set<T,Syn>::contains(const T& t) const {
   return _s.contains(t);
}

template <typename T, typename Syn>
T& set<T,Syn>::find(const T& t) const {
   return _s.find(t);
}

/*
 * Size.
 */
template <typename T, typename Syn>
unsigned long set<T,Syn>::size() const {
   return _s.size();
}

/*
 * Return iterators over elements.
 */
template <typename T, typename Syn>
auto_ptr< iterator<T> > set<T,Syn>::iter_create() const {
   return auto_ptr< iterator<T> >(new set_iterator<T,Syn>(*this));
}
  
template <typename T, typename Syn>
auto_ptr< iterator<T> > set<T,Syn>::iter_reverse_create() const {
   return auto_ptr< iterator<T> >(new set_iterator_reverse<T,Syn>(*this));
}

} /* namespace collections */

#endif
