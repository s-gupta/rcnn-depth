/*
 * Multisets (thread-safe).
 *
 * Multisets allow duplicate items.
 */
#ifndef COLLECTIONS__MULTISET_HH
#define COLLECTIONS__MULTISET_HH

#include "collections/abstract/collection.hh"
#include "collections/abstract/multiset.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "collections/set.hh"
#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "functors/comparable_functors.hh"
#include "lang/exceptions/ex_not_found.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"

namespace collections {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::pointers::auto_collection;
using concurrent::threads::synchronization::locks::auto_read_lock;
using concurrent::threads::synchronization::locks::auto_write_lock;
using concurrent::threads::synchronization::synchronizables::unsynchronized;
using functors::comparable_functor;
using functors::compare_functors;
using lang::exceptions::ex_not_found;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

/*
 * Declare classes for iterators over multisets.
 */
template <typename T, typename Syn>
class multiset_iterator;

template <typename T, typename Syn>
class multiset_iterator_reverse;

/*
 * Multisets.
 */
template <typename T, typename Syn = unsynchronized>
class multiset : public abstract::multiset<T>,
                 protected Syn {
public:
   /*
    * Friend classes.
    */
   friend class multiset_iterator<T,Syn>;
   friend class multiset_iterator_reverse<T,Syn>;

   /*
    * Define the iterator types.
    */
   typedef multiset_iterator<T,Syn> iterator_t;
   typedef multiset_iterator_reverse<T,Syn> iterator_reverse_t;

   /*
    * Constructor.
    * Return an empty multiset.
    */
   multiset(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Constructor.
    * Create a multiset from a collection.
    */
   explicit multiset(
      const collection<T>&, 
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );
   
   /*
    * Copy constructor.
    */
   multiset(const multiset<T,Syn>&);

   /*
    * Destructor.
    */
   virtual ~multiset();

   /*
    * Add element(s) to the multiset.
    * Return a reference to the multiset.
    */
   multiset<T,Syn>& add(T&);
   multiset<T,Syn>& add(const collection<T>&);

   /*
    * Remove element(s) from the multiset.
    * Return a reference to the multiset.
    */
   multiset<T,Syn>& remove(const T&);
   multiset<T,Syn>& remove(const collection<T>&);

   /*
    * Remove all instances of the element(s) from the multiset.
    * Return a reference to the multiset.
    */
   multiset<T,Syn>& remove_all(const T&);
   multiset<T,Syn>& remove_all(const collection<T>&);

   /*
    * Remove all element(s) from the multiset.
    * Return a reference to the multiset.
    */
   multiset<T,Syn>& clear();
    
   /*
    * Search.
    * Throw an exception (ex_not_found) when attempting to find an element not
    * contained in the multiset.
    */
   bool contains(const T&) const;
   T& find(const T&) const;

   /*
    * Search.
    * Return a list of all element(s) matching the given element.
    */
   list<T> find_all(const T&) const;

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
   /************************************************************************
    * Multiset data structures.
    ************************************************************************/

   /*
    * Bin containing equivalent (as defined by comparison functor) set items.
    */
   class equiv_bin {
   public:
      /*
       * Constructor.
       */
      explicit equiv_bin(T& t)
       : items()
      { items.add(t); }

      /*
       * Copy constructor.
       */
      equiv_bin(const equiv_bin& bin)
       : items(bin.items)
      { }

      /*
       * Destructor.
       */
      ~equiv_bin() { /* do nothing */ }

      /*
       * Equivalent items.
       */
      list<T> items;
   };

   /*
    * Comparison functor on bins.
    * Call the given compare functor on representative elements.
    */
   class equiv_bin_compare_functor : public comparable_functor<equiv_bin> {
   public:
      /*
       * Constructor.
       */
      explicit equiv_bin_compare_functor(const comparable_functor<T>& f)
       : _f(f) { }
       
      /*
       * Copy constructor.
       */
      explicit equiv_bin_compare_functor(const equiv_bin_compare_functor& f)
       : _f(f._f) { }

      /*
       * Comparison function.
       */
      int operator()(const equiv_bin& bin0, const equiv_bin& bin1) const {
         return _f(bin0.items.head(), bin1.items.head());
      }

   protected:
      const comparable_functor<T>& _f;
   };

   /*
    * Multiset data.
    */
   const equiv_bin_compare_functor              _f_compare; /* comparison function */
   auto_collection< equiv_bin, set<equiv_bin> > _set;       /* set of equivalence classes */
   unsigned long                                _size;      /* multiset size */

   /************************************************************************
    * Multiset helper functions.
    ************************************************************************/
    
   /*
    * Add an element to the multiset.
    */
   void add_item(T&);

   /*
    * Remove an element from the multiset.
    */
   void remove_item(const T&);

   /*
    * Remove all instances of an element from the multiset.
    */
   void remove_item_all(const T&);
};

/*
 * Multiset iterator.
 */
template <typename T, typename Syn = unsynchronized>
class multiset_iterator : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit multiset_iterator(const multiset<T,Syn>&);

   /*
    * Copy constructor.
    */
   multiset_iterator(const multiset_iterator<T,Syn>&);

   /*
    * Destructor.
    */
   virtual ~multiset_iterator();

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
   /*
    * Iterator data.
    */
   const multiset<T,Syn>&                                        _s;         /* multiset being iterated over */
   auto_read_lock<const Syn>                                     _rlock;     /* read lock on set */
   auto_ptr< set_iterator<typename multiset<T,Syn>::equiv_bin> > _bin_iter;  /* iterator over bins */
   auto_ptr< list_iterator<T> >                                  _item_iter; /* iterator over items in current bin */

   /*
    * Update bin and item iterators for next item.
    */
   void next_item();
};

/*
 * Multiset reverse iterator.
 */
template <typename T, typename Syn = unsynchronized>
class multiset_iterator_reverse : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit multiset_iterator_reverse(const multiset<T,Syn>&);

   /*
    * Copy constructor.
    */
   multiset_iterator_reverse(const multiset_iterator_reverse<T,Syn>&);

   /*
    * Destructor.
    */
   virtual ~multiset_iterator_reverse();

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
   /*
    * Iterator data.
    */
   const multiset<T,Syn>&                                                _s;         /* multiset being iterated over */
   auto_read_lock<const Syn>                                             _rlock;     /* read lock on set */
   auto_ptr< set_iterator_reverse<typename multiset<T,Syn>::equiv_bin> > _bin_iter;  /* iterator over bins */
   auto_ptr< list_iterator_reverse<T> >                                  _item_iter; /* iterator over items in current bin */

   /*
    * Update bin and item iterators for next item.
    */
   void next_item();
};

/***************************************************************************
 * Multiset iterator implementation.
 ***************************************************************************/

/*
 * Constructors.
 * Initialize position.
 */
template <typename T, typename Syn>
multiset_iterator<T,Syn>::multiset_iterator(
   const multiset<T,Syn>& s)
 : _s(s),
   _rlock(_s),
   _bin_iter(
      new set_iterator<typename multiset<T,Syn>::equiv_bin>(*(_s._set))
   ),
   _item_iter(
      (_bin_iter->has_next()) ?
         (new list_iterator<T>(_bin_iter->next().items)) : NULL
   )
{
   this->next_item();
}

template <typename T, typename Syn>
multiset_iterator_reverse<T,Syn>::multiset_iterator_reverse(
   const multiset<T,Syn>& s)
 : _s(s),
   _rlock(_s),
   _bin_iter(
      new set_iterator_reverse<typename multiset<T,Syn>::equiv_bin>(*(_s._set))
   ),
   _item_iter(
      (_bin_iter->has_next()) ?
         (new list_iterator_reverse<T>(_bin_iter->next().items)) : NULL
   )
{
   this->next_item();
}

/*
 * Copy constructors.
 */
template <typename T, typename Syn>
multiset_iterator<T,Syn>::multiset_iterator(
   const multiset_iterator<T,Syn>& i)
 : _s(i._s),
   _rlock(_s),
   _bin_iter(
      new set_iterator<typename multiset<T,Syn>::equiv_bin>(*(i._bin_iter))
   ),
   _item_iter(
      (i._item_iter.get() == NULL) ?
         NULL : (new list_iterator<T>(*(i._item_iter)))
   )
{ }

template <typename T, typename Syn>
multiset_iterator_reverse<T,Syn>::multiset_iterator_reverse(
   const multiset_iterator_reverse<T,Syn>& i)
 : _s(i._s),
   _rlock(_s),
   _bin_iter(
      new set_iterator_reverse<typename multiset<T,Syn>::equiv_bin>(
         *(i._bin_iter)
      )
   ),
   _item_iter(
      (i._item_iter.get() == NULL) ?
         NULL : (new list_iterator_reverse<T>(*(i._item_iter)))
   )
{ }

/*
 * Destructors.
 * Do nothing as the multiset is automatically unlocked upon destruction
 * of the read lock.
 */
template <typename T, typename Syn>
multiset_iterator<T,Syn>::~multiset_iterator() {
   /* do nothing */
}

template <typename T, typename Syn>
multiset_iterator_reverse<T,Syn>::~multiset_iterator_reverse() {
   /* do nothing */
}

/*
 * Check if there is another item available.
 */
template <typename T, typename Syn>
bool multiset_iterator<T,Syn>::has_next() const {
   return (_item_iter.get() != NULL);
}
   
template <typename T, typename Syn>
bool multiset_iterator_reverse<T,Syn>::has_next() const {
   return (_item_iter.get() != NULL);
}

/*
 * Return the next item.
 * Throw an exception (ex_not_found) if there are no more items.
 */
template <typename T, typename Syn>
T& multiset_iterator<T,Syn>::next() {
   if (_item_iter.get() != NULL) {
      T& t = _item_iter->next();
      this->next_item();
      return t;
   } else {
      throw ex_not_found(
         "iterator attempted to read past end of multiset"
      );
   }
}

template <typename T, typename Syn>
T& multiset_iterator_reverse<T,Syn>::next() {
   if (_item_iter.get() != NULL) {
      T& t = _item_iter->next();
      this->next_item();
      return t;
   } else {
      throw ex_not_found(
         "reverse iterator attempted to read past start of multiset"
      );
   }
}

/*
 * Update bin and item iterators for next item.
 */
template <typename T, typename Syn>
void multiset_iterator<T,Syn>::next_item() {
   if (_item_iter.get() != NULL) {
      /* check if all items in current bin visited */
      while (!(_item_iter->has_next())) {
         /* move to next bin */
         if (_bin_iter->has_next()) {
            _item_iter.reset(
               new list_iterator<T>(_bin_iter->next().items)
            );
         } else {
            _item_iter.reset(NULL);
            break;
         }
      }
   }
}

template <typename T, typename Syn>
void multiset_iterator_reverse<T,Syn>::next_item() {
   if (_item_iter.get() != NULL) {
      /* check if all items in current bin visited */
      while (!(_item_iter->has_next())) {
         /* move to next bin */
         if (_bin_iter->has_next()) {
            _item_iter.reset(
               new list_iterator_reverse<T>(_bin_iter->next().items)
            );
         } else {
            _item_iter.reset(NULL);
            break;
         }
      }
   }
}

/***************************************************************************
 * Multiset implementation.
 ***************************************************************************/
 
/*
 * Constructor.
 * Return an empty multiset.
 */
template <typename T, typename Syn>
multiset<T,Syn>::multiset(const comparable_functor<T>& f)
 : Syn(),
   _f_compare(f),
   _set(new set<equiv_bin>(_f_compare)),
   _size(0)
{ }

/*
 * Constructor.
 * Create a multiset from a collection.
 */
template <typename T, typename Syn>
multiset<T,Syn>::multiset(
   const collection<T>&         c, 
   const comparable_functor<T>& f)
 : Syn(),
   _f_compare(f),
   _set(new set<equiv_bin>(_f_compare)),
   _size(0)
{
   auto_ptr< iterator<T> > i = c.iter_create();
   while (i->has_next())
      this->add_item(i->next());
}

/*
 * Copy constructor.
 */
template <typename T, typename Syn>
multiset<T,Syn>::multiset(const multiset<T,Syn>& s)
 : Syn(),
   _f_compare(s.f_compare),
   _set(new set<equiv_bin>(_f_compare)),
   _size(0)
{
   typename multiset<T,Syn>::iterator_t i(s);
   while (i.has_next())
      this->add_item(i.next());
}

/*
 * Destructor.
 */
template <typename T, typename Syn>
multiset<T,Syn>::~multiset() {
   /* do nothing */
}

/*
 * Add an element to the multiset.
 */
template <typename T, typename Syn>
void multiset<T,Syn>::add_item(T& t) {
   auto_ptr<equiv_bin> bin(new equiv_bin(t));
   if (_set->contains(*bin)) {
      equiv_bin& bin_match = _set->find(*bin);
      bin_match.items.add(t);
   } else {
      _set->add(*bin);
      bin.release();
   }
   _size++;
}

/*
 * Add element(s) to the multiset.
 * Return a reference to the multiset.
 */
template <typename T, typename Syn>
multiset<T,Syn>& multiset<T,Syn>::add(T& t) {
   auto_write_lock<const Syn> wlock(*this);
   this->add_item(t);
   return *this;
}

template <typename T, typename Syn>
multiset<T,Syn>& multiset<T,Syn>::add(const collection<T>& c) {
   list<T> lst(c);
   auto_write_lock<const Syn> wlock(*this);
   for (typename list<T>::iterator_t i(lst); i.has_next(); )
      this->add_item(i.next());
   return *this;
}

/*
 * Remove an element from the multiset.
 */
template <typename T, typename Syn>
void multiset<T,Syn>::remove_item(const T& t) {
   equiv_bin bin(const_cast<T&>(t));   /* const_cast is safe here */
   if (_set->contains(bin)) {
      /* remove one instance of the element */
      equiv_bin& bin_match = _set->find(bin);
      bin_match.items.remove_head();
      _size--;
      /* remove bin if empty after element removal */
      if (bin_match.items.is_empty()) {
         _set->remove(bin_match);
         delete &bin_match;
      }
   }
}

/*
 * Remove element(s) from the multiset.
 * Return a reference to the multiset.
 */
template <typename T, typename Syn>
multiset<T,Syn>& multiset<T,Syn>::remove(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   this->remove_item(t);
   return *this;
}

template <typename T, typename Syn>
multiset<T,Syn>& multiset<T,Syn>::remove(const collection<T>& c) {
   list<T> lst(c);
   auto_write_lock<const Syn> wlock(*this);
   for (typename list<T>::iterator_t i(lst); i.has_next(); )
      this->remove_item(i.next());
   return *this;
}

/*
 * Remove all instances of an element from the multiset.
 */
template <typename T, typename Syn>
void multiset<T,Syn>::remove_item_all(const T& t) {
   equiv_bin bin(const_cast<T&>(t));   /* const_cast is safe here */
   if (_set->contains(bin)) {
      equiv_bin& bin_match = _set->find(bin);
      _set->remove(bin_match);
      _size -= bin_match.items.size();
      delete &bin_match;
   }
}

/*
 * Remove all instances of the element(s) from the multiset.
 * Return a reference to the multiset.
 */
template <typename T, typename Syn>
multiset<T,Syn>& multiset<T,Syn>::remove_all(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   this->remove_item_all(t);
   return *this;
}

template <typename T, typename Syn>
multiset<T,Syn>& multiset<T,Syn>::remove_all(const collection<T>& c) {
   list<T> lst(c);
   auto_write_lock<const Syn> wlock(*this);
   for (typename list<T>::iterator_t i(lst); i.has_next(); )
      this->remove_item_all(i.next());
   return *this;
}

/*
 * Remove all element(s) from the multiset.
 * Return a reference to the multiset.
 */
template <typename T, typename Syn>
multiset<T,Syn>& multiset<T,Syn>::clear() {
   auto_write_lock<const Syn> wlock(*this);
   _set.reset(new set<equiv_bin>(_f_compare));
   _size = 0;
   return *this;
}

/*
 * Check if an element is in the multiset.
 */
template <typename T, typename Syn>
bool multiset<T,Syn>::contains(const T& t) const {
   equiv_bin bin(const_cast<T&>(t));   /* const_cast is safe here */
   auto_read_lock<const Syn> rlock(*this);
   return _set->contains(bin);
}

/*
 * Find an element in the multiset that is equal to the given element.
 * Throw an exception (ex_not_found) if no such element is in the set.
 */
template <typename T, typename Syn>
T& multiset<T,Syn>::find(const T& t) const {
   equiv_bin bin(const_cast<T&>(t));   /* const_cast is safe here */
   auto_read_lock<const Syn> rlock(*this);
   equiv_bin& bin_match = _set->find(bin);
   return bin_match.items.head();
}

/*
 * Return a list of all element(s) matching the given element.
 */
template <typename T, typename Syn>
list<T> multiset<T,Syn>::find_all(const T& t) const {
   equiv_bin bin(const_cast<T&>(t));   /* const_cast is safe here */
   auto_read_lock<const Syn> rlock(*this);
   if (_set->contains(bin)) {
      equiv_bin& bin_match = _set->find(bin);
      return bin_match.items;
   } else {
      return list<T>();
   }
}

/*
 * Get number of elements (counting duplicates) in multiset.
 */
template <typename T, typename Syn>
unsigned long multiset<T,Syn>::size() const {
   auto_read_lock<const Syn> rlock(*this);
   return _size;
}

/*
 * Return iterator over elements.
 */
template <typename T, typename Syn>
auto_ptr< iterator<T> > multiset<T,Syn>::iter_create() const {
   return auto_ptr< iterator<T> >(new multiset_iterator<T,Syn>(*this));
}

/*
 * Return reverse iterator over elements.
 */
template <typename T, typename Syn>
auto_ptr< iterator<T> > multiset<T,Syn>::iter_reverse_create() const {
   return auto_ptr< iterator<T> >(new multiset_iterator_reverse<T,Syn>(*this));
}

} /* namespace collections */

#endif
