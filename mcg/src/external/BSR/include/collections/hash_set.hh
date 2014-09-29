/*
 * Hash sets (thread-safe).
 *
 * Hash sets do not allow duplicate items.
 * Adding an item already in the set has no effect.
 *
 * Hash sets can be automatically created for objects which implement the 
 * equalable and hashable interfaces and hash to objects which implement the 
 * comparable interface.  Hash sets for other types may be constructed by 
 * specifying appropriate functors at creation time.
 *
 * NOTE: Non-const objects returned by the hash functor are owned by the 
 * hash set and will be deleted upon destruction of the hash set.
 */
#ifndef COLLECTIONS__HASH_SET_HH
#define COLLECTIONS__HASH_SET_HH

#include "collections/abstract/collection.hh"
#include "collections/abstract/hash_set.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "collections/set.hh"
#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "functors/equalable_functors.hh"
#include "functors/comparable_functors.hh"
#include "functors/hashable_functors.hh"
#include "functors/iterable_functors.hh"
#include "functors/mappable_functors.hh"
#include "lang/exceptions/ex_not_found.hh"
#include "lang/iterators/iterator.hh"
#include "lang/null.hh"
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
using functors::equalable_functor;
using functors::equal_functors;
using functors::comparable_functor;
using functors::compare_functors;
using functors::hash_functors;
using functors::iter_functors;
using functors::iter_functor_delete;
using functors::mappable_functor;
using lang::exceptions::ex_not_found;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

/*
 * Declare classes for iterators over hash sets.
 */
template <typename T, typename H, typename Syn>
class hash_set_iterator;

template <typename T, typename H, typename Syn>
class hash_set_iterator_reverse;

/*
 * Hash sets.
 */
template <typename T, typename H, typename Syn = unsynchronized>
class hash_set : public abstract::hash_set<T,H>, 
                 protected Syn {
public:
   /*
    * Friend classes.
    */
   friend class hash_set_iterator<T,H,Syn>;
   friend class hash_set_iterator_reverse<T,H,Syn>;

   /*
    * Define the iterator types.
    */
   typedef hash_set_iterator<T,H,Syn> iterator_t;
   typedef hash_set_iterator_reverse<T,H,Syn> iterator_reverse_t;

   /*
    * Constructor.
    * Return an empty hash set.
    * Optionally specify custom functors for hash set operations.
    */
   explicit hash_set(
      const equalable_functor<T>&        = equal_functors<T>::f_equal(), 
      const comparable_functor<H>&       = compare_functors<H>::f_compare(),
      const mappable_functor<const T,H>& = (hash_functors<T,H>::f_hash())
   ); 

   /*
    * Constructor.
    * Create a hash set from the given collection.
    * Optionally specify custom functors for hash set operations.
    */
   explicit hash_set(
      const collection<T>&, 
      const equalable_functor<T>&        = equal_functors<T>::f_equal(), 
      const comparable_functor<H>&       = compare_functors<H>::f_compare(),
      const mappable_functor<const T,H>& = (hash_functors<T,H>::f_hash())
   );
      
   /*
    * Copy constructor.
    */
   hash_set(const hash_set<T,H,Syn>&);

   /*
    * Destructor.
    */
   virtual ~hash_set();

   /*
    * Add element(s) to the hash set.
    * Return a reference to the set.
    */
   hash_set<T,H,Syn>& add(T&);
   hash_set<T,H,Syn>& add(const collection<T>&);
   
   /*
    * Remove element(s) from the hash set.
    * Return a reference to the set.
    */
   hash_set<T,H,Syn>& remove(const T&);
   hash_set<T,H,Syn>& remove(const collection<T>&);
   
   /*
    * Remove all element(s) with the given hash from the set.
    * Return a reference to the set.
    */
   hash_set<T,H,Syn>& remove_hash(const H&);

   /*
    * Remove all element(s) from the hash set.
    * Return a reference to the set.
    */
   hash_set<T,H,Syn>& clear();
   
   /*
    * Search.
    * Return an element in the set matching the given element.
    * Throw an exception (ex_not_found) when attempting to find an element not
    * contained in the set.
    */
   bool contains(const T&) const;
   T& find(const T&) const;
   
   /*
    * Search.
    * Find an element with the given hash.
    * Throw an exception (ex_not_found) when attempting to find an element not
    * contained in the set.
    */
   bool contains_hash(const H&) const;
   T& find_hash(const H&) const;
   
   /*
    * Search.
    * Find all elements in the set with the same hash as the given element.
    */
   list<T> find_bin(const T&) const;

   /*
    * Search.
    * Find all elements in the set with the given hash.
    */
   list<T> find_bin_hash(const H&) const;

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
    * Hash set data structures.
    ************************************************************************/

   /*
    * Bin in the hash set.
    */
   class hash_bin {
   public:
      /*
       * Constructor.
       */
      explicit hash_bin(H& h)
       : hash_value(h),
         items()
      { }

      /*
       * Destructor.
       */
      virtual ~hash_bin() { /* do nothing */ }

      /*
       * Bin data.
       */
      H& hash_value;    /* hash value of items in the bin */
      list<T> items;    /* items sharing the hash value */
      
   protected:
      /*
       * Copy constructor.
       * WARNING: Does not create a deep copy of the hash value.
       */
      explicit hash_bin(const hash_bin& bin)
       : hash_value(bin.hash_value),
         items(bin.items)
      { }
   };
   
   /*
    * Auto hash bin.
    * An auto_hash_bin is a bin which owns its hash value and deletes the 
    * hash value (iff it is a non-const object) upon the bin's destruction.
    */
   class auto_hash_bin : public hash_bin {
   public:
      /*
       * Constructor.
       */
      explicit auto_hash_bin(H& h) : hash_bin(h) { }

      /*
       * Destructor.
       * Destroy the hash value (iff non-const object).
       */
      virtual ~auto_hash_bin() {
         const iter_functor_delete<H>& f_delete = iter_functors<H>::f_delete();
         f_delete(this->hash_value);
      }
      
   private:
      /*
       * Copy constructor.
       * WARNING: Should not be used - does not create a deep copy of the 
       * owned hash value.  It is declared private to prevent use.
       */
      explicit auto_hash_bin(const auto_hash_bin& bin) : hash_bin(bin) { }
   }; 
   
   /*
    * Comparison functor on hash bins.
    */
   class hash_bin_compare_functor : public comparable_functor<hash_bin> {
   public:
      /*
       * Constructor.
       */
      explicit hash_bin_compare_functor(const comparable_functor<H>& f)
       : _f(f) { }         

      /*
       * Copy constructor.
       */
      explicit hash_bin_compare_functor(const hash_bin_compare_functor& f)
       : _f(f._f) { }

      /*
       * Comparison function.
       */
      int operator()(const hash_bin& bin0, const hash_bin& bin1) const {
         return _f(bin0.hash_value, bin1.hash_value);
      }
      
   protected:
      const comparable_functor<H>& _f;
   };
  
   /*
    * Hash set data.
    */
   const equalable_functor<T>&                _f_equal;    /* equality function */
   const hash_bin_compare_functor             _f_compare;  /* comparison function */
   const mappable_functor<const T,H>&         _f_hash;     /* hash function */
   auto_collection< hash_bin, set<hash_bin> > _set;        /* hash set bins */
   unsigned long                              _size;       /* hash set size */

   /************************************************************************
    * Hash set helper functions.
    ************************************************************************/

   /*
    * Add an element to the set.
    */
   void add_item(T&);

   /*
    * Remove an element from the set.
    */
   void remove_item(const T&);
};

/*
 * Hash set iterator.
 */
template <typename T, typename H, typename Syn = unsynchronized>
class hash_set_iterator : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit hash_set_iterator(const hash_set<T,H,Syn>&);

   /*
    * Copy constructor.
    */
   hash_set_iterator(const hash_set_iterator<T,H,Syn>&);

   /*
    * Destructor.
    */
   virtual ~hash_set_iterator();

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
   const hash_set<T,H,Syn>&                                       _h;         /* hash set being iterated over */
   auto_read_lock<const Syn>                                      _rlock;     /* read lock on hash set */
   auto_ptr< set_iterator<typename hash_set<T,H,Syn>::hash_bin> > _bin_iter;  /* iterator over bins */
   auto_ptr< list_iterator<T> >                                   _item_iter; /* iterator over items in current bin */

   /*
    * Update bin and item iterators for next item.
    */
   void next_item();
};

/*
 * Hash set reverse iterator.
 */
template <typename T, typename H, typename Syn = unsynchronized>
class hash_set_iterator_reverse : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit hash_set_iterator_reverse(const hash_set<T,H,Syn>&);

   /*
    * Copy constructor.
    */
   hash_set_iterator_reverse(const hash_set_iterator_reverse<T,H,Syn>&);

   /*
    * Destructor.
    */
   virtual ~hash_set_iterator_reverse();

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
   const hash_set<T,H,Syn>&                                               _h;          /* hash set being iterated over */
   auto_read_lock<const Syn>                                              _rlock;      /* read lock on hash set */
   auto_ptr< set_iterator_reverse<typename hash_set<T,H,Syn>::hash_bin> > _bin_iter;   /* iterator over bins */
   auto_ptr< list_iterator_reverse<T> >                                   _item_iter;  /* iterator over items in current bin */

   /*
    * Update bin and item iterators for next item.
    */
   void next_item();
};
   
/***************************************************************************
 * Hash set iterator implementation.
 ***************************************************************************/

/*
 * Constructors.
 * Initialize position.
 */
template <typename T, typename H, typename Syn>
hash_set_iterator<T,H,Syn>::hash_set_iterator(
   const hash_set<T,H,Syn>& h)
 : _h(h),
   _rlock(_h),
   _bin_iter(
      new set_iterator<typename hash_set<T,H,Syn>::hash_bin>(*(_h._set))
   ),
   _item_iter(
      (_bin_iter->has_next()) ?
         (new list_iterator<T>(_bin_iter->next().items)) : NULL
   )
{
   this->next_item();
}

template <typename T, typename H, typename Syn>
hash_set_iterator_reverse<T,H,Syn>::hash_set_iterator_reverse(
   const hash_set<T,H,Syn>& h)
 : _h(h),
   _rlock(_h),
   _bin_iter(
      new set_iterator_reverse<typename hash_set<T,H,Syn>::hash_bin>(
         *(_h._set)
      )
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
template <typename T, typename H, typename Syn>
hash_set_iterator<T,H,Syn>::hash_set_iterator(
   const hash_set_iterator<T,H,Syn>& i)
 : _h(i._h),
   _rlock(_h),
   _bin_iter(
      new set_iterator<typename hash_set<T,H,Syn>::hash_bin>(*(i._bin_iter))
   ),
   _item_iter(
      (i._item_iter.get() == NULL) ?
         NULL : (new list_iterator<T>(*(i._item_iter)))
   )
{ }

template <typename T, typename H, typename Syn>
hash_set_iterator_reverse<T,H,Syn>::hash_set_iterator_reverse(
   const hash_set_iterator_reverse<T,H,Syn>& i)
 : _h(i._h),
   _rlock(_h),
   _bin_iter(
      new set_iterator_reverse<typename hash_set<T,H,Syn>::hash_bin>(
         *(i._bin_iter))
   ),
   _item_iter(
      (i._item_iter.get() == NULL) ?
         NULL : (new list_iterator_reverse<T>(*(i._item_iter)))
   )
{ }

/*
 * Destructors.
 * Do nothing as the hash set is automatically unlocked upon destruction
 * of the read lock.
 */
template <typename T, typename H, typename Syn>
hash_set_iterator<T,H,Syn>::~hash_set_iterator() {
   /* do nothing */
}

template <typename T, typename H, typename Syn>
hash_set_iterator_reverse<T,H,Syn>::~hash_set_iterator_reverse() {
   /* do nothing */
}

/*
 * Check if there is another item available.
 */
template <typename T, typename H, typename Syn>
bool hash_set_iterator<T,H,Syn>::has_next() const {
   return (_item_iter.get() != NULL);
}
   
template <typename T, typename H, typename Syn>
bool hash_set_iterator_reverse<T,H,Syn>::has_next() const {
   return (_item_iter.get() != NULL);
}

/*
 * Return the next item.
 * Throw an exception (ex_not_found) if there are no more items.
 */
template <typename T, typename H, typename Syn>
T& hash_set_iterator<T,H,Syn>::next() {
   if (_item_iter.get() != NULL) {
      T& t = _item_iter->next();
      this->next_item();
      return t;
   } else {
      throw ex_not_found(
         "iterator attempted to read past end of hash_set"
      );
   }
}

template <typename T, typename H, typename Syn>
T& hash_set_iterator_reverse<T,H,Syn>::next() {
   if (_item_iter.get() != NULL) {
      T& t = _item_iter->next();
      this->next_item();
      return t;
   } else {
      throw ex_not_found(
         "reverse iterator attempted to read past start of hash_set"
      );
   }
}

/*
 * Update bin and item iterators for next item.
 */
template <typename T, typename H, typename Syn>
void hash_set_iterator<T,H,Syn>::next_item() {
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

template <typename T, typename H, typename Syn>
void hash_set_iterator_reverse<T,H,Syn>::next_item() {
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
 * Hash set implementation.
 ***************************************************************************/
 
/*
 * Constructor.
 * Return an empty hash set.
 * Use the specified functors for hash set operations.
 */
template <typename T, typename H, typename Syn>
hash_set<T,H,Syn>::hash_set(
   const equalable_functor<T>&        f_eq, 
   const comparable_functor<H>&       f_c,
   const mappable_functor<const T,H>& f_h) 
 : Syn(),
   _f_equal(f_eq), 
   _f_compare(f_c), 
   _f_hash(f_h), 
   _set(new set<hash_bin>(_f_compare)),
   _size(0)
{ }

/*
 * Constructor.
 * Create a hash set from the given collection.
 * Use the specified functors for hash set operations.
 */
template <typename T, typename H, typename Syn>
hash_set<T,H,Syn>::hash_set(
   const collection<T>&               c, 
   const equalable_functor<T>&        f_eq, 
   const comparable_functor<H>&       f_c,
   const mappable_functor<const T,H>& f_h) 
 : Syn(),
   _f_equal(f_eq), 
   _f_compare(f_c), 
   _f_hash(f_h), 
   _set(new set<hash_bin>(_f_compare)),
   _size(0)
{
   auto_ptr< iterator<T> > i = c.iter_create();
   while (i->has_next())
      this->add_item(i->next());
}

/*
 * Copy constructor.
 */
template <typename T, typename H, typename Syn>
hash_set<T,H,Syn>::hash_set(const hash_set<T,H,Syn>& h)
 : Syn(), 
   _f_equal(h.f_equal),
   _f_compare(h.f_compare),
   _f_hash(h.f_hash),
   _set(new set<hash_bin>(_f_compare)),
   _size(0)
{
   typename hash_set<T,H,Syn>::iterator_t i(h);
   while (i.has_next())
      this->add_item(i.next());
}
     
/*
 * Destructor.
 */
template <typename T, typename H, typename Syn>
hash_set<T,H,Syn>::~hash_set() {
   /* do nothing */
}

/*
 * Add an element to the set.
 */
template <typename T, typename H, typename Syn>
void hash_set<T,H,Syn>::add_item(T& t) {
   auto_ptr<hash_bin> bin(new auto_hash_bin(_f_hash(t)));
   if (_set->contains(*bin)) {
      hash_bin& bin_match = _set->find(*bin);
      bool found = false;
      for (typename list<T>::iterator_t i(bin_match.items);
          ((i.has_next()) && (!found)); )
         found = _f_equal(t, i.next());
      if (!found) {
         bin_match.items.add(t);
         _size++;
      }
   } else {
      bin->items.add(t);
      _set->add(*bin);
      bin.release();
      _size++;
   }
}

/*
 * Add an element to the hash set.
 * Return a reference to the hash set.
 */
template <typename T, typename H, typename Syn>
hash_set<T,H,Syn>& hash_set<T,H,Syn>::add(T& t) {
   auto_write_lock<const Syn> wlock(*this);
   this->add_item(t);
   return *this;
}

/*
 * Add multiple elements to the hash set.
 * Return a reference to the hash set.
 */
template <typename T, typename H, typename Syn>
hash_set<T,H,Syn>& hash_set<T,H,Syn>::add(const collection<T>& c) {
   list<T> lst(c);
   auto_write_lock<const Syn> wlock(*this);
   for (typename list<T>::iterator_t i(lst); i.has_next(); )
      this->add_item(i.next());
   return *this;
}

/*
 * Remove an element from the set.
 */
template <typename T, typename H, typename Syn>
void hash_set<T,H,Syn>::remove_item(const T& t) {
   auto_hash_bin bin(_f_hash(t));
   if (_set->contains(bin)) {
      /* remove the element */
      hash_bin& bin_match = _set->find(bin);
      unsigned long n_items = bin_match.items.size();
      for (unsigned long n = 0; n < n_items; n++) {
         T& item = bin_match.items.remove_head();
         if (_f_equal(t, item)) {
            _size--;
            break;
         } else {
            bin_match.items.append(item);
         }
      }
      /* remove bin if empty after element removal */
      if (bin_match.items.is_empty()) {
         _set->remove(bin_match);
         delete &bin_match;
      }
   }
}

/*
 * Remove an element from the hash set.
 * Return a reference to the hash set.
 */
template <typename T, typename H, typename Syn>
hash_set<T,H,Syn>& hash_set<T,H,Syn>::remove(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   this->remove_item(t);
   return *this;
}

/*
 * Remove multiple elements from the hash set.
 * Return a reference to the hash set.
 */
template <typename T, typename H, typename Syn>
hash_set<T,H,Syn>& hash_set<T,H,Syn>::remove(const collection<T>& c) {
   list<T> lst(c);
   auto_write_lock<const Syn> wlock(*this);
   for (typename list<T>::iterator_t i(lst); i.has_next(); )
      this->remove_item(i.next());
   return *this;
}

/*
 * Remove all element(s) with the given hash from the set.
 * Return a reference to the set.
 */
template <typename T, typename H, typename Syn>
hash_set<T,H,Syn>& hash_set<T,H,Syn>::remove_hash(const H& h) {
   auto_write_lock<const Syn> wlock(*this);
   hash_bin bin(const_cast<H&>(h)); /* const_cast is safe here */
   if (_set->contains(bin)) {
      hash_bin& bin_match = _set->find(bin);
      _set->remove(bin_match);
      _size -= bin_match.items.size();
      delete &bin_match;
   }
   return *this;
}

/*
 * Remove all element(s) from the hash set.
 * Return a reference to the set.
 */
template <typename T, typename H, typename Syn>
hash_set<T,H,Syn>& hash_set<T,H,Syn>::clear() {
   auto_write_lock<const Syn> wlock(*this);
   _set.reset(new set<hash_bin>(_f_compare));
   _size = 0;
   return *this;
}

/*
 * Check if an element is in the hash set.
 */
template <typename T, typename H, typename Syn>
bool hash_set<T,H,Syn>::contains(const T& t) const {
   bool found = false;
   auto_read_lock<const Syn> rlock(*this);
   auto_hash_bin bin(_f_hash(t));
   if (_set->contains(bin)) {
      hash_bin& bin_match = _set->find(bin);
      for (typename list<T>::iterator_t i(bin_match.items);
           ((i.has_next()) && (!found)); )
      {
         found = _f_equal(t, i.next());
      }
   }
   return found;
}

/*
 * Find an element in the hash set that is equal to the given element.
 * Throw an exception (ex_not_found) if no such element is in the set.
 */
template <typename T, typename H, typename Syn>
T& hash_set<T,H,Syn>::find(const T& t) const {
   auto_read_lock<const Syn> rlock(*this);
   auto_hash_bin bin(_f_hash(t));
   hash_bin& bin_match = _set->find(bin);
   for (typename list<T>::iterator_t i(bin_match.items); i.has_next(); ) {
      T& item = i.next();
      if (_f_equal(t, item))
         return item;
   }
   throw ex_not_found(
      "item not found in hash_set"
   );
}

/*
 * Check if the set contains an element with the given hash.
 */
template <typename T, typename H, typename Syn>
bool hash_set<T,H,Syn>::contains_hash(const H& h) const {
   auto_read_lock<const Syn> rlock(*this);
   hash_bin bin(const_cast<H&>(h)); /* const_cast is safe here */
   return _set->contains(bin);
}

/*
 * Find an element with the given hash.
 * Throw an exception (ex_not_found) if no such element is in the set.
 */
template <typename T, typename H, typename Syn>
T& hash_set<T,H,Syn>::find_hash(const H& h) const {
   auto_read_lock<const Syn> rlock(*this);
   hash_bin bin(const_cast<H&>(h)); /* const_cast is safe here */
   hash_bin& bin_match = _set->find(bin);
   if (bin_match.items.is_empty()) {
      throw ex_not_found(
         "item not found in hash_set"
      );
   } else {
      return bin_match.items.head();
   }
}

/*
 * Find all elements in the set with the same hash as the given element.
 */
template <typename T, typename H, typename Syn>
list<T> hash_set<T,H,Syn>::find_bin(const T& t) const {
   auto_read_lock<const Syn> rlock(*this);
   auto_hash_bin bin(_f_hash(t));
   if (_set->contains(bin)) {
      hash_bin& bin_match = _set->find(bin);
      return bin_match.items;
   } else {
      return list<T>();
   }
}

/*
 * Find all elements in the set with the given hash.
 */
template <typename T, typename H, typename Syn>
list<T> hash_set<T,H,Syn>::find_bin_hash(const H& h) const {
   auto_read_lock<const Syn> rlock(*this);
   hash_bin bin(const_cast<H&>(h)); /* const_cast is safe here */
   if (_set->contains(bin)) {
      hash_bin& bin_match = _set->find(bin);
      return bin_match.items;
   } else {
      return list<T>();
   }
}

/*
 * Get number of elements in hash set.
 */
template <typename T, typename H, typename Syn>
unsigned long hash_set<T,H,Syn>::size() const {
   auto_read_lock<const Syn> rlock(*this);
   return _size;
}

/*
 * Return iterator over elements.
 */
template <typename T, typename H, typename Syn>
auto_ptr< iterator<T> > hash_set<T,H,Syn>::iter_create() const {
   return auto_ptr< iterator<T> >(
      new hash_set_iterator<T,H,Syn>(*this)
   );
}

/*
 * Return reverse iterator over elements.
 */
template <typename T, typename H, typename Syn>
auto_ptr< iterator<T> > hash_set<T,H,Syn>::iter_reverse_create() const {
   return auto_ptr< iterator<T> >(
      new hash_set_iterator_reverse<T,H,Syn>(*this)
   );
}

} /* namespace collections */

#endif
