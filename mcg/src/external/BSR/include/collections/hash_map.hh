/*
 * Hash maps (thread-safe).
 *
 * Hash maps do not allow duplicate items in the domain.
 *
 * Hash maps can be automatically created for domain objects which implement
 * the equalable and hashable interfaces and hash to objects which implement
 * the comparable interface.  Hash maps for other types may be constructed by 
 * specifying appropriate functors at creation time.
 *
 * NOTE: Non-const objects returned by the hash functor are owned by the 
 * hash map and will be deleted upon destruction of the hash map.
 */
#ifndef COLLECTIONS__HASH_MAP_HH
#define COLLECTIONS__HASH_MAP_HH

#include "collections/abstract/collection.hh"
#include "collections/abstract/hash_map.hh"
#include "collections/abstract/map.hh"
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
#include "lang/exceptions/ex_invalid_argument.hh"
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
using lang::exceptions::ex_invalid_argument;
using lang::exceptions::ex_not_found;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

/*
 * Declare classes for iterators over hash maps.
 */
template <typename T, typename U, typename H, typename Syn>
class hash_map_iterator;

template <typename T, typename U, typename H, typename Syn>
class hash_map_iterator_reverse;

/*
 * Hash maps.
 */
template <typename T, typename U, typename H, typename Syn = unsynchronized>
class hash_map : public abstract::hash_map<T,U,H>,
                 protected Syn {
public:
   /*
    * Friend classes.
    */
   friend class hash_map_iterator<T,U,H,Syn>;
   friend class hash_map_iterator_reverse<T,U,H,Syn>;

   /*
    * Define the iterator types.
    */
   typedef hash_map_iterator<T,U,H,Syn> iterator_t;
   typedef hash_map_iterator_reverse<T,U,H,Syn> iterator_reverse_t;

   /*
    * Constructor.
    * Return an empty hash map.
    * Optionally specify custom functors for hash map operations.
    */
   explicit hash_map(
      const equalable_functor<T>&        = equal_functors<T>::f_equal(), 
      const comparable_functor<H>&       = compare_functors<H>::f_compare(),
      const mappable_functor<const T,H>& = (hash_functors<T,H>::f_hash())
   );
   
   /*
    * Constructor.
    * Create a hash map containing the mappings t -> u for corresponding
    * elements of the given collections.  Optionally specify custom functors
    * for hash map operations.
    */
   explicit hash_map(
      const collection<T>&,
      const collection<U>&,
      const equalable_functor<T>&        = equal_functors<T>::f_equal(), 
      const comparable_functor<H>&       = compare_functors<H>::f_compare(),
      const mappable_functor<const T,H>& = (hash_functors<T,H>::f_hash())
   );

   /*
    * Constructor.
    * Create a hash map from the given map.
    * Optionally specify custom functors for hash map operations.
    */
   explicit hash_map(
      const abstract::map<T,U>&, 
      const equalable_functor<T>&        = equal_functors<T>::f_equal(), 
      const comparable_functor<H>&       = compare_functors<H>::f_compare(),
      const mappable_functor<const T,H>& = (hash_functors<T,H>::f_hash())
   );
      
   /*
    * Copy constructor.
    */
   hash_map(const hash_map<T,U,H,Syn>&);

   /*
    * Destructor.
    */
   virtual ~hash_map();
   
   /*
    * Add the mapping t -> u, overwriting any existing mapping for t.
    * Return a reference to the map.
    */
   hash_map<T,U,H>& add(T& /* t */, U& /* u */);

   /*
    * Add the mapping t -> u for corresponding elements of the collections.
    * Overwriting previous mappings if they exist.
    * The collections must have the same size.
    * Return a reference to the map.
    */
   hash_map<T,U,H>& add(const collection<T>&, const collection<U>&);

   /*
    * Add all mappings in the given map to the current map.
    * Overwriting previous mappings if they exist.
    * Return a reference to the map.
    */
   hash_map<T,U,H>& add(const abstract::map<T,U>&);
   
   /*
    * Remove element(s) and their corresponding images from the map.
    * Return a reference to the map.
    */
   hash_map<T,U,H>& remove(const T&);
   hash_map<T,U,H>& remove(const collection<T>&);
   
   /*
    * Remove all element(s) with the given hash and their images from the map.
    * Return a reference to the map.
    */
   hash_map<T,U,H>& remove_hash(const H&);

   /*
    * Clear the map, resetting it to the empty map.
    * Return a reference to the map.
    */
   hash_map<T,U,H>& clear();

   /*
    * Search.
    * Return the element in the domain matching the given element.
    * Throw an exception (ex_not_found) when attempting to find an element not
    * contained in the map.
    */
   bool contains(const T&) const;
   T& find(const T&) const;
   
   /*
    * Search.
    * Find an element in the domain with the given hash.
    * Throw an exception (ex_not_found) when attempting to find an element not
    * contained in the set.
    */
   bool contains_hash(const H&) const;
   T& find_hash(const H&) const;

   /*
    * Search.
    * Find all elements in the domain of the map with the same hash as the 
    * given element.
    */
   list<T> find_bin(const T&) const;

   /*
    * Search.
    * Find all elements in the domain of the map with the given hash.
    */
   list<T> find_bin_hash(const H&) const;

   /*
    * Search.
    * Return the image of the given domain element under the map.
    * Throw an exception (ex_not_found) when attempting to find the image of
    * an element not contained in the map.
    */
   U& find_image(const T&) const;

   /*
    * Size.
    * Return the number of elements in the domain.
    */
   unsigned long size() const;
   
   /*
    * Return a list of elements in the domain of the map and a list of their 
    * corresponding images in the range.  Hence, range() does not necessarily 
    * return a unique list of elements as two or more elements of the domain 
    * might have the same image.
    */
   list<T> domain() const;
   list<U> range() const;

   /*
    * Return lists of elements in the domain and range.
    * This method allows one to atomically capture the contents of a hash map
    * if there is a potential for the hash map to change between calls of the
    * above domain() and range() methods.
    */
   void contents(
      auto_ptr< collections::list<T> >&,  /* domain */
      auto_ptr< collections::list<U> >&   /* range */
   ) const;

   /*
    * Return iterators over elements in the domain.
    */
   auto_ptr< iterator<T> > iter_create() const;
   auto_ptr< iterator<T> > iter_reverse_create() const;

protected:
   /************************************************************************
    * Hash map data structures.
    ************************************************************************/

   /*
    * Bin in the hash map.
    */
   class hash_bin {
   public:
      /*
       * Constructor.
       */
      explicit hash_bin(H& h)
       : hash_value(h),
         items_domain(),
         items_range()
      { }

      /*
       * Destructor.
       */
      virtual ~hash_bin() { /* do nothing */ }

      /*
       * Bin data.
       */
      H& hash_value;          /* hash value of domain items in the bin */
      list<T> items_domain;   /* domain items sharing the hash value */
      list<U> items_range;    /* corresponding items in the range */

   protected:
      /*
       * Copy constructor.
       * WARNING: Does not create a deep copy of the hash value.
       */
      explicit hash_bin(const hash_bin& bin)
       : hash_value(bin.hash_value), 
         items_domain(bin.items_domain),
         items_range(bin.items_range)
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
    * Hash map data.
    */
   const equalable_functor<T>&                _f_equal;    /* equality function */
   const hash_bin_compare_functor             _f_compare;  /* comparison function */
   const mappable_functor<const T,H>&         _f_hash;     /* hash function */
   auto_collection< hash_bin, set<hash_bin> > _set;        /* hash map bins */
   unsigned long                              _size;       /* hash map size */

   /************************************************************************
    * Hash map helper functions.
    ************************************************************************/

   /*
    * Add a mapping to the hash map.
    */
   void add_mapping(T&, U&);

   /*
    * Remove an element and its image from the hash map.
    */
   void remove_item(const T&);
};

/*
 * Hash map iterator.
 */
template <typename T, typename U, typename H, typename Syn = unsynchronized>
class hash_map_iterator : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit hash_map_iterator(const hash_map<T,U,H,Syn>&);

   /*
    * Copy constructor.
    */
   hash_map_iterator(const hash_map_iterator<T,U,H,Syn>&);

   /*
    * Destructor.
    */
   virtual ~hash_map_iterator();

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
   const hash_map<T,U,H,Syn>&                                       _h;         /* hash map being iterated over */
   auto_read_lock<const Syn>                                        _rlock;     /* read lock on hash map */
   auto_ptr< set_iterator<typename hash_map<T,U,H,Syn>::hash_bin> > _bin_iter;  /* iterator over bins */
   auto_ptr< list_iterator<T> >                                     _item_iter; /* iterator over items in current bin */

   /*
    * Update bin and item iterators for next item.
    */
   void next_item();
};

/*
 * Hash map reverse iterator.
 */
template <typename T, typename U, typename H, typename Syn = unsynchronized>
class hash_map_iterator_reverse : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit hash_map_iterator_reverse(const hash_map<T,U,H,Syn>&);

   /*
    * Copy constructor.
    */
   hash_map_iterator_reverse(const hash_map_iterator_reverse<T,U,H,Syn>&);

   /*
    * Destructor.
    */
   virtual ~hash_map_iterator_reverse();

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
   const hash_map<T,U,H,Syn>&                                               _h;          /* hash map being iterated over */
   auto_read_lock<const Syn>                                                _rlock;      /* read lock on hash map */
   auto_ptr< set_iterator_reverse<typename hash_map<T,U,H,Syn>::hash_bin> > _bin_iter;   /* iterator over bins */
   auto_ptr< list_iterator_reverse<T> >                                     _item_iter;  /* iterator over items in current bin */

   /*
    * Update bin and item iterators for next item.
    */
   void next_item();
};
   
/***************************************************************************
 * Hash map iterator implementation.
 ***************************************************************************/

/*
 * Constructors.
 * Initialize position.
 */
template <typename T, typename U, typename H, typename Syn>
hash_map_iterator<T,U,H,Syn>::hash_map_iterator(
   const hash_map<T,U,H,Syn>& h)
 : _h(h),
   _rlock(_h),
   _bin_iter(
      new set_iterator<typename hash_map<T,U,H,Syn>::hash_bin>(*(_h._set))
   ),
   _item_iter(
      (_bin_iter->has_next()) ?
         (new list_iterator<T>(_bin_iter->next().items_domain)) : NULL
   )
{
   this->next_item();
}

template <typename T, typename U, typename H, typename Syn>
hash_map_iterator_reverse<T,U,H,Syn>::hash_map_iterator_reverse(
   const hash_map<T,U,H,Syn>& h)
 : _h(h),
   _rlock(_h),
   _bin_iter(
      new set_iterator_reverse<typename hash_map<T,U,H,Syn>::hash_bin>(
         *(_h._set)
      )
   ),
   _item_iter(
      (_bin_iter->has_next()) ?
         (new list_iterator_reverse<T>(_bin_iter->next().items_domain)) : NULL
   )
{
   this->next_item();
}

/*
 * Copy constructors.
 */
template <typename T, typename U, typename H, typename Syn>
hash_map_iterator<T,U,H,Syn>::hash_map_iterator(
   const hash_map_iterator<T,U,H,Syn>& i)
 : _h(i._h),
   _rlock(_h),
   _bin_iter(
      new set_iterator<typename hash_map<T,U,H,Syn>::hash_bin>(*(i._bin_iter))
   ),
   _item_iter(
      (i._item_iter.get() == NULL) ? 
         NULL : (new list_iterator<T>(*(i._item_iter)))
   )
{ }

template <typename T, typename U, typename H, typename Syn>
hash_map_iterator_reverse<T,U,H,Syn>::hash_map_iterator_reverse(
   const hash_map_iterator_reverse<T,U,H,Syn>& i)
 : _h(i._h),
   _rlock(_h),
   _bin_iter(
      new set_iterator_reverse<typename hash_map<T,U,H,Syn>::hash_bin>(
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
 * Do nothing as the hash map is automatically unlocked upon destruction
 * of the read lock.
 */
template <typename T, typename U, typename H, typename Syn>
hash_map_iterator<T,U,H,Syn>::~hash_map_iterator() {
   /* do nothing */
}

template <typename T, typename U, typename H, typename Syn>
hash_map_iterator_reverse<T,U,H,Syn>::~hash_map_iterator_reverse() {
   /* do nothing */
}

/*
 * Check if there is another item available.
 */
template <typename T, typename U, typename H, typename Syn>
bool hash_map_iterator<T,U,H,Syn>::has_next() const {
   return (_item_iter.get() != NULL);
}
   
template <typename T, typename U, typename H, typename Syn>
bool hash_map_iterator_reverse<T,U,H,Syn>::has_next() const {
   return (_item_iter.get() != NULL);
}

/*
 * Return the next item.
 * Throw an exception (ex_not_found) if there are no more items.
 */
template <typename T, typename U, typename H, typename Syn>
T& hash_map_iterator<T,U,H,Syn>::next() {
   if (_item_iter.get() != NULL) {
      T& t = _item_iter->next();
      this->next_item();
      return t;
   } else {
      throw ex_not_found(
         "iterator attempted to read past end of hash_map"
      );
   }
}

template <typename T, typename U, typename H, typename Syn>
T& hash_map_iterator_reverse<T,U,H,Syn>::next() {
   if (_item_iter.get() != NULL) {
      T& t = _item_iter->next();
      this->next_item();
      return t;
   } else {
      throw ex_not_found(
         "reverse iterator attempted to read past start of hash_map"
      );
   }
}

/*
 * Update bin and item iterators for next item.
 */
template <typename T, typename U, typename H, typename Syn>
void hash_map_iterator<T,U,H,Syn>::next_item() {
   if (_item_iter.get() != NULL) {
      /* check if all items in current bin visited */
      while (!(_item_iter->has_next())) {
         /* move to next bin */
         if (_bin_iter->has_next()) {
            _item_iter.reset(
               new list_iterator<T>(_bin_iter->next().items_domain)
            );
         } else {
            _item_iter.reset(NULL);
            break;
         }
      }
   }
}

template <typename T, typename U, typename H, typename Syn>
void hash_map_iterator_reverse<T,U,H,Syn>::next_item() {
   if (_item_iter.get() != NULL) {
      /* check if all items in current bin visited */
      while (!(_item_iter->has_next())) {
         /* move to next bin */
         if (_bin_iter->has_next()) {
            _item_iter.reset(
               new list_iterator_reverse<T>(_bin_iter->next().items_domain)
            );
         } else {
            _item_iter.reset(NULL);
            break;
         }
      }
   }
}

/***************************************************************************
 * Hash map implementation.
 ***************************************************************************/

/*
 * Constructor.
 * Return an empty hash map.
 * Optionally specify custom functors for hash map operations.
 */
template <typename T, typename U, typename H, typename Syn>
hash_map<T,U,H,Syn>::hash_map(
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
 * Create a hash map containing the mappings t -> u for corresponding
 * elements of the given collections.  Optionally specify custom functors
 * for hash map operations.
 */
template <typename T, typename U, typename H, typename Syn>
hash_map<T,U,H,Syn>::hash_map(
   const collection<T>&               c_t,
   const collection<U>&               c_u,
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
   this->add(c_t, c_u);
}

/*
 * Constructor.
 * Create a hash map from the given collection.
 * Optionally specify custom functors for hash map operations.
 */
template <typename T, typename U, typename H, typename Syn>
hash_map<T,U,H,Syn>::hash_map(
   const abstract::map<T,U>&          m, 
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
   this->add(m);
}
    
/*
 * Copy constructor.
 */
template <typename T, typename U, typename H, typename Syn>
hash_map<T,U,H,Syn>::hash_map(const hash_map<T,U,H,Syn>& m)
 : Syn(), 
   _f_equal(m.f_equal),
   _f_compare(m.f_compare),
   _f_hash(m.f_hash),
   _set(new set<hash_bin>(_f_compare)),
   _size(0)
{
   this->add(m);
}

/*
 * Destructor.
 */
template <typename T, typename U, typename H, typename Syn>
hash_map<T,U,H,Syn>::~hash_map() {
   /* do nothing */
}

/*
 * Add a mapping to the hash map.
 */
template <typename T, typename U, typename H, typename Syn>
void hash_map<T,U,H,Syn>::add_mapping(T& t, U& u) {
   auto_ptr<hash_bin> bin(new auto_hash_bin(_f_hash(t)));
   if (_set->contains(*bin)) {
      hash_bin& bin_match = _set->find(*bin);
      bool found = false;
      unsigned long n_items = bin_match.items_domain.size();
      for (unsigned long n = 0; (n < n_items) && (!found); n++) {
         T& item_domain = bin_match.items_domain.remove_head();
         U& item_range  = bin_match.items_range.remove_head();
         if (_f_equal(t, item_domain)) {
            bin_match.items_domain.append(t);
            bin_match.items_range.append(u);
            found = true;
         } else {
            bin_match.items_domain.append(item_domain);
            bin_match.items_range.append(item_range);
         }
      }
      if (!found) {
         bin_match.items_domain.append(t);
         bin_match.items_range.append(u);
         _size++;
      }
   } else {
      bin->items_domain.append(t);
      bin->items_range.append(u);
      _set->add(*bin);
      bin.release();
      _size++;
   }
}

/*
 * Add the mapping t -> u, overwriting any existing mapping for t.
 * Return a reference to the map.
 */
template <typename T, typename U, typename H, typename Syn>
hash_map<T,U,H>& hash_map<T,U,H,Syn>::add(T& t, U& u) {
   auto_write_lock<const Syn> wlock(*this);
   this->add_mapping(t, u);
   return *this;
}

/*
 * Add the mapping t -> u for corresponding elements of the collections.
 * Overwriting previous mappings if they exist.
 * The collections must have the same size.
 * Return a reference to the map.
 */
template <typename T, typename U, typename H, typename Syn>
hash_map<T,U,H>& hash_map<T,U,H,Syn>::add(
   const collection<T>& c_t, 
   const collection<U>& c_u)
{
   list<T> lst_t(c_t);
   list<U> lst_u(c_u);
   if (c_t.size() != c_u.size())
      throw ex_invalid_argument(
         "collections being add to hash_map must have the same size"
      );
   typename list<T>::iterator_t i_t(lst_t);
   typename list<U>::iterator_t i_u(lst_u);
   auto_write_lock<const Syn> wlock(*this);
   while (i_t.has_next())
      this->add_mapping(i_t.next(), i_u.next());
   return *this;
}

/*
 * Add all mappings in the given map to the current map.
 * Overwriting previous mappings if they exist.
 * Return a reference to the map.
 */
template <typename T, typename U, typename H, typename Syn>
hash_map<T,U,H>& hash_map<T,U,H,Syn>::add(const abstract::map<T,U>& m) {
   auto_ptr< list<T> > domain;
   auto_ptr< list<U> > range;
   m.contents(domain, range);
   return this->add(*domain, *range);
}

/*
 * Remove an element and its image from the hash map.
 */
template <typename T, typename U, typename H, typename Syn>
void hash_map<T,U,H,Syn>::remove_item(const T& t) {
   auto_hash_bin bin(_f_hash(t));
   if (_set->contains(bin)) {
      /* remove the element */
      hash_bin& bin_match = _set->find(bin);
      unsigned long n_items = bin_match.items_domain.size();
      for (unsigned long n = 0; n < n_items; n++) {
         T& item_domain = bin_match.items_domain.remove_head();
         U& item_range  = bin_match.items_range.remove_head();
         if (_f_equal(t, item_domain)) {
            _size--;
            break;
         } else {
            bin_match.items_domain.append(item_domain);
            bin_match.items_range.append(item_range);
         }
      }
      /* remove bin if empty after element removal */
      if (bin_match.items_domain.is_empty()) {
         _set->remove(bin_match);
         delete &bin_match;
      }
   }
}

/*
 * Remove element(s) and their corresponding images from the map.
 * Return a reference to the map.
 */
template <typename T, typename U, typename H, typename Syn>
hash_map<T,U,H>& hash_map<T,U,H,Syn>::remove(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   this->remove_item(t);
   return *this;
}

template <typename T, typename U, typename H, typename Syn>
hash_map<T,U,H>& hash_map<T,U,H,Syn>::remove(const collection<T>& c) {
   list<T> lst(c);
   auto_write_lock<const Syn> wlock(*this);
   for (typename list<T>::iterator_t i(lst); i.has_next(); )
      this->remove_item(i.next());
   return *this;
}

/*
 * Remove all element(s) with the given hash and their images from the map.
 * Return a reference to the map.
 */
template <typename T, typename U, typename H, typename Syn>
hash_map<T,U,H>& hash_map<T,U,H,Syn>::remove_hash(const H& h) {
   auto_write_lock<const Syn> wlock(*this);
   hash_bin bin(const_cast<H&>(h)); /* const_cast is safe here */
   if (_set->contains(bin)) {
      hash_bin& bin_match = _set->find(bin);
      _set->remove(bin_match);
      _size -= bin_match.items_domain.size();
      delete &bin_match;
   }
   return *this;
}

/*
 * Clear the map, resetting it to the empty map.
 * Return a reference to the map.
 */
template <typename T, typename U, typename H, typename Syn>
hash_map<T,U,H>& hash_map<T,U,H,Syn>::clear() {
   auto_write_lock<const Syn> wlock(*this);
   _set.reset(new set<hash_bin>(_f_compare));
   _size = 0;
   return *this;
}

/*
 * Search.
 * Return the element in the domain matching the given element.
 * Throw an exception (ex_not_found) when attempting to find an element not
 * contained in the map.
 */
template <typename T, typename U, typename H, typename Syn>
bool hash_map<T,U,H,Syn>::contains(const T& t) const {
   bool found = false;
   auto_read_lock<const Syn> rlock(*this);
   auto_hash_bin bin(_f_hash(t));
   if (_set->contains(bin)) {
      hash_bin& bin_match = _set->find(bin);
      for (typename list<T>::iterator_t i(bin_match.items_domain);
           ((i.has_next()) && (!found)); )
      {
         found = _f_equal(t, i.next());
      }
   }
   return found;
}

template <typename T, typename U, typename H, typename Syn>
T& hash_map<T,U,H,Syn>::find(const T& t) const {
   auto_read_lock<const Syn> rlock(*this);
   auto_hash_bin bin(_f_hash(t));
   hash_bin& bin_match = _set->find(bin);
   for (typename list<T>::iterator_t i(bin_match.items_domain);
        i.has_next(); )
   {
      T& item = i.next();
      if (_f_equal(t, item))
         return item;
   }
   throw ex_not_found(
      "item not found in hash_map"
   );
}

/*
 * Search.
 * Find an element in the domain with the given hash.
 * Throw an exception (ex_not_found) when attempting to find an element not
 * contained in the set.
 */
template <typename T, typename U, typename H, typename Syn>
bool hash_map<T,U,H,Syn>::contains_hash(const H& h) const {
   auto_read_lock<const Syn> rlock(*this);
   hash_bin bin(const_cast<H&>(h)); /* const_cast is safe here */
   return _set->contains(bin);
}

template <typename T, typename U, typename H, typename Syn>
T& hash_map<T,U,H,Syn>::find_hash(const H& h) const {
   auto_read_lock<const Syn> rlock(*this);
   hash_bin bin(const_cast<H&>(h)); /* const_cast is safe here */
   hash_bin& bin_match = _set->find(bin);
   if (bin_match.items_domain.is_empty()) {
      throw ex_not_found(
         "item not found in hash_map"
      );
   } else {
      return bin_match.items_domain.head();
   }
}

/*
 * Search.
 * Find all elements in the domain of the map with the same hash as the 
 * given element.
 */
template <typename T, typename U, typename H, typename Syn>
list<T> hash_map<T,U,H,Syn>::find_bin(const T& t) const {
   auto_read_lock<const Syn> rlock(*this);
   auto_hash_bin bin(_f_hash(t));
   if (_set->contains(bin)) {
      hash_bin& bin_match = _set->find(bin);
      return bin_match.items_domain;
   } else {
      return list<T>();
   }
}

/*
 * Search.
 * Find all elements in the domain of the map with the given hash.
 */
template <typename T, typename U, typename H, typename Syn>
list<T> hash_map<T,U,H,Syn>::find_bin_hash(const H& h) const {
   auto_read_lock<const Syn> rlock(*this);
   hash_bin bin(const_cast<H&>(h)); /* const_cast is safe here */
   if (_set->contains(bin)) {
      hash_bin& bin_match = _set->find(bin);
      return bin_match.items_domain;
   } else {
      return list<T>();
   }
}

/*
 * Search.
 * Return the image of the given domain element under the map.
 * Throw an exception (ex_not_found) when attempting to find the image of an
 * element not contained in the map.
 */
template <typename T, typename U, typename H, typename Syn>
U& hash_map<T,U,H,Syn>::find_image(const T& t) const {
   auto_read_lock<const Syn> rlock(*this);
   auto_hash_bin bin(_f_hash(t));
   hash_bin& bin_match = _set->find(bin);
   typename list<T>::iterator_t i_domain(bin_match.items_domain);
   typename list<U>::iterator_t i_range(bin_match.items_range);
   while (i_domain.has_next()) {
      T& item_domain = i_domain.next();
      U& item_range  = i_range.next();
      if (_f_equal(t, item_domain))
         return item_range;
   }
   throw ex_not_found(
      "item not found in hash_map"
   );
}

/*
 * Size.
 * Return the number of elements in the domain.
 */
template <typename T, typename U, typename H, typename Syn>
unsigned long hash_map<T,U,H,Syn>::size() const {
   auto_read_lock<const Syn> rlock(*this);
   return _size;
}

/*
 * Return a list of elements in the domain of the map and a list of their 
 * corresponding images in the range.  Hence, range() does not necessarily 
 * return a unique list of elements as two or more elements of the domain 
 * might have the same image.
 */
template <typename T, typename U, typename H, typename Syn>
list<T> hash_map<T,U,H,Syn>::domain() const {
   list<T> lst;
   auto_read_lock<const Syn> rlock(*this);
   typename set<hash_bin>::iterator_t i(*_set);
   while (i.has_next())
      lst.add(i.next().items_domain);
   return lst;
}

template <typename T, typename U, typename H, typename Syn>
list<U> hash_map<T,U,H,Syn>::range() const {
   list<U> lst;
   auto_read_lock<const Syn> rlock(*this);
   typename set<hash_bin>::iterator_t i(*_set);
   while (i.has_next())
      lst.add(i.next().items_range);
   return lst;
}

/*
 * Return lists of elements in the domain and range.
 * This method allows one to atomically capture the contents of a hash map
 * if there is a potential for the hash map to change between calls of the
 * above domain() and range() methods.
 */
template <typename T, typename U, typename H, typename Syn>
void hash_map<T,U,H,Syn>::contents(
   auto_ptr< collections::list<T> >& domain,
   auto_ptr< collections::list<U> >& range) const
{
   domain.reset(new list<T>());
   range.reset(new list<U>());
   auto_read_lock<const Syn> rlock(*this);
   typename set<hash_bin>::iterator_t i(*_set);
   while (i.has_next()) {
      hash_bin& bin = i.next();
      domain->add(bin.items_domain);
      range->add(bin.items_range);
   }
}

/*
 * Return iterators over elements in the domain.
 */
template <typename T, typename U, typename H, typename Syn>
auto_ptr< iterator<T> > hash_map<T,U,H,Syn>::iter_create() const {
   return auto_ptr< iterator<T> >(
      new hash_map_iterator<T,U,H,Syn>(*this)
   );
}

template <typename T, typename U, typename H, typename Syn>
auto_ptr< iterator<T> > hash_map<T,U,H,Syn>::iter_reverse_create() const {
   return auto_ptr< iterator<T> >(
      new hash_map_iterator_reverse<T,U,H,Syn>(*this)
   );
}

} /* namespace collections */

#endif
