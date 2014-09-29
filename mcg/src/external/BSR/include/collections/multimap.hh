/*
 * Multimaps (thread-safe).
 *
 * Multimaps allow duplicate items in the domain.
 */
#ifndef COLLECTIONS__MULTIMAP_HH
#define COLLECTIONS__MULTIMAP_HH

#include "collections/abstract/collection.hh"
#include "collections/abstract/map.hh"
#include "collections/abstract/multimap.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "collections/set.hh"
#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "functors/comparable_functors.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
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
using lang::exceptions::ex_invalid_argument;
using lang::exceptions::ex_not_found;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

/*
 * Declare classes for iterators over multimaps.
 */
template <typename T, typename U, typename Syn>
class multimap_iterator;

template <typename T, typename U, typename Syn>
class multimap_iterator_reverse;

/*
 * Multimaps.
 */
template <typename T, typename U, typename Syn = unsynchronized>
class multimap : public abstract::multimap<T,U>,
                 protected Syn {
public:
   /*
    * Friend classes.
    */
   friend class multimap_iterator<T,U,Syn>;
   friend class multimap_iterator_reverse<T,U,Syn>;

   /*
    * Define the iterator types.
    */
   typedef multimap_iterator<T,U,Syn> iterator_t;
   typedef multimap_iterator_reverse<T,U,Syn> iterator_reverse_t;

   /*
    * Constructor.
    * Return an empty multimap.
    * Optionally specify the comparison function to use for domain elements.
    */
   multimap(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Constructor.
    * Create a multimap containing the mappings t -> u for corresponding 
    * elements of the given collections.  Optionally specify the comparison
    * function to use for domain elements.
    */
   explicit multimap(
      const collection<T>&, 
      const collection<U>&,
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Constructor.
    * Create a multimap with the same contents as the given map.
    * Optionally specify the comparison function to use for domain elements.
    */
   explicit multimap(
      const abstract::map<T,U>&,
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );
   
   /*
    * Constructor.
    * Create a multimap with the same contents as the given multimap.
    * Optionally specify the comparison function to use for domain elements.
    */
   explicit multimap(
      const abstract::multimap<T,U>&,
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Copy constructor.
    */
   multimap(const multimap<T,U,Syn>&);

   /*
    * Destructor.
    */
   virtual ~multimap();

   /*
    * Add the mapping t -> u.
    * Return a reference to the multimap.
    */
   multimap<T,U,Syn>& add(T& /* t */, U& /* u */);

   /*
    * Add the mapping t -> u for corresponding elements of the collections.
    * The collections must have the same size.
    * Return a reference to the multimap.
    */
   multimap<T,U,Syn>& add(const collection<T>&, const collection<U>&);

   /*
    * Add all mappings in the given map/multimap to the current multimap.
    * Return a reference to the multimap.
    */
   multimap<T,U,Syn>& add(const abstract::map<T,U>&);
   multimap<T,U,Syn>& add(const abstract::multimap<T,U>&);

   /*
    * Remove element(s) and their corresponding images from the multimap.
    * Return a reference to the multimap.
    */
   multimap<T,U,Syn>& remove(const T&);
   multimap<T,U,Syn>& remove(const collection<T>&);

   /*
    * Remove all instances of the element(s) from the multimap.
    * Return a reference to the multimap.
    */
   multimap<T,U,Syn>& remove_all(const T&);
   multimap<T,U,Syn>& remove_all(const collection<T>&);

   /*
    * Clear the multimap, resetting it to the empty multimap.
    * Return a reference to the multimap.
    */
   multimap<T,U,Syn>& clear();

   /*
    * Search.
    * Return an element in the domain matching the given element.
    * Throw an exception (ex_not_found) when attempting to find an element not 
    * contained in the map.
    */
   bool contains(const T&) const;
   T& find(const T&) const;

   /*
    * Search.
    * Return a list of all element(s) in the domain matching the given element.
    */
   list<T> find_all(const T&) const;

   /*
    * Search.
    * Return an image of the given domain element under the multimap.
    * Throw an exception (ex_not_found) when attempting to find the image of
    * an element not contained in the multimap.
    */
   U& find_image(const T&) const;

   /*
    * Search.
    * Return a list of all element(s) in the range which are images of the 
    * given domain element.
    */
   list<U> find_images(const T&) const;

   /*
    * Size.
    * Return the number of mappings.
    */
   unsigned long size() const;
   
   /*
    * Return a list of elements in the domain of the multimap and a list of 
    * their corresponding images in the range.  All mappings t -> u in the 
    * multimap appear as a pair of elements (t, u) at corresponding positions
    * in the lists.  Hence, both the domain and range list may have repeated 
    * elements.
    */
   list<T> domain() const;
   list<U> range() const;
   
   /*
    * Return lists of elements in the domain and range.
    * This method allows one to atomically capture the contents of a multimap
    * if there is a potential for the multimap to change between calls of the
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
    * Multimap data structures.
    ************************************************************************/

   /*
    * Bin containing mappings with equivalent (as defined by the comparison
    * functor) domain items.
    */
   class equiv_bin {
   public:
      /*
       * Constructor.
       * Initialize the bin to contain the given domain element.
       * This should only be used for search/comparison to other bins.
       */
      explicit equiv_bin(T& t)
       : items_domain(),
         items_range()
      {
         items_domain.add(t);
      }

      /*
       * Constructor.
       * Initialize the bin to contain the given mapping.
       */
      explicit equiv_bin(T& t, U& u)
       : items_domain(),
         items_range()
      {
         items_domain.add(t);
         items_range.add(u);
      }

      /*
       * Copy constructor.
       */
      equiv_bin(const equiv_bin& bin)
       : items_domain(bin.items_domain),
         items_range(bin.items_range)
      { }

      /*
       * Destructor.
       */
      ~equiv_bin() { /* do nothing */ }

      /*
       * Equivalent items.
       */
      list<T> items_domain;
      list<U> items_range;
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
         return _f(bin0.items_domain.head(), bin1.items_domain.head());
      }

   protected:
      const comparable_functor<T>& _f;
   };

   /*
    * Multimap data.
    */
   const equiv_bin_compare_functor              _f_compare; /* comparison function */
   auto_collection< equiv_bin, set<equiv_bin> > _set;       /* set of equivalence classes */
   unsigned long                                _size;      /* multimap size */

   /************************************************************************
    * Multiset helper functions.
    ************************************************************************/
 
   /*
    * Add a mapping to the multimap.
    */
   void add_mapping(T&, U&);

   /*
    * Remove an element and its image from the multimap.
    */
   void remove_item(const T&);

   /*
    * Remove all instances of an element (and its images) from the multimap.
    */
   void remove_item_all(const T&);
};

/*
 * Multimap iterator.
 */
template <typename T, typename U, typename Syn = unsynchronized>
class multimap_iterator : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit multimap_iterator(const multimap<T,U,Syn>&);

   /*
    * Copy constructor.
    */
   multimap_iterator(const multimap_iterator<T,U,Syn>&);

   /*
    * Destructor.
    */
   virtual ~multimap_iterator();

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
   const multimap<T,U,Syn>&                                        _m;         /* multimap being iterated over */
   auto_read_lock<const Syn>                                       _rlock;     /* read lock on set */
   auto_ptr< set_iterator<typename multimap<T,U,Syn>::equiv_bin> > _bin_iter;  /* iterator over bins */
   auto_ptr< list_iterator<T> >                                    _item_iter; /* iterator over items in current bin */

   /*
    * Update bin and item iterators for next item.
    */
   void next_item();
};

/*
 * Multimap reverse iterator.
 */
template <typename T, typename U, typename Syn = unsynchronized>
class multimap_iterator_reverse : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit multimap_iterator_reverse(const multimap<T,U,Syn>&);

   /*
    * Copy constructor.
    */
   multimap_iterator_reverse(const multimap_iterator_reverse<T,U,Syn>&);

   /*
    * Destructor.
    */
   virtual ~multimap_iterator_reverse();

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
   const multimap<T,U,Syn>&                                                _m;         /* multimap being iterated over */
   auto_read_lock<const Syn>                                               _rlock;     /* read lock on set */
   auto_ptr< set_iterator_reverse<typename multimap<T,U,Syn>::equiv_bin> > _bin_iter;  /* iterator over bins */
   auto_ptr< list_iterator_reverse<T> >                                    _item_iter; /* iterator over items in current bin */

   /*
    * Update bin and item iterators for next item.
    */
   void next_item();
};

/***************************************************************************
 * Multimap iterator implementation.
 ***************************************************************************/

/*
 * Constructors.
 * Initialize position.
 */
template <typename T, typename U, typename Syn>
multimap_iterator<T,U,Syn>::multimap_iterator(
   const multimap<T,U,Syn>& m)
 : _m(m),
   _rlock(_m),
   _bin_iter(
      new set_iterator<typename multimap<T,U,Syn>::equiv_bin>(*(_m._set))
   ),
   _item_iter(
      (_bin_iter->has_next()) ?
         (new list_iterator<T>(_bin_iter->next().items_domain)) : NULL
   )
{
   this->next_item();
}

template <typename T, typename U, typename Syn>
multimap_iterator_reverse<T,U,Syn>::multimap_iterator_reverse(
   const multimap<T,U,Syn>& m)
 : _m(m),
   _rlock(_m),
   _bin_iter(
      new set_iterator_reverse<typename multimap<T,U,Syn>::equiv_bin>(
         *(_m._set)
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
template <typename T, typename U, typename Syn>
multimap_iterator<T,U,Syn>::multimap_iterator(
   const multimap_iterator<T,U,Syn>& i)
 : _m(i._m),
   _rlock(_m),
   _bin_iter(
      new set_iterator<typename multimap<T,U,Syn>::equiv_bin>(*(i._bin_iter))
   ),
   _item_iter(
      (i._item_iter.get() == NULL) ?
         NULL : (new list_iterator<T>(*(i._item_iter)))
   )
{ }

template <typename T, typename U, typename Syn>
multimap_iterator_reverse<T,U,Syn>::multimap_iterator_reverse(
   const multimap_iterator_reverse<T,U,Syn>& i)
 : _m(i._m),
   _rlock(_m),
   _bin_iter(
      new set_iterator_reverse<typename multimap<T,U,Syn>::equiv_bin>(
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
 * Do nothing as the multimap is automatically unlocked upon destruction
 * of the read lock.
 */
template <typename T, typename U, typename Syn>
multimap_iterator<T,U,Syn>::~multimap_iterator() {
   /* do nothing */
}

template <typename T, typename U, typename Syn>
multimap_iterator_reverse<T,U,Syn>::~multimap_iterator_reverse() {
   /* do nothing */
}

/*
 * Check if there is another item available.
 */
template <typename T, typename U, typename Syn>
bool multimap_iterator<T,U,Syn>::has_next() const {
   return (_item_iter.get() != NULL);
}
   
template <typename T, typename U, typename Syn>
bool multimap_iterator_reverse<T,U,Syn>::has_next() const {
   return (_item_iter.get() != NULL);
}

/*
 * Return the next item.
 * Throw an exception (ex_not_found) if there are no more items.
 */
template <typename T, typename U, typename Syn>
T& multimap_iterator<T,U,Syn>::next() {
   if (_item_iter.get() != NULL) {
      T& t = _item_iter->next();
      this->next_item();
      return t;
   } else {
      throw ex_not_found(
         "iterator attempted to read past end of multimap"
      );
   }
}

template <typename T, typename U, typename Syn>
T& multimap_iterator_reverse<T,U,Syn>::next() {
   if (_item_iter.get() != NULL) {
      T& t = _item_iter->next();
      this->next_item();
      return t;
   } else {
      throw ex_not_found(
         "reverse iterator attempted to read past start of multimap"
      );
   }
}

/*
 * Update bin and item iterators for next item.
 */
template <typename T, typename U, typename Syn>
void multimap_iterator<T,U,Syn>::next_item() {
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

template <typename T, typename U, typename Syn>
void multimap_iterator_reverse<T,U,Syn>::next_item() {
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
 * Multimap implementation.
 ***************************************************************************/

/*
 * Constructor.
 * Return an empty multimap.
 * Specify the comparison function to use for domain elements.
 */
template <typename T, typename U, typename Syn>
multimap<T,U,Syn>::multimap(const comparable_functor<T>& f)
 : Syn(),
   _f_compare(f),
   _set(new set<equiv_bin>(_f_compare)),
   _size(0)
{ }

/*
 * Constructor.
 * Create a multimap containing the mappings t -> u for corresponding 
 * elements of the given collections.  Specify the comparison function
 * to use for domain elements.
 */
template <typename T, typename U, typename Syn>
multimap<T,U,Syn>::multimap(
   const collection<T>&         c_t, 
   const collection<U>&         c_u,
   const comparable_functor<T>& f)
 : Syn(),
   _f_compare(f),
   _set(new set<equiv_bin>(_f_compare)),
   _size(0)
{
   this->add(c_t, c_u);
}

/*
 * Constructor.
 * Create a multimap with the same contents as the given map.
 * Specify the comparison function to use for domain elements.
 */
template <typename T, typename U, typename Syn>
multimap<T,U,Syn>::multimap(
   const abstract::map<T,U>&    m,
   const comparable_functor<T>& f)
 : Syn(),
   _f_compare(f),
   _set(new set<equiv_bin>(_f_compare)),
   _size(0)
{
   this->add(m);
}

/*
 * Constructor.
 * Create a multimap with the same contents as the given multimap.
 * Specify the comparison function to use for domain elements.
 */
template <typename T, typename U, typename Syn>
multimap<T,U,Syn>::multimap(
   const abstract::multimap<T,U>& m,
   const comparable_functor<T>&   f)
 : Syn(),
   _f_compare(f),
   _set(new set<equiv_bin>(_f_compare)),
   _size(0)
{
   this->add(m);
}

/*
 * Copy constructor.
 */
template <typename T, typename U, typename Syn>
multimap<T,U,Syn>::multimap(const multimap<T,U,Syn>& m)
 : Syn(),
   _f_compare(m._f_compare),
   _set(new set<equiv_bin>(_f_compare)),
   _size(0)
{
   this->add(m);
}

/*
 * Destructor.
 */
template <typename T, typename U, typename Syn>
multimap<T,U,Syn>::~multimap() {
   /* do nothing */
}

/*
 * Add a mapping to the multimap.
 */
template <typename T, typename U, typename Syn>
void multimap<T,U,Syn>::add_mapping(T& t, U& u) {
   auto_ptr<equiv_bin> bin(new equiv_bin(t, u));
   if (_set->contains(*bin)) {
      equiv_bin& bin_match = _set->find(*bin);
      bin_match.items_domain.add(t);
      bin_match.items_range.add(u);
   } else {
      _set->add(*bin);
      bin.release();
   }
   _size++;
}

/*
 * Add the mapping t -> u.
 * Return a reference to the multimap.
 */
template <typename T, typename U, typename Syn>
multimap<T,U,Syn>& multimap<T,U,Syn>::add(T& t, U& u) {
   auto_write_lock<const Syn> wlock(*this);
   this->add_mapping(t, u);
   return *this;
}

/*
 * Add the mapping t -> u for corresponding elements of the collections.
 * The collections must have the same size.
 * Return a reference to the multimap.
 */
template <typename T, typename U, typename Syn>
multimap<T,U,Syn>& multimap<T,U,Syn>::add(
   const collection<T>& c_t,
   const collection<U>& c_u) 
{
   list<T> lst_t(c_t);
   list<U> lst_u(c_u);
   if (c_t.size() != c_u.size())
      throw ex_invalid_argument(
         "collections being add to multimap must have the same size"
      );
   typename list<T>::iterator_t i_t(lst_t);
   typename list<U>::iterator_t i_u(lst_u);
   auto_write_lock<const Syn> wlock(*this);
   while (i_t.has_next())
      this->add_mapping(i_t.next(), i_u.next());
   return *this;
}

/*
 * Add all mappings in the given map/multimap to the current multimap.
 * Return a reference to the multimap.
 */
template <typename T, typename U, typename Syn>
multimap<T,U,Syn>& multimap<T,U,Syn>::add(const abstract::map<T,U>& m) {
   auto_ptr< list<T> > domain;
   auto_ptr< list<U> > range;
   m.contents(domain, range);
   return this->add(*domain, *range);
}

template <typename T, typename U, typename Syn>
multimap<T,U,Syn>& multimap<T,U,Syn>::add(const abstract::multimap<T,U>& m) {
   auto_ptr< list<T> > domain;
   auto_ptr< list<U> > range;
   m.contents(domain, range);
   return this->add(*domain, *range);
}

/*
 * Remove an element and its image from the multimap.
 */
template <typename T, typename U, typename Syn>
void multimap<T,U,Syn>::remove_item(const T& t) {
   equiv_bin bin(const_cast<T&>(t));   /* const_cast is safe here */
   if (_set->contains(bin)) {
      /* remove one instance of the element */
      equiv_bin& bin_match = _set->find(bin);
      bin_match.items_domain.remove_head();
      bin_match.items_range.remove_head();
      _size--;
      /* remove bin if empty after element removal */
      if (bin_match.items_domain.is_empty()) {
         _set->remove(bin_match);
         delete &bin_match;
      }
   }
}

/*
 * Remove element(s) and their corresponding images from the multimap.
 * Return a reference to the multimap.
 */
template <typename T, typename U, typename Syn>
multimap<T,U,Syn>& multimap<T,U,Syn>::remove(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   this->remove_item(t);
   return *this;
}

template <typename T, typename U, typename Syn>
multimap<T,U,Syn>& multimap<T,U,Syn>::remove(const collection<T>& c) {
   list<T> lst(c);
   auto_write_lock<const Syn> wlock(*this);
   for (typename list<T>::iterator_t i(lst); i.has_next(); )
      this->remove_item(i.next());
   return *this;
}

/*
 * Remove all instances of an element (and its images) from the multimap.
 */
template <typename T, typename U, typename Syn>
void multimap<T,U,Syn>::remove_item_all(const T& t) {
   equiv_bin bin(const_cast<T&>(t));   /* const_cast is safe here */
   if (_set->contains(bin)) {
      equiv_bin& bin_match = _set->find(bin);
      _set->remove(bin_match);
      _size -= bin_match.items_domain.size();
      delete &bin_match;
   }
}

/*
 * Remove all instances of the element(s) from the multimap.
 * Return a reference to the multimap.
 */
template <typename T, typename U, typename Syn>
multimap<T,U,Syn>& multimap<T,U,Syn>::remove_all(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   this->remove_item_all(t);
   return *this;
}

template <typename T, typename U, typename Syn>
multimap<T,U,Syn>& multimap<T,U,Syn>::remove_all(const collection<T>& c) {
   list<T> lst(c);
   auto_write_lock<const Syn> wlock(*this);
   for (typename list<T>::iterator_t i(lst); i.has_next(); )
      this->remove_item_all(i.next());
   return *this;
}

/*
 * Clear the multimap, resetting it to the empty multimap.
 * Return a reference to the multimap.
 */
template <typename T, typename U, typename Syn>
multimap<T,U,Syn>& multimap<T,U,Syn>::clear() {
   auto_write_lock<const Syn> wlock(*this);
   _set.reset(new set<equiv_bin>(_f_compare));
   _size = 0;
   return *this;
}

/*
 * Search.
 * Return an element in the domain matching the given element.
 * Throw an exception (ex_not_found) when attempting to find an element not 
 * contained in the map.
 */
template <typename T, typename U, typename Syn>
bool multimap<T,U,Syn>::contains(const T& t) const {
   equiv_bin bin(const_cast<T&>(t));   /* const_cast is safe here */
   auto_read_lock<const Syn> rlock(*this);
   return _set->contains(bin);
}

template <typename T, typename U, typename Syn>
T& multimap<T,U,Syn>::find(const T& t) const {
   equiv_bin bin(const_cast<T&>(t));   /* const_cast is safe here */
   auto_read_lock<const Syn> rlock(*this);
   equiv_bin& bin_match = _set->find(bin);
   return bin_match.items_domain.head();
}

/*
 * Search.
 * Return a list of all element(s) in the domain matching the given element.
 */
template <typename T, typename U, typename Syn>
list<T> multimap<T,U,Syn>::find_all(const T& t) const {
   equiv_bin bin(const_cast<T&>(t));   /* const_cast is safe here */
   auto_read_lock<const Syn> rlock(*this);
   if (_set->contains(bin)) {
      equiv_bin& bin_match = _set->find(bin);
      return bin_match.items_domain;
   } else {
      return list<T>();
   }
}

/*
 * Search.
 * Return an image of the given domain element under the multimap.
 * Throw an exception (ex_not_found) when attempting to find the image of an 
 * element not contained in the multimap.
 */
template <typename T, typename U, typename Syn>
U& multimap<T,U,Syn>::find_image(const T& t) const {
   equiv_bin bin(const_cast<T&>(t));   /* const_cast is safe here */
   auto_read_lock<const Syn> rlock(*this);
   equiv_bin& bin_match = _set->find(bin);
   return bin_match.items_range.head();
}

/*
 * Search.
 * Return a list of all element(s) in the range which are images of the 
 * given domain element.
 */
template <typename T, typename U, typename Syn>
list<U> multimap<T,U,Syn>::find_images(const T& t) const {
   equiv_bin bin(const_cast<T&>(t));   /* const_cast is safe here */
   auto_read_lock<const Syn> rlock(*this);
   if (_set->contains(bin)) {
      equiv_bin& bin_match = _set->find(bin);
      return bin_match.items_range;
   } else {
      return list<U>();
   }
}

/*
 * Size.
 * Return the number of mappings.
 */
template <typename T, typename U, typename Syn>
unsigned long multimap<T,U,Syn>::size() const {
   auto_read_lock<const Syn> rlock(*this);
   return _size;
}

/*
 * Return a list of elements in the domain of the multimap and a list of 
 * their corresponding images in the range.  All mappings t -> u in the 
 * multimap appear as a pair of elements (t, u) at corresponding positions
 * in the lists.  Hence, both the domain and range list may have repeated 
 * elements.
 */
template <typename T, typename U, typename Syn>
list<T> multimap<T,U,Syn>::domain() const {
   list<T> lst;
   auto_read_lock<const Syn> rlock(*this);
   typename set<equiv_bin>::iterator_t i(*_set);
   while (i.has_next())
      lst.add(i.next().items_domain);
   return lst;
}

template <typename T, typename U, typename Syn>
list<U> multimap<T,U,Syn>::range() const {
   list<U> lst;
   auto_read_lock<const Syn> rlock(*this);
   typename set<equiv_bin>::iterator_t i(*_set);
   while (i.has_next())
      lst.add(i.next().items_range);
   return lst;
}

/*
 * Return lists of elements in the domain and range.
 * This method allows one to atomically capture the contents of a multimap
 * if there is a potential for the multimap to change between calls of the
 * above domain() and range() methods.
 */
template <typename T, typename U, typename Syn>
void multimap<T,U,Syn>::contents(
   auto_ptr< collections::list<T> >& domain,
   auto_ptr< collections::list<U> >& range) const
{
   domain.reset(new list<T>());
   range.reset(new list<U>());
   auto_read_lock<const Syn> rlock(*this);
   typename set<equiv_bin>::iterator_t i(*_set);
   while (i.has_next()) {
      equiv_bin& bin = i.next();
      domain->add(bin.items_domain);
      range->add(bin.items_range);
   }
}

/*
 * Return iterators over elements in the domain.
 */
template <typename T, typename U, typename Syn>
auto_ptr< iterator<T> > multimap<T,U,Syn>::iter_create() const {
   return auto_ptr< iterator<T> >(
      new multimap_iterator<T,U,Syn>(*this)
   );
}

template <typename T, typename U, typename Syn>
auto_ptr< iterator<T> > multimap<T,U,Syn>::iter_reverse_create() const {
   return auto_ptr< iterator<T> >(
      new multimap_iterator_reverse<T,U,Syn>(*this)
   );
}

} /* namespace collections */

#endif
