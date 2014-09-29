/*
 * Maps (thread-safe).
 *
 * Maps do not allow duplicate items in the domain.  Adding a mapping for an
 * item already in the domain overwrites the old mapping for that item.
 */
#ifndef COLLECTIONS__MAP_HH
#define COLLECTIONS__MAP_HH

#include "collections/abstract/collection.hh"
#include "collections/abstract/map.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "collections/set.hh"
#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "functors/comparable_functors.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
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
using functors::comparable_functor;
using functors::compare_functors;
using lang::exceptions::ex_invalid_argument;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

/*
 * Declare classes for iterators over elements in the domain of the map.
 */
template <typename T, typename U, typename Syn>
class map_iterator;

template <typename T, typename U, typename Syn>
class map_iterator_reverse;

/*
 * Maps.
 */
template <typename T, typename U, typename Syn = unsynchronized>
class map : public abstract::map<T,U>,
            protected Syn {
public:
   /*
    * Friend classes.
    */
   friend class map_iterator<T,U,Syn>;
   friend class map_iterator_reverse<T,U,Syn>;
   
   /*
    * Define the iterator types.
    */
   typedef map_iterator<T,U,Syn> iterator_t;
   typedef map_iterator_reverse<T,U,Syn> iterator_reverse_t;

   /*
    * Constructor.
    * Return an empty map.
    * Optionally specify the comparison function to use for domain elements.
    */
   map(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Constructor.
    * Create a map containing the mappings t -> u for corresponding elements 
    * of the given collections.  Optionally specify the comparison function 
    * to use for domain elements.
    */
   explicit map(
      const collection<T>&, 
      const collection<U>&,
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Constructor.
    * Create a map with the same contents as the given map.
    * Optionally specify the comparison function to use for domain elements.
    */
   explicit map(
      const abstract::map<T,U>&,
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Copy constructor.
    */
   map(const map<T,U,Syn>&);

   /*
    * Destructor.
    */
   virtual ~map();
   
   /*
    * Add the mapping t -> u, overwriting any existing mapping for t.
    * Return a reference to the map.
    */
   map<T,U,Syn>& add(T& /* t */, U& /* u */);

   /*
    * Add the mapping t -> u for corresponding elements of the collections.
    * Overwriting previous mappings if they exist.
    * The collections must have the same size.
    * Return a reference to the map.
    */
   map<T,U,Syn>& add(const collection<T>&, const collection<U>&);

   /*
    * Add all mappings in the given map to the current map.
    * Overwriting previous mappings if they exist.
    * Return a reference to the map.
    */
   map<T,U,Syn>& add(const abstract::map<T,U>&);
   
   /*
    * Remove element(s) and their corresponding images from the map.
    * Return a reference to the map.
    */
   map<T,U,Syn>& remove(const T&);
   map<T,U,Syn>& remove(const collection<T>&);

   /*
    * Clear the map, resetting it to the empty map.
    * Return a reference to the map.
    */
   map<T,U,Syn>& clear();

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
    * This method allows one to atomically capture the contents of a map if 
    * there is a potential for the map to change between calls of the above 
    * domain() and range() methods.
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
    * Map data structures.
    ************************************************************************/

   /*
    * Image of an element under a map.
    */
   class map_image {
   public:
      /*
       * Constructor.
       */
      explicit map_image(U& item) : u(item) { }

      /*
       * Copy constructor.
       */
      map_image(const map_image& m) : u(m.u) { }

      /*
       * Destructor.
       */
      ~map_image() { /* do nothing */ }
     
      /*
       * Data.
       */
      U& u;
   };
   
   /*
    * Mapping between elements.
    */
   class mapping {
   public:
      /*
       * Constructor.
       * Create a t -> NULL mapping.
       */
      explicit mapping(T& t_item)
       : t(t_item),
         img(NULL)
      { }
      
      /*
       * Constructor.
       * Create a t -> u mapping.
       */
      explicit mapping(T& t_item, U& u_item)
       : t(t_item),
         img(new map_image(u_item))
      { }

      /*
       * Copy constructor.
       */
      mapping(const mapping& m)
       : t(m.t),
         img((m.img.get() != NULL) ? (new map_image(*(m.img))) : NULL)
      { }

      /*
       * Destructor.
       */
      ~mapping() { /* do nothing */ }

      /*
       * Data.
       */
      T& t;
      auto_ptr<map_image> img;
   };

   /*
    * Comparison functor on mappings.
    * Call the given compare functor on the corresponding domain elements.
    */
   class mapping_compare_functor : public comparable_functor<mapping> {
   public:
      /*
       * Constructor.
       */
      explicit mapping_compare_functor(const comparable_functor<T>& f)
       : _f(f)
      { }
      
      /*
       * Copy constructor.
       */
      explicit mapping_compare_functor(const mapping_compare_functor& f)
       : _f(f._f)
      { }

      /*
       * Comparison function.
       */
      int operator()(const mapping& m0, const mapping& m1) const {
         return _f(m0.t, m1.t);
      }

   protected:
      const comparable_functor<T>& _f;
   };

   /*
    * Map data.
    */
   const mapping_compare_functor            _f_compare; /* comparison functor */
   auto_collection< mapping, set<mapping> > _set;       /* set of mappings */

   /************************************************************************
    * Map helper functions.
    ************************************************************************/
    
   /*
    * Add a mapping to the map.
    */
   void add_mapping(T&, U&);

   /*
    * Remove an element and its image from the map.
    */
   void remove_item(const T&);
};

/*
 * Map iterator.
 */
template <typename T, typename U, typename Syn = unsynchronized>
class map_iterator : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit map_iterator(const map<T,U,Syn>&);

   /*
    * Copy constructor.
    */
   map_iterator(const map_iterator<T,U,Syn>&);

   /*
    * Destructor.
    */
   virtual ~map_iterator();

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
   const map<T,U,Syn>&                          _m;         /* map being iterated over */
   auto_read_lock<const Syn>                    _rlock;     /* read lock on map */
   set_iterator<typename map<T,U,Syn>::mapping> _map_iter;  /* iterator over mappings */
};

/*
 * Map reverse iterator.
 */
template <typename T, typename U, typename Syn = unsynchronized>
class map_iterator_reverse : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit map_iterator_reverse(const map<T,U,Syn>&);

   /*
    * Copy constructor.
    */
   map_iterator_reverse(const map_iterator_reverse<T,U,Syn>&);

   /*
    * Destructor.
    */
   virtual ~map_iterator_reverse();

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
   const map<T,U,Syn>&                                  _m;          /* map being iterated over */
   auto_read_lock<const Syn>                            _rlock;      /* read lock on map */
   set_iterator_reverse<typename map<T,U,Syn>::mapping> _map_iter;   /* iterator over mappings */
};

/***************************************************************************
 * Map iterator implementation.
 ***************************************************************************/

/*
 * Constructors.
 * Initialize iterator over mappings.
 */
template <typename T, typename U, typename Syn>
map_iterator<T,U,Syn>::map_iterator(const map<T,U,Syn>& m)
 : _m(m),
   _rlock(_m),
   _map_iter(*(_m._set))
{ }

template <typename T, typename U, typename Syn>
map_iterator_reverse<T,U,Syn>::map_iterator_reverse(const map<T,U,Syn>& m)
 : _m(m),
   _rlock(_m),
   _map_iter(*(_m._set))
{ }

/*
 * Copy constructors.
 */
template <typename T, typename U, typename Syn>
map_iterator<T,U,Syn>::map_iterator(
   const map_iterator<T,U,Syn>& i)
 : _m(i._m),
   _rlock(_m),
   _map_iter(i._map_iter)
{ }

template <typename T, typename U, typename Syn>
map_iterator_reverse<T,U,Syn>::map_iterator_reverse(
   const map_iterator_reverse<T,U,Syn>& i)
 : _m(i._m),
   _rlock(_m),
   _map_iter(i._map_iter)
{ }

/*
 * Destructors.
 * Do nothing as the map is automatically unlocked upon destruction of the 
 * read lock.
 */
template <typename T, typename U, typename Syn>
map_iterator<T,U,Syn>::~map_iterator() {
   /* do nothing */
}

template <typename T, typename U, typename Syn>
map_iterator_reverse<T,U,Syn>::~map_iterator_reverse() {
   /* do nothing */
}

/*
 * Check if there is another item available.
 */
template <typename T, typename U, typename Syn>
bool map_iterator<T,U,Syn>::has_next() const {
   return _map_iter.has_next();
}
   
template <typename T, typename U, typename Syn>
bool map_iterator_reverse<T,U,Syn>::has_next() const {
   return _map_iter.has_next();
}

/*
 * Return the next item.
 * Throw an exception (ex_not_found) if there are no more items.
 */
template <typename T, typename U, typename Syn>
T& map_iterator<T,U,Syn>::next() {
   return _map_iter.next().t;
}
   
template <typename T, typename U, typename Syn>
T& map_iterator_reverse<T,U,Syn>::next() {
   return _map_iter.next().t;
}

/***************************************************************************
 * Map implementation.
 ***************************************************************************/

/*
 * Constructor.
 * Return an empty map that uses the given comparison function for domain 
 * elements.
 */
template <typename T, typename U, typename Syn>
map<T,U,Syn>::map(const comparable_functor<T>& f)
 : Syn(),
   _f_compare(f),
   _set(new set<mapping>(_f_compare))
{ }

/*
 * Constructor.
 * Create a map containing the mappings t -> u for corresponding elements 
 * of the given collections.  Specify the comparison function to use for 
 * domain elements.
 */
template <typename T, typename U, typename Syn>
map<T,U,Syn>::map(
   const collection<T>& c_t, 
   const collection<U>& c_u,
   const comparable_functor<T>& f)
 : Syn(),
   _f_compare(f),
   _set(new set<mapping>(_f_compare))
{
   this->add(c_t, c_u);
}

/*
 * Constructor.
 * Create a map with the same contents as the given map.
 * Specify the comparison function to use for domain elements.
 */
template <typename T, typename U, typename Syn>
map<T,U,Syn>::map(
   const abstract::map<T,U>& m,
   const comparable_functor<T>& f)
 : Syn(),
   _f_compare(f),
   _set(new set<mapping>(_f_compare))
{
   this->add(m);
}

/*
 * Copy constructor.
 */
template <typename T, typename U, typename Syn>
map<T,U,Syn>::map(const map<T,U,Syn>& m)
 : Syn(),
   _f_compare(m._f_compare),
   _set(new set<mapping>(_f_compare))
{
   typename set<mapping>::iterator_t i(*(m._set));
   while (i.has_next()) {
      mapping& mppng = i.next();
      this->add_mapping(mppng.t, mppng.img->u);
   }
}

/*
 * Destructor.
 */
template <typename T, typename U, typename Syn>
map<T,U,Syn>::~map() {
   /* do nothing */
}

/*
 * Add a mapping to the map.
 */
template <typename T, typename U, typename Syn>
void map<T,U,Syn>::add_mapping(T& t, U& u) {
   auto_ptr<mapping> m(new mapping(t,u));
   if (_set->contains(*m)) {
      mapping& m_match = _set->find(*m);
      _set->remove(m_match);
      delete &m_match; 
   }
   _set->add(*m);
   m.release();
}

/*
 * Add the mapping t -> u, overwriting any existing mapping for t.
 * Return a reference to the map.
 */
template <typename T, typename U, typename Syn>
map<T,U,Syn>& map<T,U,Syn>::add(T& t, U& u) {
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
template <typename T, typename U, typename Syn>
map<T,U,Syn>& map<T,U,Syn>::add(
   const collection<T>& c_t,
   const collection<U>& c_u)
{
   list<T> lst_t(c_t);
   list<U> lst_u(c_u);
   if (c_t.size() != c_u.size())
      throw ex_invalid_argument(
         "collections being add to map must have the same size"
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
template <typename T, typename U, typename Syn>
map<T,U,Syn>& map<T,U,Syn>::add(const abstract::map<T,U>& m) {
   auto_ptr< list<T> > domain;
   auto_ptr< list<U> > range;
   m.contents(domain, range);
   return this->add(*domain, *range);
}

/*
 * Remove an element and its image from the map.
 */
template <typename T, typename U, typename Syn>
void map<T,U,Syn>::remove_item(const T& t) {
   mapping m(const_cast<T&>(t)); /* const_cast is safe here */
   if (_set->contains(m)) {
      mapping& m_match = _set->find(m);
      _set->remove(m_match);
      delete &m_match; 
   }
}

/*
 * Remove element(s) and their corresponding images from the map.
 * Return a reference to the map.
 */
template <typename T, typename U, typename Syn>
map<T,U,Syn>& map<T,U,Syn>::remove(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   this->remove_item(t);
   return *this;
}

template <typename T, typename U, typename Syn>
map<T,U,Syn>& map<T,U,Syn>::remove(const collection<T>& c_t) {
   list<T> lst(c_t);
   auto_write_lock<const Syn> wlock(*this);
   for (typename list<T>::iterator_t i(lst); i.has_next(); )
      this->remove_item(i.next());
   return *this;
}

/*
 * Clear the map, resetting it to the empty map.
 * Return a reference to the map.
 */
template <typename T, typename U, typename Syn>
map<T,U,Syn>& map<T,U,Syn>::clear() {
   auto_write_lock<const Syn> wlock(*this);
   _set.reset(new set<mapping>(_f_compare));
   return *this;
}

/*
 * Search.
 * Return the element in the domain matching the given element.
 * Throw an exception (ex_not_found) when attempting to find an element not 
 * contained in the map.
 */
template <typename T, typename U, typename Syn>
bool map<T,U,Syn>::contains(const T& t) const {
   mapping m(const_cast<T&>(t)); /* const_cast is safe here */
   auto_read_lock<const Syn> rlock(*this);
   return _set->contains(m);
}

template <typename T, typename U, typename Syn>
T& map<T,U,Syn>::find(const T& t) const {
   mapping m(const_cast<T&>(t)); /* const_cast is safe here */
   auto_read_lock<const Syn> rlock(*this);
   mapping& m_match = _set->find(m);
   return m_match.t;
}

/*
 * Search.
 * Return the image of the given domain element under the map.
 * Throw an exception (ex_not_found) when attempting to find the image of an
 * element not contained in the map.
 */
template <typename T, typename U, typename Syn>
U& map<T,U,Syn>::find_image(const T& t) const {
   mapping m(const_cast<T&>(t)); /* const_cast is safe here */
   auto_read_lock<const Syn> rlock(*this);
   mapping& m_match = _set->find(m);
   return m_match.img->u;
}

/*
 * Size.
 * Return the number of elements in the domain.
 */
template <typename T, typename U, typename Syn>
unsigned long map<T,U,Syn>::size() const {
   auto_read_lock<const Syn> rlock(*this);
   return _set->size();
}

/*
 * Return a list of elements in the domain of the map and a list of their 
 * corresponding images in the range.  Hence, range() does not necessarily 
 * return a unique list of elements as two or more elements of the domain 
 * might have the same image.
 */
template <typename T, typename U, typename Syn>
list<T> map<T,U,Syn>::domain() const {
   list<T> lst;
   auto_read_lock<const Syn> rlock(*this);
   typename set<mapping>::iterator_t i(*_set);
   while (i.has_next())
      lst.add(i.next().t);
   return lst;
}

template <typename T, typename U, typename Syn>
list<U> map<T,U,Syn>::range() const {
   list<U> lst;
   auto_read_lock<const Syn> rlock(*this);
   typename set<mapping>::iterator_t i(*_set);
   while (i.has_next())
      lst.add(i.next().img->u);
   return lst;
}

/*
 * Return lists of elements in the domain and range.
 */
template <typename T, typename U, typename Syn>
void map<T,U,Syn>::contents(
   auto_ptr< collections::list<T> >& domain,
   auto_ptr< collections::list<U> >& range) const
{
   domain.reset(new list<T>());
   range.reset(new list<U>());
   auto_read_lock<const Syn> rlock(*this);
   typename set<mapping>::iterator_t i(*_set);
   while (i.has_next()) {
      mapping& m = i.next();
      domain->add(m.t);
      range->add(m.img->u);
   }
}

/*
 * Return iterators over elements in the domain.
 */
template <typename T, typename U, typename Syn>
auto_ptr< iterator<T> > map<T,U,Syn>::iter_create() const {
   return auto_ptr< iterator<T> >(new map_iterator<T,U,Syn>(*this));
}

template <typename T, typename U, typename Syn>
auto_ptr< iterator<T> > map<T,U,Syn>::iter_reverse_create() const {
   return auto_ptr< iterator<T> >(new map_iterator_reverse<T,U,Syn>(*this));
}

} /* namespace collections */

#endif
