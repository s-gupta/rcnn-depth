/*
 * Array lists (thread-safe).
 * An array list is a list stored as an array.
 */
#ifndef COLLECTIONS__ARRAY_LIST_HH
#define COLLECTIONS__ARRAY_LIST_HH

#include "collections/abstract/collection.hh"
#include "collections/abstract/array.hh"
#include "collections/abstract/list.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "config/safety.hh"
#include "functors/comparable_functors.hh"
#include "io/serialization/serial_input_stream.hh"
#include "io/serialization/serial_output_stream.hh"
#include "io/serialization/serializers.hh"
#include "lang/array.hh"
#include "lang/exceptions/ex_index_out_of_bounds.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "lang/exceptions/ex_not_found.hh"
#include "lang/iterators/iterator.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/math.hh"
#include "math/random/generators/rand_gen.hh"

/* 
 * Enable/disable bounds checking for array lists.
 */
#ifdef CONFIG__SAFETY__CHECK_BOUNDS
   #define COLLECTIONS__ARRAY_LIST__CHECK_BOUNDS (true)
#else
   #define COLLECTIONS__ARRAY_LIST__CHECK_BOUNDS (false)
#endif

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
using io::serialization::serial_input_stream;
using io::serialization::serial_output_stream;
using io::serialization::serializer;
using io::serialization::serializers;
using lang::exceptions::ex_index_out_of_bounds;
using lang::exceptions::ex_invalid_argument;
using lang::exceptions::ex_not_found;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;
using math::random::generators::rand_gen;

/*
 * Declare classes for iterators over array lists.
 */
template <typename T, typename Syn>
class array_list_iterator;

template <typename T, typename Syn>
class array_list_iterator_reverse;

/*
 * Array lists.
 */
template <typename T, typename Syn = unsynchronized>
class array_list : public abstract::array<T>, 
                   public abstract::list<T>,
                   protected Syn {
public:
   /*
    * Friend classes.
    */
   friend class array_list_iterator<T,Syn>;
   friend class array_list_iterator_reverse<T,Syn>;

   /*
    * Define the iterator types.
    */
   typedef array_list_iterator<T,Syn> iterator_t;
   typedef array_list_iterator_reverse<T,Syn> iterator_reverse_t;

   /*
    * Constructor.
    * Create an empty array list.
    */
   array_list();

   /*
    * Constructor.
    * Create an array list from a collection.
    */
   explicit array_list(const collection<T>&);
   
   /*
    * Copy constructor.
    */
   array_list(const array_list<T,Syn>&);

   /*
    * Destructor.
    */
   virtual ~array_list();

   /*
    * Serialize.
    */
   void serialize(
      serial_output_stream&,
      const serializer<T>& = serializers<T>::s_default()
   ) const;

   /*
    * Deserialize.
    */
   static auto_collection< T, array_list<T,Syn> > deserialize(
      serial_input_stream&,
      const serializer<T>& = serializers<T>::s_default()
   );

   /*
    * Element reference.
    */
   T& operator[](unsigned long) const;
   
   /*
    * Element replacement.
    * Return a reference to the element that was replaced.
    */
   T& replace(unsigned long, T&);
   
   /*
    * Subarray operations.
    * Add elements in the specified range to the given collection.
    */
   void subarray(
      unsigned long,                      /* start index */
      unsigned long,                      /* end index   */
      collection<T>&
   ) const;
   
   void subarray(
      unsigned long,                      /* start index */ 
      unsigned long,                      /* step size   */
      unsigned long,                      /* end index   */
      collection<T>&
   ) const;
   
   void subarray(
      const lang::array<unsigned long>&,  /* indices of desired elements */
      collection<T>&
   ) const;

   /*
    * Return head element(s).
    * Throw an exception (ex_not_found) if the array list doesn't contain
    * enough elements.
    */
   T& head() const;
   void head(unsigned long, collection<T>&) const;

   /*
    * Return tail element(s). 
    * Throw an exception (ex_not_found) if the array list doesn't contain
    * enough elements.
    */
   T& tail() const;
   void tail(unsigned long, collection<T>&) const;

   /*
    * Collection interface for adding element(s) to tail of array list.
    * Return a reference to the array list.
    */
   array_list<T,Syn>& add(T&);
   array_list<T,Syn>& add(const collection<T>&);
    
   /*
    * Addition of element(s) to head of array list.
    * Return a reference to the array list.
    */
   array_list<T,Syn>& prepend(T&);
   array_list<T,Syn>& prepend(const collection<T>&);
   
   /*
    * Addition of element(s) to tail of array list.
    * Return a reference to the array list.
    */
   array_list<T,Syn>& append(T&);
   array_list<T,Syn>& append(const collection<T>&);
   
   /*
    * Remove and return head element(s).
    * Throw an exception (ex_not_found) if the array list doesn't contain
    * enough elements.
    */
   T& remove_head();
   void remove_head(unsigned long, collection<T>&);

   /*
    * Remove and return tail element(s).
    * Throw an exception (ex_not_found) if the array list doesn't contain
    * enough elements.
    */
   T& remove_tail();
   void remove_tail(unsigned long, collection<T>&);

   /*
    * Remove all element(s) from the array list.
    * Return a reference to the array list.
    */
   array_list<T,Syn>& clear();

   /*
    * Reverse array list.
    * Return a reference to the array list.
    */
   array_list<T,Syn>& reverse();

   /*
    * Size.
    */
   unsigned long size() const;

   /*
    * Get/set expand factor.
    * Allocated array storage automatically grows/shrinks by this factor as
    * elements are added/removed.
    */
   double expand_factor() const;
   double expand_factor(double);
   
   /*
    * Return iterators over elements.
    */
   auto_ptr< iterator<T> > iter_create() const;
   auto_ptr< iterator<T> > iter_reverse_create() const;

   /*
    * Randomly permute the elements of the array list.
    * Return an index array mapping resulting position --> original position.
    *
    * The generator from which to draw values used in computing the random
    * permutation may optionally be specified.  If unspecified, a generator 
    * is created and used for this purpose.
    */
   using abstract::array<T>::randperm;
   
   lang::array<unsigned long> randperm(rand_gen<>&);
   
   /*
    * Randomly permute the elements of the specified subarray within the array.
    * Return an index array mapping resulting position --> original position
    * (within the subarray, so index 0 corresponds to start of subarray).
    */
   using abstract::array<T>::randperm_subarray;

   lang::array<unsigned long> randperm_subarray(
      unsigned long,    /* start index */ 
      unsigned long,    /* end index   */
      rand_gen<>&       /* generator   */
   );

   /*
    * Sort elements in ascending order according to the given comparison
    * functor.
    *
    * Sorting uses an O(n + (n/p)*log(n/p)) expected time algorithm, where
    * n is the array size and p is the number of available processors.
    */
   void sort(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Sort the array list.
    * Return an index array mapping sorted position --> original position.
    */
   lang::array<unsigned long> sort_idx(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );
  
   /*
    * Sort the specified subarray within the array list.
    */
   void sort_subarray(
      unsigned long,    /* start index */ 
      unsigned long,    /* end index   */
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Sort the specified subarray within the array list.
    * Return an index array mapping sorted position --> original position
    * (within the subarray, so index 0 corresponds to start of subarray).
    */
   lang::array<unsigned long> sort_subarray_idx(
      unsigned long,    /* start index */ 
      unsigned long,    /* end index   */
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Remove duplicate elements (as defined by the given comparison functor)
    * and place the remaining unique elements in sorted order.
    */
   void unique(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Remove duplicate elements (as defined by the given comparison functor)
    * and place the remaining unique elements in sorted order.  In addition,
    * return an index array containing the original positions of the unique
    * elements.
    */
   lang::array<unsigned long> unique_idx(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

protected:
   /************************************************************************
    * Array list data structures.
    ************************************************************************/
    
   /*
    * Node in the array list.
    */
   class array_list_node {
   public:
      /*
       * Constructors.
       */
      explicit array_list_node(T& item) : t(item) { }
      explicit array_list_node(const array_list_node& node) : t(node.t) { }

      /*
       * Destructor.
       */
      ~array_list_node() { /* do nothing */ }

      /*
       * Node data.
       */
      T& t;
   };

   /*
    * Array list data.
    */
   unsigned long _size;                            /* number of elements in the array */
   unsigned long _head_space;                      /* extra space allocated at head of array */
   unsigned long _tail_space;                      /* extra space allocated at tail of array */
   lang::array< auto_ptr<array_list_node> > _data; /* array elements */

   /*
    * Parameters for managing extra space.
    */
   unsigned long _empty_space_max;  /* max extra head/tail space before shrinking array */
   double _expand_factor;           /* factor by which we grow/shrink array */
   double _empty_factor;            /* factor determining amount of empty head/tail space allowed */
                                    /* (computed as (_expand_factor + 1 - 1/(1+_expand_factor)) */
                                 
   /*
    * Default expand factor.
    */
   static const double default_expand_factor;
  
   /************************************************************************
    * Array list node compare functor.
    ************************************************************************/

   /*
    * Comparison functor on pointers to list nodes.
    */
   class ptr_node_compare_functor
    : public comparable_functor< auto_ptr<array_list_node> > {
   public:
      /*
       * Constructor.
       */
      explicit ptr_node_compare_functor(const comparable_functor<T>& f)
       : _f(f) { }

      /*
       * Copy constructor.
       */
      explicit ptr_node_compare_functor(const ptr_node_compare_functor& f)
       : _f(f._f) { }
      
      /*
       * Comparison function.
       */
      int operator()(
         const auto_ptr<array_list_node>& p0,
         const auto_ptr<array_list_node>& p1) const
      { return _f(p0->t, p1->t); }
      
   protected:
      const comparable_functor<T>& _f;
   };

   /************************************************************************
    * Array list helper functions.
    ************************************************************************/
   
   /*
    * Resize the array list to have the given amount of head and tail space.
    */
   void resize(unsigned long, unsigned long);
   
   /*
    * Expand the head of the array list by the given number of elements.
    * Allocate more space and copy array as needed.
    */
   void resize_expand_head(unsigned long);

   /*
    * Expand the tail of the array list by the given number of elements.
    * Allocate more space and copy array as needed.
    */
   void resize_expand_tail(unsigned long);

   /*
    * Shrink the head of the array list by the given number of elements.
    * Reallocate less space and copy array as needed.
    */
   void resize_shrink_head(unsigned long);

   /*
    * Shrink the tail of the array list by the given number of elements.
    * Reallocate less space and copy array as needed.
    */
   void resize_shrink_tail(unsigned long);
};

/*
 * Array list iterator.
 */
template <typename T, typename Syn = unsynchronized>
class array_list_iterator : public iterator<T> {
public:
   /*
    * Constructor.
    * Seek to the head of the array list.
    */
   explicit array_list_iterator(const array_list<T,Syn>&);

   /*
    * Copy constructor.
    */
   array_list_iterator(const array_list_iterator<T,Syn>&);

   /*
    * Destructor.
    */
   virtual ~array_list_iterator();

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
   const array_list<T,Syn>& _arr;   /* array list being iterated over */
   unsigned long _curr;             /* current position in array list */
};

/*
 * Array_list reverse iterator.
 */
template <typename T, typename Syn = unsynchronized>
class array_list_iterator_reverse : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit array_list_iterator_reverse(const array_list<T,Syn>&);

   /*
    * Copy constructor.
    */
   array_list_iterator_reverse(const array_list_iterator_reverse<T,Syn>&);
   
   /*
    * Destructor.
    */
   virtual ~array_list_iterator_reverse();
   
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
   const array_list<T,Syn>& _arr;   /* array list being iterated over */
   unsigned long _curr;             /* current position in array list */
};

/***************************************************************************
 * Array list iterator implementation.
 ***************************************************************************/

/*
 * Constructors.
 * Lock array list and initialize position.
 */
template <typename T, typename Syn>
array_list_iterator<T,Syn>::array_list_iterator(
   const array_list<T,Syn>& a)
 : _arr(a)
{
   _arr.read_lock();
   _curr = 0;
}

template <typename T, typename Syn>
array_list_iterator_reverse<T,Syn>::array_list_iterator_reverse(
   const array_list<T,Syn>& a)
 : _arr(a)
{
   _arr.read_lock();
   _curr = _arr._size;
}

/*
 * Copy constructors.
 */
template <typename T, typename Syn>
array_list_iterator<T,Syn>::array_list_iterator(
   const array_list_iterator<T,Syn>& i) 
 : _arr(i._arr),
   _curr(i._curr)
{
   _arr.read_lock();
}

template <typename T, typename Syn>
array_list_iterator_reverse<T,Syn>::array_list_iterator_reverse(
   const array_list_iterator_reverse<T,Syn>& i)
 : _arr(i._arr),
   _curr(i._curr)
{
   _arr.read_lock();
}

/*
 * Destructors.
 * Unlock the array list.
 */
template <typename T, typename Syn>
array_list_iterator<T,Syn>::~array_list_iterator() {
   _arr.read_unlock();
}

template <typename T, typename Syn>
array_list_iterator_reverse<T,Syn>::~array_list_iterator_reverse() {
   _arr.read_unlock();
}

/*
 * Check if there is another item available.
 */
template <typename T, typename Syn>
bool array_list_iterator<T,Syn>::has_next() const {
   return (_curr < _arr._size);
}

template <typename T, typename Syn>
bool array_list_iterator_reverse<T,Syn>::has_next() const {
   return (_curr > 0);
}

/*
 * Return the next item.
 * Throw an exception (ex_not_found) if there are no more items.
 */
template <typename T, typename Syn>
T& array_list_iterator<T,Syn>::next() {
   if (_curr < _arr._size) {
      T& t = _arr._data[_arr._head_space + _curr]->t;
      _curr++;
      return t;
   } else {
      throw ex_not_found(
         "iterator attempted to read past tail of array_list"
      );
   }
}

template <typename T, typename Syn>
T& array_list_iterator_reverse<T,Syn>::next() {
   if (_curr > 0) {
      _curr--;
      return _arr._data[_arr._head_space + _curr]->t;
   } else {
      throw ex_not_found(
         "reverse iterator attempted to read past head of array_list"
      );
   }
}

/***************************************************************************
 * Array list implementation.
 ***************************************************************************/

/*
 * Default expand factor.
 */
template <typename T, typename Syn>
const double array_list<T,Syn>::default_expand_factor = 0.2;

/*
 * Default constructor.
 * Create an empty list.
 */
template <typename T, typename Syn>
array_list<T,Syn>::array_list()
 : Syn(),
   _size(0),
   _head_space(0),
   _tail_space(0),
   _data(_head_space + _size + _tail_space),
   _empty_space_max(0),
   _expand_factor(array_list<T,Syn>::default_expand_factor), 
   _empty_factor(
      _expand_factor + double(1) - (double(1)/(double(1) + _expand_factor))
   )
{ }

/*
 * Constructor.
 * Create an array list from a collection.
 */
template <typename T, typename Syn>
array_list<T,Syn>::array_list(const collection<T>& c)
 : Syn(),
   _size(0),
   _head_space(0),
   _tail_space(0),
   _data(_head_space + _size + _tail_space),
   _empty_space_max(0),
   _expand_factor(array_list<T,Syn>::default_expand_factor), 
   _empty_factor(
      _expand_factor + double(1) - (double(1)/(double(1) + _expand_factor))
   ) 
{
   /* add elements of collection */
   this->add(c);
}

/*
 * Copy constructor.
 */
template <typename T, typename Syn>
array_list<T,Syn>::array_list(const array_list<T,Syn>& a)
 : Syn(),
   _size(0),
   _head_space(0),
   _tail_space(0),
   _data(_head_space + _size + _tail_space),
   _empty_space_max(0),
   _expand_factor(array_list<T,Syn>::default_expand_factor), 
   _empty_factor(
      _expand_factor + double(1) - (double(1)/(double(1) + _expand_factor))
   )
{
   auto_read_lock<const Syn> rlock(a);
   /* copy size information */
   _size            = a._size;
   _head_space      = a._head_space;
   _tail_space      = a._tail_space;
   /* copy extra space parameters */
   _empty_space_max = a._empty_space_max;
   _expand_factor   = a._expand_factor;
   _empty_factor    = a._empty_factor;
   /* resize array */
   _data.resize(_head_space + _size + _tail_space);
   /* copy elements */
   for (unsigned long n = _head_space; n < (_size + _head_space); n++)
      _data[n].reset(new array_list_node(a._data[n]->t));
}

/*
 * Destructor.
 */
template <typename T, typename Syn>
array_list<T,Syn>::~array_list() {
   /* do nothing */
}

/*
 * Serialize.
 */
template <typename T, typename Syn>
void array_list<T,Syn>::serialize(
   serial_output_stream& s, const serializer<T>& slzr) const
{
   auto_read_lock<const Syn> rlock(*this);
   /* serialize parameters */
   s << _expand_factor;
   s << _empty_factor;
   /* serialize data */
   s << _size;
   for (unsigned long n = _head_space; n < (_size + _head_space); n++)
      slzr.serialize(s, _data[n]->t);
}

/*
 * Deserialize.
 */
template <typename T, typename Syn>
auto_collection< T, array_list<T,Syn> > array_list<T,Syn>::deserialize(
   serial_input_stream& s, const serializer<T>& slzr)
{
   auto_collection< T, array_list<T,Syn> > a(new array_list<T,Syn>());
   /* deserialize parameters */
   s >> a->_expand_factor;
   s >> a->_empty_factor;
   /* deserialize data */
   s >> a->_size;
   a->_data.resize(a->_head_space + a->_size + a->_tail_space);
   for (unsigned long n = a->_head_space; n < (a->_size + a->_head_space); n++)
   {
      auto_ptr<T> t = slzr.deserialize(s);
      a->_data[n].reset(new array_list_node(*t));
      t.release();
   }
   return a;
}

/*
 * Element reference.
 */
template <typename T, typename Syn>
T& array_list<T,Syn>::operator[](unsigned long n) const {
   auto_read_lock<const Syn> rlock(*this);
   #ifdef COLLECTIONS__ARRAY_LIST__CHECK_BOUNDS
      /* perform bounds check */
      if (n >= _size)
         throw ex_index_out_of_bounds(
            "array_list index out of bounds", 
            n
         );
   #endif
   /* retrieve element */
   return _data[_head_space + n]->t;
}
 
/*
 * Element replacement.
 * Replace the element at the given index.
 * Return a reference to the element that was replaced.
 */
template <typename T, typename Syn>
T& array_list<T,Syn>::replace(unsigned long n, T& t) {
   auto_write_lock<const Syn> wlock(*this);
   #ifdef COLLECTIONS__ARRAY_LIST__CHECK_BOUNDS
      /* perform bounds check */
      if (n >= _size)
         throw ex_index_out_of_bounds(
            "array_list index out of bounds", 
            n
         );
   #endif
   /* replace element */
   T& t_old = _data[_head_space + n]->t;
   _data[_head_space + n].reset(new array_list_node(t));
   return t_old;
}

/*
 * Add subarray of elements in range [start:end] to given collection.
 */
template <typename T, typename Syn>
void array_list<T,Syn>::subarray(
   unsigned long start,
   unsigned long end,
   collection<T>& c) const 
{
   this->subarray(start, 1, end, c);
}

/*
 * Add subarray of elements in range [start:step:end] to the given collection.
 */
template <typename T, typename Syn>
void array_list<T,Syn>::subarray(
   unsigned long start,
   unsigned long step,
   unsigned long end,
   collection<T>& c) const 
{
   list<T> l;  /* list of elements in range */
   {
      /* lock array list and obtain elements in range */
      auto_read_lock<const Syn> rlock(*this);
      #ifdef COLLECTIONS__ARRAY_LIST__CHECK_BOUNDS
         /* check array bounds */
         if (start >= _size)
            throw ex_index_out_of_bounds(
               "array_list start index out of bounds", 
               start
            );
         else if (end >= _size)
            throw ex_index_out_of_bounds(
               "array_list end index out of bounds", 
               end
            );
      #endif
      /* check step size */
      if (step == 0)
         throw ex_invalid_argument(
            "array_list step must be nonzero"
         );
      /* grab subarray */
      for (unsigned long n = start; n <= end; n += step)
         l.add(_data[_head_space + n]->t);
   }
   c.add(l);
}

/*
 * Add the subarray with elements taken from the given indices to the
 * given collection.
 */
template <typename T, typename Syn>
void array_list<T,Syn>::subarray(
   const lang::array<unsigned long>& index_arr, 
   collection<T>& c) const 
{
   list<T> l;  /* list of elements at given indices */
   {
      /* lock array list and obtain elements in range */
      auto_read_lock<const Syn> rlock(*this);
      unsigned long n_indices = index_arr.size();
      for (unsigned long n = 0; n < n_indices; n++) {
         /* get index */
         unsigned long index = index_arr[n];
         #ifdef COLLECTIONS__ARRAY_LIST__CHECK_BOUNDS
            /* check index bounds */
            if (index >= _size)
               throw ex_index_out_of_bounds(
                  "array_list index out of bounds", 
                  index
               );
         #endif
         /* add the element */
         l.add(_data[_head_space + index]->t);
      }
   }
   c.add(l);
}

/*
 * Return head element.
 * Throw an exception (ex_not_found) if the array list is empty.
 */
template <typename T, typename Syn>
T& array_list<T,Syn>::head() const {
   auto_read_lock<const Syn> rlock(*this);
   if (_size == 0)
      throw ex_not_found(
         "attempt to fetch head element of empty array_list"
      );
   return _data[_head_space]->t;
}

/*
 * Add the first n elements of the array list to the given collection.
 * Throw an exception (ex_not_found) if the array list doesn't contain at 
 * least n elements.
 */
template <typename T, typename Syn>
void array_list<T,Syn>::head(unsigned long n, collection<T>& c) const {
   list<T> l;
   {
      /* lock array list and obtain first n elements */
      auto_read_lock<const Syn> rlock(*this);
      if (_size < n)
         throw ex_not_found(
            "attempt to fetch too many head elements of array_list"
         );
      for (unsigned long i = 0; i < n; i++)
         l.add(_data[_head_space + i]->t);
   }
   c.add(l);
}

/*
 * Return tail element.
 * Throw an exception (ex_not_found) if the array list is empty.
 */
template <typename T, typename Syn>
T& array_list<T,Syn>::tail() const {
   auto_read_lock<const Syn> rlock(*this);
   if (_size == 0)
      throw ex_not_found(
         "attempt to fetch tail element of empty array_list"
      );
   return _data[_head_space + _size-1]->t;
}

/*
 * Add the last n elements of the array list to the given collection.
 * Throw an exception (ex_not_found) if the array list doesn't contain at
 * least n elements.
 */
template <typename T, typename Syn>
void array_list<T,Syn>::tail(unsigned long n, collection<T>& c) const {
   list<T> l;
   {
      /* lock array list and obtain last n elements */
      auto_read_lock<const Syn> rlock(*this);
      if (_size < n)
         throw ex_not_found(
            "attempt to fetch too many tail elements of array_list"
         );
      for (unsigned long i = _size-n; i < _size; i++)
         l.add(_data[_head_space + i]->t);
   }
   c.add(l);
}

/*
 * Collection interface to add an element to the tail of the array list.
 * Return a reference to the array list.
 */
template <typename T, typename Syn>
array_list<T,Syn>& array_list<T,Syn>::add(T& t) {
   return this->append(t);
}

/*
 * Add all elements in a collection to the tail of the array list.
 * Return a reference to the array list.
 */
template <typename T, typename Syn>
array_list<T,Syn>& array_list<T,Syn>::add(const collection<T>& c) {
   return this->append(c);
}

/*
 * Add an element to the head of the array list.
 * Return a reference to the array list.
 */
template <typename T, typename Syn>
array_list<T,Syn>& array_list<T,Syn>::prepend(T& t) {
   auto_write_lock<const Syn> wlock(*this);
   this->resize_expand_head(1);
   _data[_head_space].reset(new array_list_node(t));
   return *this;
}

/*
 * Add multiple elements to the head of the array list.
 * Return a reference to the array list.
 */
template <typename T, typename Syn>
array_list<T,Syn>& array_list<T,Syn>::prepend(const collection<T>& c) {
   list<T> lst(c);
   auto_write_lock<const Syn> wlock(*this);
   this->resize_expand_head(lst.size());
   unsigned long n = 0;
   for (typename list<T>::iterator_t i(lst); i.has_next(); n++)
      _data[_head_space + n].reset(new array_list_node(i.next()));
   return *this;
}

/*
 * Add an element to the tail of the array list.
 * Return a reference to the array list.
 */
template <typename T, typename Syn>
array_list<T,Syn>& array_list<T,Syn>::append(T& t) {
   auto_write_lock<const Syn> wlock(*this);
   this->resize_expand_tail(1);
   _data[_head_space + _size-1].reset(new array_list_node(t));
   return *this;
}
 
/*
 * Add multiple elements to the tail of the array list.
 * Return a reference to the array list.
 */
template <typename T, typename Syn>
array_list<T,Syn>& array_list<T,Syn>::append(const collection<T>& c) {
   list<T> lst(c);
   auto_write_lock<const Syn> wlock(*this);
   unsigned long l_size = lst.size();
   this->resize_expand_tail(l_size);
   unsigned long n = (_size - l_size);
   for (typename list<T>::iterator_t i(lst); i.has_next(); n++)
      _data[_head_space + n].reset(new array_list_node(i.next()));
   return *this;
}

/*
 * Remove and return the head element of the array list.
 * Throw an exception (ex_not_found) if the array list is empty.
 */
template <typename T, typename Syn>
T& array_list<T,Syn>::remove_head() {
   auto_write_lock<const Syn> wlock(*this);
   if (_size == 0)
      throw ex_not_found(
         "attempt to remove head element of empty array_list"
      );
   T& t = _data[_head_space]->t;
   this->resize_shrink_head(1);
   return t;
}

/*
 * Remove and return the first n elements of the array list.
 * Throw an exception (ex_not_found) if the array list contains less
 * than n elements.
 */
template <typename T, typename Syn>
void array_list<T,Syn>::remove_head(unsigned long n, collection<T>& c) {
   list<T> l;
   {
      /* lock array list, obtain and remove first n elements */
      auto_write_lock<const Syn> wlock(*this);
      if (_size < n)
         throw ex_not_found(
            "attempt to remove too many head elements of array_list"
         );
      for (unsigned long i = 0; i < n; i++)
         l.add(_data[_head_space + i]->t);
      this->resize_shrink_head(n);
   }
   c.add(l);
}

/*
 * Remove and return the tail element of the array list.
 * Throw an exception (ex_not_found) if the array list is empty.
 */
template <typename T, typename Syn>
T& array_list<T,Syn>::remove_tail() {
   auto_write_lock<const Syn> wlock(*this);
   if (_size == 0)
      throw ex_not_found(
         "attempt to remove tail element of empty array_list"
      );
   T& t = _data[_head_space + _size-1]->t;
   this->resize_shrink_tail(1);
   return t;
}

/*
 * Remove and return the last n elements of the array list.
 * Throw an exception (ex_not_found) if the array list contains less
 * than n elements.
 */
template <typename T, typename Syn>
void array_list<T,Syn>::remove_tail(unsigned long n, collection<T>& c) {
   list<T> l;
   {
      /* lock array list, obtain and remove last n elements */
      auto_write_lock<const Syn> wlock(*this);
      if (_size < n)
         throw ex_not_found(
            "attempt to remove too many tail elements of array_list"
         );
      for (unsigned long i = _size-n; i < _size; i++)
         l.add(_data[_head_space + i]->t);
      this->resize_shrink_tail(n);
   }
   c.add(l);
}

/*
 * Remove all element(s) from the array list.
 * Return a reference to the array list.
 */
template <typename T, typename Syn>
array_list<T,Syn>& array_list<T,Syn>::clear() {
   auto_write_lock<const Syn> wlock(*this);
   _size = 0;
   _head_space = 0;
   _tail_space = 0;
   _data.resize(0);
   _empty_space_max = 0;
   return *this;
}

/*
 * Reverse the order of elements in the array list.
 * Return a reference to the array list.
 */
template <typename T, typename Syn>
array_list<T,Syn>& array_list<T,Syn>::reverse() {
   auto_write_lock<const Syn> wlock(*this);
   for (unsigned long start = _head_space, end = _head_space + _size;
        start < end;
        start++)
   {
      end--;
      auto_ptr<array_list_node> temp = _data[start];
      _data[start] = _data[end];
      _data[end]   = temp;
   }
   return *this;
}

/*
 * Get number of elements in array list.
 */
template <typename T, typename Syn>
unsigned long array_list<T,Syn>::size() const {
   auto_read_lock<const Syn> rlock(*this);
   return _size;
}

/*
 * Get expand factor.
 */
template <typename T, typename Syn>
double array_list<T,Syn>::expand_factor() const {
   auto_read_lock<const Syn> rlock(*this);
   return _expand_factor;
}

/*
 * Set expand factor.
 * The expand factor must be > 0.
 * Return the new expand factor.
 */
template <typename T, typename Syn>
double array_list<T,Syn>::expand_factor(double e) {
   /* check argument */
   if (e <= double(0))
      throw ex_invalid_argument(
         "array_list expand factor must be > 0"
      );
   /* set expand factor */
   auto_write_lock<const Syn> wlock(*this);
   _expand_factor = e;
   _empty_factor  = e + double(1) - double(1)/(double(1) + e);
   return e;
}

/*
 * Return iterator over elements.
 */
template <typename T, typename Syn>
auto_ptr< iterator<T> > array_list<T,Syn>::iter_create() const {
   return auto_ptr< iterator<T> >(
      new array_list_iterator<T,Syn>(*this)
   );
}

/*
 * Return reverse iterator over elements.
 */
template <typename T, typename Syn>
auto_ptr< iterator<T> > array_list<T,Syn>::iter_reverse_create() const {
   return auto_ptr< iterator<T> >(
      new array_list_iterator_reverse<T,Syn>(*this)
   );
}

/*
 * Randomly permute the elements of the array list using the given random
 * number generator.
 *
 * Return an index array mapping resulting position --> original position.
 */
template <typename T, typename Syn>
lang::array<unsigned long> array_list<T,Syn>::randperm(rand_gen<>& r) {
   auto_write_lock<const Syn> wlock(*this);
   if (_size > 0) {
      return _data.randperm_subarray(
         _head_space,
         _head_space + _size - 1,
         r
      );
   } else {
      return lang::array<unsigned long>(0);
   }
}

/*
 * Randomly permute the elements of the specified subarray.
 * Return an index array mapping resulting position --> original position
 * (within the subarray, so index 0 corresponds to start of subarray).
 */
template <typename T, typename Syn>
lang::array<unsigned long> array_list<T,Syn>::randperm_subarray(
   unsigned long start,
   unsigned long end,
   rand_gen<>& r)
{
   auto_write_lock<const Syn> wlock(*this);
   #ifdef COLLECTIONS__ARRAY_LIST__CHECK_BOUNDS
      /* check array bounds */
      if (start >= _size)
         throw ex_index_out_of_bounds(
            "array list start index out of bounds", 
            start
         );
      else if (end >= _size)
         throw ex_index_out_of_bounds(
            "array list end index out of bounds", 
            end
         );
   #endif
   /* randperm */
   return _data.randperm_subarray(
      _head_space + start, 
      _head_space + end,
      r
   );
}

/*
 * Sort elements in ascending order according to the given comparison functor.
 *
 * Sorting uses an O(n + (n/p)*log(n/p)) expected time algorithm, where
 * n is the array size and p is the number of available processors.
 */
template <typename T, typename Syn>
void array_list<T,Syn>::sort(
   const comparable_functor<T>& f)
{
   const ptr_node_compare_functor f_compare(f);
   auto_write_lock<const Syn> wlock(*this);
   if (_size > 1) {
      _data.sort_subarray(
         _head_space, 
         _head_space + _size - 1, 
         f_compare
      );
   }
}

/*
 * Sort the array list.
 * Return an index array mapping sorted position --> original position.
 */
template <typename T, typename Syn>
lang::array<unsigned long> array_list<T,Syn>::sort_idx(
   const comparable_functor<T>& f)
{
   const ptr_node_compare_functor f_compare(f);
   auto_write_lock<const Syn> wlock(*this);
   if (_size > 0) {
      return _data.sort_subarray_idx(
         _head_space, 
         _head_space + _size - 1,
         f_compare
      );
   } else {
      return lang::array<unsigned long>(0);
   }
}

/*
 * Sort the specified subarray within the array list.
 */
template <typename T, typename Syn>
void array_list<T,Syn>::sort_subarray(
   unsigned long start,
   unsigned long end,
   const comparable_functor<T>& f)
{
   const ptr_node_compare_functor f_compare(f);
   auto_write_lock<const Syn> wlock(*this);
   #ifdef COLLECTIONS__ARRAY_LIST__CHECK_BOUNDS
      /* check array bounds */
      if (start >= _size)
         throw ex_index_out_of_bounds(
            "array list start index out of bounds", 
            start
         );
      else if (end >= _size)
         throw ex_index_out_of_bounds(
            "array list end index out of bounds", 
            end
         );
   #endif
   /* sort */
   _data.sort_subarray(
      _head_space + start, 
      _head_space + end,
      f_compare
   );
}

/*
 * Sort the specified subarray within the array list.
 * Return an index array mapping sorted position --> original position
 * (within the subarray, so index 0 corresponds to start of subarray).
 */
template <typename T, typename Syn>
lang::array<unsigned long> array_list<T,Syn>::sort_subarray_idx(
   unsigned long start,
   unsigned long end,
   const comparable_functor<T>& f)
{
   const ptr_node_compare_functor f_compare(f);
   auto_write_lock<const Syn> wlock(*this);
   #ifdef COLLECTIONS__ARRAY_LIST__CHECK_BOUNDS
      /* check array bounds */
      if (start >= _size)
         throw ex_index_out_of_bounds(
            "array list start index out of bounds", 
            start
         );
      else if (end >= _size)
         throw ex_index_out_of_bounds(
            "array list end index out of bounds", 
            end
         );
   #endif
   /* sort */
   return _data.sort_subarray_idx(
      _head_space + start, 
      _head_space + end,
      f_compare
   );
}

/*
 * Remove duplicate elements (as defined by the given comparison functor)
 * and place the remaining unique elements in sorted order.
 */
template <typename T, typename Syn>
void array_list<T,Syn>::unique(
   const comparable_functor<T>& f)
{
   const ptr_node_compare_functor f_compare(f);
   auto_write_lock<const Syn> wlock(*this);
   if (_size > 1) {
      /* sort */
      _data.sort_subarray(
         _head_space, 
         _head_space + _size - 1, 
         f_compare
      );
      /* move unique elements */
      unsigned long n_unique = 0;
      for (unsigned long n = 1; n < _size; n++) {
         if (f_compare(_data[_head_space+n], _data[_head_space+n_unique]) != 0)
            _data[_head_space + (++n_unique)] = _data[_head_space + n];
      }
      /* resize array */
      this->resize_shrink_tail(_size - 1 - n_unique);
   }
}

/*
 * Remove duplicate elements (as defined by the given comparison functor)
 * and place the remaining unique elements in sorted order.  In addition,
 * return an index array containing the original positions of the unique
 * elements.
 */
template <typename T, typename Syn>
lang::array<unsigned long> array_list<T,Syn>::unique_idx(
   const comparable_functor<T>& f)
{
   const ptr_node_compare_functor f_compare(f);
   auto_write_lock<const Syn> wlock(*this);
   if (_size > 1) {
      /* sort */
      lang::array<unsigned long> idx = _data.sort_subarray_idx(
         _head_space, 
         _head_space + _size - 1,
         f_compare
      );
      /* move unique elements */
      unsigned long n_unique = 0;
      for (unsigned long n = 1; n < _size; n++) {
         if (f_compare(_data[_head_space+n], _data[_head_space+n_unique]) != 0)
         {
            ++n_unique;
            _data[_head_space + n_unique] = _data[_head_space + n];
            idx[n_unique]                 = idx[n];
         }
      }
      /* resize array */
      this->resize_shrink_tail(_size - 1 - n_unique);
      idx.resize(n_unique + 1);
      return idx;
   } else {
      return lang::array<unsigned long>(_size);
   }
}

/***************************************************************************
 * Array list helper function implementation.
 ***************************************************************************/

/*
 * Resize the array list to have the given amount of head and tail space.
 */
template <typename T, typename Syn>
void array_list<T,Syn>::resize(
   unsigned long head_space,
   unsigned long tail_space)
{
   /* copy current array elements */
   lang::array< auto_ptr<array_list_node> > data(_size);
   for (unsigned long n = 0; n < _size; n++)
      data[n] = _data[_head_space + n];
   /* resize current array */
   _data.resize(head_space + _size + tail_space);
   /* place elements back into resized array */
   for (unsigned long n = 0; n < _size; n++)
      _data[head_space + n] = data[n];
   /* update head and tail space */
   _head_space = head_space;
   _tail_space = tail_space;
}

/*
 * Expand the head of the array list by the given number of elements.
 * Allocate more space and copy array as needed.
 */
template <typename T, typename Syn>
void array_list<T,Syn>::resize_expand_head(unsigned long n) {
   if (n > _head_space) {
      /* insufficient space - calculate expansion sizes */
      unsigned long new_size              = _size + n;
      unsigned long new_head_space_excess = static_cast<unsigned long>(
         math::ceil((static_cast<double>(new_size))*(_expand_factor))
      );
      _empty_space_max = static_cast<unsigned long>(
         math::ceil((static_cast<double>(new_size))*(_empty_factor))
      );
      /* resize array */
      this->resize(n + new_head_space_excess, _tail_space);
   }
   /* update array */
   _size       += n;
   _head_space -= n;
}

/*
 * Expand the tail of the array list by the given number of elements.
 * Allocate more space and copy array as needed.
 */
template <typename T, typename Syn>
void array_list<T,Syn>::resize_expand_tail(unsigned long n) {
  if (n > _tail_space) {
      /* insufficient space - calculate expansion sizes */
      unsigned long new_size              = _size + n;
      unsigned long new_tail_space_excess = static_cast<unsigned long>(
         math::ceil((static_cast<double>(new_size))*(_expand_factor))
      );
      _empty_space_max = static_cast<unsigned long>(
         math::ceil((static_cast<double>(new_size))*(_empty_factor))
      );
      /* resize array */
      this->resize(_head_space, n + new_tail_space_excess);
   }
   /* update array */
   _size       += n;
   _tail_space -= n;
}

/*
 * Shrink the head of the array list by the given number of elements.
 * Reallocate less space and copy array as needed.
 */
template <typename T, typename Syn>
void array_list<T,Syn>::resize_shrink_head(unsigned long n) {
   /* delete array list nodes */
   for (unsigned long i = 0; i < n; i++)
      _data[_head_space + i].reset();
   /* update size */
   _size       -= n;
   _head_space += n;
   /* shrink array if needed */
   if (_head_space > _empty_space_max) {
      /* compute new head space, empty space limit */
      unsigned long new_head_space = static_cast<unsigned long>(
         math::ceil((static_cast<double>(_size))*(_expand_factor))
      );
      _empty_space_max = static_cast<unsigned long>(
         math::ceil((static_cast<double>(_size))*(_empty_factor))
      );
      /* resize array */
      this->resize(new_head_space, _tail_space);
   }
}

/*
 * Shrink the tail of the array list by the given number of elements.
 * Reallocate less space and copy array as needed.
 */
template <typename T, typename Syn>
void array_list<T,Syn>::resize_shrink_tail(unsigned long n) {
   /* delete array list nodes */
   for (unsigned long i = (_size - n); i < _size; i++)
      _data[_head_space + i].reset();
   /* update size */
   _size       -= n;
   _tail_space += n;
   /* shrink array if needed */
   if (_tail_space > _empty_space_max) {
      /* compute new tail space, empty space limit */
      unsigned long new_tail_space = static_cast<unsigned long>(
         math::ceil((static_cast<double>(_size))*(_expand_factor))
      );
      _empty_space_max = static_cast<unsigned long>(
         math::ceil((static_cast<double>(_size))*(_empty_factor))
      );
      /* resize array */
      this->resize(_head_space, new_tail_space);
   }
}

} /* namespace collections */

#endif
