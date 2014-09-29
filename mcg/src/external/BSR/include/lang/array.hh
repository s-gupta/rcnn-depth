/*
 * Arrays (thread-safe).
 */
#ifndef LANG__ARRAY_HH
#define LANG__ARRAY_HH

#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_read_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/synchronizable.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "concurrent/threads/child_thread.hh"
#include "concurrent/threads/runnable.hh"
#include "config/safety.hh"
#include "functors/comparable_functors.hh"
#include "functors/filterable_functors.hh"
#include "functors/foldable_functors.hh"
#include "functors/iterable_functors.hh"
#include "functors/mappable_functors.hh"
#include "io/serialization/serial_input_stream.hh"
#include "io/serialization/serial_output_stream.hh"
#include "io/serialization/serializers.hh"
#include "lang/exceptions/ex_index_out_of_bounds.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "lang/exceptions/ex_not_found.hh"
#include "lang/iterators/iterator.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/random/generators/rand_gen.hh"
#include "math/random/generators/rand_gen_uniform.hh"

/* 
 * Enable/disable bounds checking for arrays.
 */
#ifdef CONFIG__SAFETY__CHECK_BOUNDS
   #define LANG__ARRAY__CHECK_BOUNDS (true)
#else
   #define LANG__ARRAY__CHECK_BOUNDS (false)
#endif

namespace lang {
/*
 * Imports.
 */
using concurrent::threads::synchronization::locks::auto_read_lock;
using concurrent::threads::synchronization::locks::auto_read_read_lock;
using concurrent::threads::synchronization::locks::auto_write_lock;
using concurrent::threads::synchronization::synchronizables::synchronizable;
using concurrent::threads::synchronization::synchronizables::unsynchronized;
using concurrent::threads::child_thread;
using concurrent::threads::runnable;
using functors::comparable_functor;
using functors::compare_functors;
using functors::filterable_functor;
using functors::foldable_functor;
using functors::iterable_functor;
using functors::mappable_functor;
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
using math::random::generators::rand_gen_uniform;

/*
 * Declare class for iterators over arrays.
 */
template <typename T, typename Syn>
class array_iterator;

template <typename T, typename Syn>
class array_iterator_reverse;

/*
 * Array class.
 */
template <typename T, typename Syn = unsynchronized>
class array : protected Syn {
public:
   /*
    * Friend classes.
    */
   friend class array_iterator<T,Syn>;
   friend class array_iterator<const T,Syn>;
   friend class array_iterator_reverse<T,Syn>;
   friend class array_iterator_reverse<const T,Syn>;
   template <typename U, typename S> friend class array;

   /*
    * Define the iterator types.
    */
   typedef array_iterator<T,Syn> iterator_t;
   typedef array_iterator<const T,Syn> const_iterator_t;
   typedef array_iterator_reverse<T,Syn> iterator_reverse_t;
   typedef array_iterator_reverse<const T,Syn> const_iterator_reverse_t;

   /*
    * Constructors.
    */
   array();
   explicit array(unsigned long /* size */);
   explicit array(unsigned long /* size */, const T& /* init value */);

   /*
    * Copy constructors.
    */
   array(const array<T,Syn>&);
   template <typename U, typename S> array(const array<U,S>&);

   /*
    * Destructor.
    */
   virtual ~array();

   /*
    * Assignment operator.
    */
   array<T,Syn>& operator=(const array<T,Syn>&);
   template <typename U, typename S> array<T,Syn>& operator=(const array<U,S>&);

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
   static auto_ptr< array<T,Syn> > deserialize(
      serial_input_stream&,
      const serializer<T>& = serializers<T>::s_default()
   );

   /*
    * Element reference.
    */
   inline T& operator[](unsigned long);
   inline T& operator()(unsigned long);
   inline const T& operator[](unsigned long) const;
   inline const T& operator()(unsigned long) const;

   /*
    * Subarray operations.
    * Return array of elements in the specified range.
    */
   array<T,Syn> subarray(
      unsigned long,                /* start index */ 
      unsigned long                 /* end index   */
   ) const;
   
   array<T,Syn> subarray(
      unsigned long,                /* start index */ 
      unsigned long,                /* step size   */
      unsigned long                 /* end index   */
   ) const;
   
   array<T,Syn> subarray(
      const array<unsigned long>&   /* indices of desired elements */
   ) const;

   /*
    * Reverse element order.
    * Return a reference to the array.
    */
   virtual array<T,Syn>& reverse();

   /*
    * Size.
    */
   inline bool is_empty() const;
   inline unsigned long size() const;

   /*
    * Resize.
    * Default initialize new elements (if enlarging).
    * Return a reference to the array.
    */
   virtual array<T,Syn>& resize(unsigned long);

   /*
    * Return iterator over elements.
    */
   auto_ptr< iterator<T> >       iter_create();
   auto_ptr< iterator<const T> > iter_create() const;
   auto_ptr< iterator<T> >       iter_reverse_create();
   auto_ptr< iterator<const T> > iter_reverse_create() const;

   /*
    * Iterate a functor over the array.
    */
   void iter(const iterable_functor<T>&);
   void iter(const iterable_functor<const T>&) const;

   /*
    * Apply a filter to the array.
    */
   array<T,Syn> filter(const filterable_functor<T>&);
   array<T,Syn> filter(const filterable_functor<const T>&) const;
   
   /*
    * Check if some/all of the elements pass a filter.
    */
   bool exists(const filterable_functor<T>&);
   bool exists(const filterable_functor<const T>&) const;
   bool for_all(const filterable_functor<T>&);
   bool for_all(const filterable_functor<const T>&) const;
   
   /*
    * Apply a fold functor to the array.
    */
   template <typename U> U& fold(const foldable_functor<T,U>&, U&);
   template <typename U> U& fold(const foldable_functor<const T,U>&, U&) const;

   /*
    * Apply a map to elements of the array.
    */
   template <typename U> array<U,Syn> map(const mappable_functor<T,U>&);
   template <typename U> array<U,Syn> map(const mappable_functor<const T,U>&) const;

   /*
    * Randomly permute the elements of the array.
    * Return an index array mapping resulting position --> original position.
    *
    * The generator from which to draw values used in computing the random
    * permutation may optionally be specified.  If unspecified, a generator 
    * is created and used for this purpose.
    */
   virtual array<unsigned long> randperm();
 
   virtual array<unsigned long> randperm(rand_gen<>&);

   /*
    * Randomly permute the elements of the specified subarray within the array.
    * Return an index array mapping resulting position --> original position
    * (within the subarray, so index 0 corresponds to start of subarray).
    */
   virtual array<unsigned long> randperm_subarray(
      unsigned long,    /* start index */ 
      unsigned long     /* end index   */
   );

   virtual array<unsigned long> randperm_subarray(
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
   virtual void sort(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Sort the array.
    * Return an index array mapping sorted position --> original position.
    */
   virtual array<unsigned long> sort_idx(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );
  
   /*
    * Sort the specified subarray within the array.
    */
   virtual void sort_subarray(
      unsigned long,    /* start index */ 
      unsigned long,    /* end index   */
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Sort the specified subarray within the array.
    * Return an index array mapping sorted position --> original position
    * (within the subarray, so index 0 corresponds to start of subarray).
    */
   virtual array<unsigned long> sort_subarray_idx(
      unsigned long,    /* start index */ 
      unsigned long,    /* end index   */
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Remove duplicate elements (as defined by the given comparison functor)
    * and place the remaining unique elements in sorted order.
    */
   virtual void unique(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Remove duplicate elements (as defined by the given comparison functor)
    * and place the remaining unique elements in sorted order.  In addition,
    * return an index array containing the original positions of the unique
    * elements.
    */
   virtual array<unsigned long> unique_idx(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

protected:
   /************************************************************************
    * Array data structure.
    ************************************************************************/
    
   unsigned long _size; /* size of array  */
   T* _data;            /* array elements */
   
   /************************************************************************
    * Array helper functions.
    ************************************************************************/
   
   /*
    * Swap the contents of two arrays.
    * Note: Arrays are not locked by this function.
    */
   static void swap(array<T,Syn>&, array<T,Syn>&);

   /*
    * Implementation of duplicate element removal.
    * The terminator space argument is used by strings (for '\0' terminator).
    * Note: Arrays are not locked by these methods.
    */
   void make_unique(
      const comparable_functor<T>&, unsigned long = 0 /* terminator space */
   );
   
   array<unsigned long> make_unique_idx(
      const comparable_functor<T>&, unsigned long = 0 /* terminator space */
   );
   
   /*
    * Runnable object for quick sorting a specified subrange of an array in
    * parallel, given the start and stop indicates, comparison functor, and
    * the number of additional processors that can be used during the sort.
    *
    * This is an O(n + (n/p)*log(n/p)) expected time algorithm, where
    * n is the array size and p is the number of available processors.
    */
   class sorter : public runnable {
   public:
      /*
       * Constructor.
       */
      explicit sorter(
         T*,                                    /* array to sort */
         unsigned long,                         /* start of subrange to sort */
         unsigned long,                         /* end of subrange to sort */
         const comparable_functor<T>&           /* comparison functor */
      );

      /*
       * Destructor.
       */
      virtual ~sorter();
      
      /*
       * Run the sort.
       */
      virtual void run();
      
   protected:
      T* _data;                                 /* array to sort */
      unsigned long _low;                       /* start of subrange to sort */
      unsigned long _high;                      /* end of subrange to sort */
      const comparable_functor<T>& _f_compare;  /* comparison functor */
   };

   /*
    * Runnable object for quick sorting and also updating an index array
    * mapping sorted position --> original position.
    */
   class sorter_idx : public sorter {
   public:
      /*
       * Constructor.
       */
      explicit sorter_idx(
         T*,                                    /* array to sort */
         unsigned long*,                        /* index array */
         unsigned long,                         /* start of subrange to sort */
         unsigned long,                         /* end of subrange to sort */
         const comparable_functor<T>&           /* comparison functor */
      );

      /*
       * Destructor.
       */
      virtual ~sorter_idx();
      
      /*
       * Run the sort.
       */
      virtual void run();
      
   protected:
      unsigned long* _idx;                      /* index array */
   };
};

/*
 * Array iterator.
 */
template <typename T, typename Syn = unsynchronized>
class array_iterator : public iterator<T> {
public:
   /*
    * Constructor.
    * Seek to the head of the array.
    */
   array_iterator(array<T,Syn>&);

   /*
    * Copy constructor.
    */
   array_iterator(const array_iterator<T,Syn>&);

   /*
    * Destructor.
    */
   virtual ~array_iterator();

   /*
    * Check if there is another item available.
    */
   bool has_next() const;

   /*
    * Return the next item.
    * Throw an exception (not found) if there are no more items.
    */
   T& next();
   
protected:
   array<T,Syn>& _arr;  /* array being iterated over */
   unsigned long _curr; /* current position in array */
};

/*
 * Array iterator for const arrays.
 */
template <typename T, typename Syn>
class array_iterator<const T,Syn> : public iterator<const T> {
public:
   /*
    * Constructor.
    * Seek to the head of the array.
    */
   array_iterator(const array<T,Syn>&);

   /*
    * Copy constructor.
    */
   array_iterator(const array_iterator<const T,Syn>&);

   /*
    * Destructor.
    */
   virtual ~array_iterator();

   /*
    * Check if there is another item available.
    */
   bool has_next() const;

   /*
    * Return the next item.
    * Throw an exception (not found) if there are no more items.
    */
   const T& next();
   
protected:
   const array<T,Syn>& _arr;  /* array being iterated over */
   unsigned long      _curr;  /* current position in array */
};

/*
 * Array reverse iterator.
 */
template <typename T, typename Syn = unsynchronized>
class array_iterator_reverse : public iterator<T> {
public:
   /*
    * Constructor.
    * Seek to the head of the array.
    */
   array_iterator_reverse(array<T,Syn>&);

   /*
    * Copy constructor.
    */
   array_iterator_reverse(const array_iterator_reverse<T,Syn>&);

   /*
    * Destructor.
    */
   virtual ~array_iterator_reverse();

   /*
    * Check if there is another item available.
    */
   bool has_next() const;

   /*
    * Return the next item.
    * Throw an exception (not found) if there are no more items.
    */
   T& next();
   
protected:
   array<T,Syn>& _arr;  /* array being iterated over */
   unsigned long _curr; /* current position in array */
};

/*
 * Array reverse iterator for const arrays.
 */
template <typename T, typename Syn>
class array_iterator_reverse<const T,Syn> : public iterator<const T> {
public:
   /*
    * Constructor.
    * Seek to the head of the array.
    */
   array_iterator_reverse(const array<T,Syn>&);

   /*
    * Copy constructor.
    */
   array_iterator_reverse(const array_iterator_reverse<const T,Syn>&);

   /*
    * Destructor.
    */
   virtual ~array_iterator_reverse();

   /*
    * Check if there is another item available.
    */
   bool has_next() const;

   /*
    * Return the next item.
    * Throw an exception (not found) if there are no more items.
    */
   const T& next();
   
protected:
   const array<T,Syn>& _arr;  /* array being iterated over */
   unsigned long      _curr;  /* current position in array */
};

/***************************************************************************
 * Array iterator implementation.
 ***************************************************************************/

/*
 * Constructors.
 * Lock array and initialize position.
 */
template <typename T, typename Syn>
array_iterator<T,Syn>::array_iterator(array<T,Syn>& a)
 : _arr(a),
   _curr(0)
{
   _arr.read_lock();
}

template <typename T, typename Syn>
array_iterator<const T,Syn>::array_iterator(const array<T,Syn>& a)
 : _arr(a),
   _curr(0)
{
   _arr.read_lock();
}

template <typename T, typename Syn>
array_iterator_reverse<T,Syn>::array_iterator_reverse(
   array<T,Syn>& a)
 : _arr(a) 
{
   _arr.read_lock();
   _curr = _arr._size;
}

template <typename T, typename Syn>
array_iterator_reverse<const T,Syn>::array_iterator_reverse(
   const array<T,Syn>& a)
 : _arr(a)
{
   _arr.read_lock();
   _curr = _arr._size;
}

/*
 * Copy constructors.
 */
template <typename T, typename Syn>
array_iterator<T,Syn>::array_iterator(
   const array_iterator<T,Syn>& i)
 : _arr(i._arr),
   _curr(i._curr)
{
   _arr.read_lock();
}

template <typename T, typename Syn>
array_iterator<const T,Syn>::array_iterator(
   const array_iterator<const T,Syn>& i)
 : _arr(i._arr),
   _curr(i._curr)
{
   _arr.read_lock();
}

template <typename T, typename Syn>
array_iterator_reverse<T,Syn>::array_iterator_reverse(
   const array_iterator_reverse<T,Syn>& i)
 : _arr(i._arr),
   _curr(i._curr)
{
   _arr.read_lock();
}

template <typename T, typename Syn>
array_iterator_reverse<const T,Syn>::array_iterator_reverse(
   const array_iterator_reverse<const T,Syn>& i)
 : _arr(i._arr),
   _curr(i._curr)
{
   _arr.read_lock();
}

/*
 * Destructors.
 * Unlock the array.
 */
template <typename T, typename Syn>
array_iterator<T,Syn>::~array_iterator() {
   _arr.read_unlock();
}

template <typename T, typename Syn>
array_iterator<const T,Syn>::~array_iterator() {
   _arr.read_unlock();
}

template <typename T, typename Syn>
array_iterator_reverse<T,Syn>::~array_iterator_reverse() {
   _arr.read_unlock();
}

template <typename T, typename Syn>
array_iterator_reverse<const T,Syn>::~array_iterator_reverse() {
   _arr.read_unlock();
}

/*
 * Check if there is another item available.
 */
template <typename T, typename Syn>
bool array_iterator<T,Syn>::has_next() const {
   return (_curr < _arr._size);
}

template <typename T, typename Syn>
bool array_iterator<const T,Syn>::has_next() const {
   return (_curr < _arr._size);
}

template <typename T, typename Syn>
bool array_iterator_reverse<T,Syn>::has_next() const {
   return (_curr > 0);
}

template <typename T, typename Syn>
bool array_iterator_reverse<const T,Syn>::has_next() const {
   return (_curr > 0);
}

/*
 * Return the next item.
 * Throw an exception (not found) if there are no more items.
 */
template <typename T, typename Syn>
T& array_iterator<T,Syn>::next() {
   if (_curr < _arr._size) {
      T& t = _arr._data[_curr];
      _curr++;
      return t;
   } else {
      throw ex_not_found(
         "iterator attempted to read past end of array"
      );
   }
}

template <typename T, typename Syn>
const T& array_iterator<const T,Syn>::next() {
   if (_curr < _arr._size) {
      const T& t = _arr._data[_curr];
      _curr++;
      return t;
   } else {
      throw ex_not_found(
         "iterator attempted to read past end of array"
      );
   }
}

template <typename T, typename Syn>
T& array_iterator_reverse<T,Syn>::next() {
   if (_curr > 0) {
      _curr--;
      return _arr._data[_curr];
   } else {
      throw ex_not_found(
         "reverse iterator attempted to read past start of array"
      );
   }
}

template <typename T, typename Syn>
const T& array_iterator_reverse<const T,Syn>::next() {
   if (_curr > 0) {
      _curr--;
      return _arr._data[_curr];
   } else {
      throw ex_not_found(
         "reverse iterator attempted to read past start of array"
      );
   }
}

/***************************************************************************
 * Array initialization.
 ***************************************************************************/

/*
 * Default array initialization function for user-defined types.
 * Do nothing as the constructor of the user-defined type will be called to
 * initialization the array contents upon array creation.
 */
template <typename T>
void array_initialize(unsigned long size, T* data) {
   /* do nothing */
}

/*
 * Array initialization function for pointer types.
 * This function guarantees that arrays of pointers are initialized to NULL.
 */
template <typename T>
void array_initialize(unsigned long size, T** data) {
   for (unsigned long n = 0; n < size; n++)
      data[n] = NULL;
}
 
/*
 * Array initialization functions for built-in types.
 * These functions guarantee that arrays are always zero-initialized.
 */
void array_initialize(unsigned long, bool*);
void array_initialize(unsigned long, char*);
void array_initialize(unsigned long, unsigned char*);
void array_initialize(unsigned long, short*);
void array_initialize(unsigned long, unsigned short*);
void array_initialize(unsigned long, int*);
void array_initialize(unsigned long, unsigned int*);
void array_initialize(unsigned long, long*);
void array_initialize(unsigned long, unsigned long*);
void array_initialize(unsigned long, long long*);
void array_initialize(unsigned long, unsigned long long*);
void array_initialize(unsigned long, float*);
void array_initialize(unsigned long, double*);
void array_initialize(unsigned long, long double*);

/***************************************************************************
 * Array implementation.
 ***************************************************************************/
 
/*
 * Swap the contents of two arrays.
 * Note: Arrays are not locked by this function.
 */
template <typename T, typename Syn>
void array<T,Syn>::swap(array<T,Syn>& a0, array<T,Syn>& a1) {
   unsigned long temp_size = a0._size;
   T*            temp_data = a0._data;
   a0._size = a1._size;
   a0._data = a1._data;
   a1._size = temp_size;
   a1._data = temp_data;
}

/*
 * Default constructor.
 * Create empty array.
 */
template <typename T, typename Syn>
array<T,Syn>::array()
 : Syn(),
   _size(0),
   _data(NULL)
{ }

/*
 * Constructor.
 * Create array of specified size.
 */
template <typename T, typename Syn>
array<T,Syn>::array(unsigned long size)
 : Syn(),
   _size(size),
   _data((size > 0) ? (new T[size]()) : NULL)
{ 
   array_initialize(_size, _data);
}

/*
 * Constructor.
 * Create array of specified size and initialize elements to the given value.
 */
template <typename T, typename Syn>
array<T,Syn>::array(unsigned long size, const T& t)
 : Syn(),
   _size(0),
   _data(NULL)
{
   /* create and initialize array */
   array<T,Syn> a(size);
   const T t_copy(t);
   for (unsigned long n = 0; n < size; n++)
      a._data[n] = t_copy;
   /* take ownership of array */
   array<T,Syn>::swap(*this, a);
}

/*
 * Copy constructor.
 */
template <typename T, typename Syn>
array<T,Syn>::array(const array<T,Syn>& a)
 : Syn(),
   _size(0),
   _data(NULL)
{
   auto_read_lock<const Syn> rlock(a);
   /* copy source array */
   array<T,Syn> a_copy(a._size);
   const T* a_data = a._data;
   for (unsigned long n = 0; n < a._size; n++)
      a_copy._data[n] = a_data[n];
   /* take ownership of copy */
   array<T,Syn>::swap(*this, a_copy);
}

/*
 * Copy constructor (type conversion).
 */
template <typename T, typename Syn>
template <typename U, typename S>
array<T,Syn>::array(const array<U,S>& a)
 : Syn(),
   _size(0),
   _data(NULL)
{
   auto_read_lock<const S> rlock(a);
   /* copy source array */
   array<T,Syn> a_copy(a._size);
   const U* a_data = a._data;
   for (unsigned long n = 0; n < a._size; n++)
      a_copy._data[n] = a_data[n];
   /* take ownership of copy */
   array<T,Syn>::swap(*this, a_copy);
}

/*
 * Destructor.
 * Delete the underlying raw array.
 */
template <typename T, typename Syn>
array<T,Syn>::~array() {
   delete [] _data;
}

/*
 * Assignment operator.
 */
template <typename T, typename Syn>
array<T,Syn>& array<T,Syn>::operator=(const array<T,Syn>& a) {
   /* copy source array */
   array<T,Syn> a_copy(a);
   /* swap copy into current array */
   auto_write_lock<const Syn> wlock(*this);
   array<T,Syn>::swap(*this, a_copy);
   return *this;
}

/*
 * Assignment operator (type conversion).
 */
template <typename T, typename Syn>
template <typename U, typename S>
array<T,Syn>& array<T,Syn>::operator=(const array<U,S>& a) {
   /* copy source array */
   array<T,Syn> a_copy(a);
   /* swap copy into current array */
   auto_write_lock<const Syn> wlock(*this);
   array<T,Syn>::swap(*this, a_copy);
   return *this;
}

/*
 * Serialize.
 */
template <typename T, typename Syn>
void array<T,Syn>::serialize(
   serial_output_stream& s, const serializer<T>& slzr) const
{
   auto_read_lock<const Syn> rlock(*this);
   s << _size;
   for (unsigned long n = 0; n < _size; n++)
      slzr.serialize(s, _data[n]);
}

/*
 * Deserialize.
 */
template <typename T, typename Syn>
auto_ptr< array<T,Syn> > array<T,Syn>::deserialize(
   serial_input_stream& s, const serializer<T>& slzr)
{
   unsigned long size = 0;
   s >> size;
   auto_ptr< array<T,Syn> > a(new array<T,Syn>(size));
   for (unsigned long n = 0; n < size; n++) {
      auto_ptr<T> t = slzr.deserialize(s);
      a->_data[n] = *t;
   }
   return a;
}

/*
 * Element reference on non-const array.
 * Return reference to non-const element.
 */
template <typename T, typename Syn>
inline T& array<T,Syn>::operator[](unsigned long n) {
   auto_read_lock<const Syn> rlock(*this);
   #ifdef LANG__ARRAY__CHECK_BOUNDS
      /* perform bounds check */
      if (n >= _size)
         throw ex_index_out_of_bounds(
            "array index out of bounds", 
            n
         );
   #endif
   /* retrieve element */
   return _data[n];
}

/*
 * Functor form of element reference on non-const array.
 * Return reference to non-const element.
 */
template <typename T, typename Syn>
inline T& array<T,Syn>::operator()(unsigned long n) {
   return this->operator[](n);
}

/*
 * Element reference on const array.
 * Return reference to const element.
 */
template <typename T, typename Syn>
inline const T& array<T,Syn>::operator[](unsigned long n) const {
   auto_read_lock<const Syn> rlock(*this);
   #ifdef LANG__ARRAY__CHECK_BOUNDS
      /* perform bounds check */
      if (n >= _size)
         throw ex_index_out_of_bounds(
            "array index out of bounds", 
            n
         );
   #endif
   /* retrieve element */
   return _data[n];
}

/*
 * Functor form of element reference on const array.
 * Return reference to const element.
 */
template <typename T, typename Syn>
inline const T& array<T,Syn>::operator()(unsigned long n) const {
   return this->operator[](n);
}

/*
 * Return subarray of elements in range [start:end].
 */
template <typename T, typename Syn>
array<T,Syn> array<T,Syn>::subarray(
   unsigned long start,
   unsigned long end) const
{
   return this->subarray(start, 1, end);
}

/*
 * Return subarray of elements in range [start:step:end].
 */
template <typename T, typename Syn>
array<T,Syn> array<T,Syn>::subarray(
   unsigned long start,
   unsigned long step,
   unsigned long end) const
{
   auto_read_lock<const Syn> rlock(*this);
   #ifdef LANG__ARRAY__CHECK_BOUNDS
      /* check array bounds */
      if (start >= _size)
         throw ex_index_out_of_bounds(
            "array start index out of bounds", 
            start
         );
      else if (end >= _size)
         throw ex_index_out_of_bounds(
            "array end index out of bounds", 
            end
         );
   #endif
   /* check step size */
   if (step == 0)
      throw ex_invalid_argument(
         "array step must be nonzero"
      );
   /* grab subarray */
   array<T,Syn> arr;
   arr._size = (start <= end) ? ((end - start)/step + 1) : 0;
   arr._data = (arr._size > 0) ? (new T[arr._size]()) : NULL;
   const T* data = _data;
   for (unsigned long n = 0; n < arr._size; n++, start += step)
      arr._data[n] = data[start];
   return arr;
}

/*
 * Return the array with elements taken from the given indices.
 */
template <typename T, typename Syn>
array<T,Syn> array<T,Syn>::subarray(
   const array<unsigned long>& index_arr) const
{
   auto_read_read_lock<const synchronizable> rrlock(*this, index_arr);
   const unsigned long* indices = index_arr._data;
   unsigned long      n_indices = index_arr._size;
   array<T,Syn> arr(n_indices);
   const T* data = _data;
   for (unsigned long n = 0; n < n_indices; n++) {
      /* get index */
      unsigned long index = indices[n];
      #ifdef LANG__ARRAY__CHECK_BOUNDS
         /* check index bounds */
         if (index >= _size)
            throw ex_index_out_of_bounds(
               "array index out of bounds", 
               index
            );
      #endif
      /* copy the element */
      arr._data[n] = data[index];
   }
   return arr;
}

/*
 * Reverse order of elements in array.
 */
template <typename T, typename Syn>
array<T,Syn>& array<T,Syn>::reverse() {
   auto_write_lock<const Syn> wlock(*this);
   for (unsigned long start = 0, end = _size; start < end; start++) {
      end--;
      T temp       = _data[start];
      _data[start] = _data[end];
      _data[end]   = temp;
   }
   return *this;
}

/*
 * Check if array is empty (size zero).
 */
template <typename T, typename Syn>
inline bool array<T,Syn>::is_empty() const {
   auto_read_lock<const Syn> rlock(*this);
   return (_size == 0);
}

/*
 * Get number of elements in array.
 */
template <typename T, typename Syn>
inline unsigned long array<T,Syn>::size() const {
   auto_read_lock<const Syn> rlock(*this);
   return _size;
}

/*
 * Resize.
 * Default initialize new elements (if enlarging).
 * Return a reference to the array.
 */
template <typename T, typename Syn>
array<T,Syn>& array<T,Syn>::resize(unsigned long size) {
   /* allocate resized array */
   array<T,Syn> a(size);
   /* lock current array */
   auto_write_lock<const Syn> wlock(*this);
   /* copy current array into part of resized array */
   unsigned long size_min = (_size < a._size) ? _size : a._size;
   for (unsigned long n = 0; n < size_min; n++)
      a._data[n] = _data[n];
   /* swap resized array into current */
   array<T,Syn>::swap(*this, a);
   return *this;
}

/*
 * Return iterator over elements.
 */
template <typename T, typename Syn>
auto_ptr< iterator<T> > array<T,Syn>::iter_create() {
   return auto_ptr< iterator<T> >(
      new array_iterator<T,Syn>(*this)
   );
}

/*
 * Return iterator over const elements.
 */
template <typename T, typename Syn>
auto_ptr< iterator<const T> > array<T,Syn>::iter_create() const {
   return auto_ptr< iterator<const T> >(
      new array_iterator<const T,Syn>(*this)
   );
}

/*
 * Return reverse iterator over elements.
 */
template <typename T, typename Syn>
auto_ptr< iterator<T> > array<T,Syn>::iter_reverse_create() {
   return auto_ptr< iterator<T> >(
      new array_iterator_reverse<T,Syn>(*this)
   );
}

/*
 * Return reverse iterator over const elements.
 */
template <typename T, typename Syn>
auto_ptr< iterator<const T> > array<T,Syn>::iter_reverse_create() const {
   return auto_ptr< iterator<const T> >(
      new array_iterator_reverse<const T,Syn>(*this)
   );
}

/*
 * Iterate a functor over a non-const array.
 */
template <typename T, typename Syn>
void array<T,Syn>::iter(const iterable_functor<T>& f) {
   auto_write_lock<const Syn> wlock(*this);
   for (unsigned long n = 0; n < _size; n++)
      f(_data[n]);
}

/*
 * Iterate a functor over a const array.
 */
template <typename T, typename Syn>
void array<T,Syn>::iter(const iterable_functor<const T>& f) const {
   auto_read_lock<const Syn> rlock(*this);
   for (unsigned long n = 0; n < _size; n++)
      f(_data[n]);
}

/*
 * Apply filter to non-const array.
 */
template <typename T, typename Syn>
array<T,Syn> array<T,Syn>::filter(const filterable_functor<T>& f) {
   auto_write_lock<const Syn> wlock(*this);
   /* compute filter results */
   array<bool> filter_result(_size, false);
   unsigned long filtered_size = 0;
   for (unsigned long n = 0; n < _size; n++) {
      filter_result[n] = f(_data[n]);
      if (filter_result[n])
         filtered_size++;
   }
   /* copy elements which passed filter */
   array<T,Syn> a(filtered_size);
   unsigned long n_filt = 0;
   for (unsigned long n = 0; n < _size; n++) {
      if (filter_result[n]) {
         a._data[n_filt] = _data[n];
         n_filt++;
      }
   }
   return a;
}

/*
 * Apply filter to const array.
 */
template <typename T, typename Syn>
array<T,Syn> array<T,Syn>::filter(const filterable_functor<const T>& f) const {
   auto_read_lock<const Syn> rlock(*this);
   /* compute filter results */
   array<bool> filter_result(_size, false);
   unsigned long filtered_size = 0;
   for (unsigned long n = 0; n < _size; n++) {
      filter_result[n] = f(_data[n]);
      if (filter_result[n])
         filtered_size++;
   }
   /* copy elements which passed filter */
   array<T,Syn> a(filtered_size);
   unsigned long n_filt = 0;
   for (unsigned long n = 0; n < _size; n++) {
      if (filter_result[n]) {
         a._data[n_filt] = _data[n];
         n_filt++;
      }
   }
   return a;
}

/*
 * Test if functor returns true on some element in non-const array.
 * Note that testing is terminated as soon as such an element is found.
 */
template <typename T, typename Syn>
bool array<T,Syn>::exists(const filterable_functor<T>& f) {
   auto_write_lock<const Syn> wlock(*this);
   bool found = false;
   for (unsigned long n = 0; (n < _size) && (!found); n++)
      found = f(_data[n]);
   return found;
}

/*
 * Test if functor returns true on some element in const array.
 * Note that testing is terminated as soon as such an element is found.
 */
template <typename T, typename Syn>
bool array<T,Syn>::exists(const filterable_functor<const T>& f) const {
   auto_read_lock<const Syn> rlock(*this);
   bool found = false;
   for (unsigned long n = 0; (n < _size) && (!found); n++)
      found = f(_data[n]);
   return found;
}
 
/*
 * Test if functor returns true on all elements in non-const array.
 * Note that testing is terminated as soon as a false value is found.
 */
template <typename T, typename Syn>
bool array<T,Syn>::for_all(const filterable_functor<T>& f) {
   auto_write_lock<const Syn> wlock(*this);
   bool all = true;
   for (unsigned long n = 0; (n < _size) && (all); n++)
      all = f(_data[n]);
   return all;
}

/*
 * Test if functor returns true on all elements in const array.
 * Note that testing is terminated as soon as a false value is found.
 */
template <typename T, typename Syn>
bool array<T,Syn>::for_all(const filterable_functor<const T>& f) const {
   auto_read_lock<const Syn> rlock(*this);
   bool all = true;
   for (unsigned long n = 0; (n < _size) && (all); n++)
      all = f(_data[n]);
   return all;
}

/*
 * Fold functor over non-const array.
 */
template <typename T, typename Syn>
template <typename U>
U& array<T,Syn>::fold(const foldable_functor<T,U>& f, U& u) {
   auto_write_lock<const Syn> wlock(*this);
   U* u_result = &u;
   for (unsigned long n = 0; n < _size; n++) {
      U& u_temp = f(_data[n], *u_result);
      u_result = &u_temp;
   }
   return *u_result;
}

/*
 * Fold functor over const array.
 */
template <typename T, typename Syn>
template <typename U>
U& array<T,Syn>::fold(const foldable_functor<const T,U>& f, U& u) const {
   auto_read_lock<const Syn> rlock(*this);
   U* u_result = &u;
   for (unsigned long n = 0; n < _size; n++) {
      U& u_temp = f(_data[n], *u_result);
      u_result = &u_temp;
   }
   return *u_result;
}

/*
 * Map elements of a non-const array.
 */
template <typename T, typename Syn>
template <typename U>
array<U,Syn> array<T,Syn>::map(const mappable_functor<T,U>& f) {
   auto_write_lock<const Syn> wlock(*this);
   array<U,Syn> a(_size);
   for (unsigned long n = 0; n < _size; n++)
      a._data[n] = f(_data[n]);
   return a;
}

/*
 * Map elements of a const array.
 */
template <typename T, typename Syn>
template <typename U>
array<U,Syn> array<T,Syn>::map(const mappable_functor<const T,U>& f) const {
   auto_read_lock<const Syn> rlock(*this);
   array<U,Syn> a(_size);
   for (unsigned long n = 0; n < _size; n++)
      a._data[n] = f(_data[n]);
   return a;
}

/*
 * Randomly permute the elements of the array.
 * Return an index array mapping resulting position --> original position.
 */
template <typename T, typename Syn>
array<unsigned long> array<T,Syn>::randperm() {
   rand_gen_uniform<> r;
   return this->randperm(r);
}

/*
 * Randomly permute the elements of the array using the given random 
 * number generator.
 *
 * Return an index array mapping resulting position --> original position.
 */
template <typename T, typename Syn>
array<unsigned long> array<T,Syn>::randperm(rand_gen<>& r) {
   auto_write_lock<const Syn> wlock(*this);
   /* generate and sort random-valued array */
   array<double> rand_arr(_size);
   for (unsigned long n = 0; n < _size; n++)
      rand_arr._data[n] = r.generate();
   array<unsigned long> idx = rand_arr.sort_idx();
   /* create array ordered by resulting index */
   array<T,Syn> a(_size);
   for (unsigned long n = 0; n < _size; n++)
      a._data[n] = _data[idx._data[n]];
   /* take ownership of reordered array */
   array<T,Syn>::swap(*this, a);
   return idx;
}

/*
 * Randomly permute the elements of the specified subarray within the array.
 * Return an index array mapping resulting position --> original position
 * (within the subarray, so index 0 corresponds to start of subarray).
 */
template <typename T, typename Syn>
array<unsigned long> array<T,Syn>::randperm_subarray(
   unsigned long start,
   unsigned long end)
{
   rand_gen_uniform<> r;
   return this->randperm_subarray(start, end, r);
}

/*
 * Randomly permute the elements of the specified subarray using the given 
 * random number generator.
 */
template <typename T, typename Syn>
array<unsigned long> array<T,Syn>::randperm_subarray(
   unsigned long start,
   unsigned long end,
   rand_gen<>& r)
{
   auto_write_lock<const Syn> wlock(*this);
   #ifdef LANG__ARRAY__CHECK_BOUNDS
      /* check array bounds */
      if (start >= _size)
         throw ex_index_out_of_bounds(
            "array start index out of bounds", 
            start
         );
      else if (end >= _size)
         throw ex_index_out_of_bounds(
            "array end index out of bounds", 
            end
         );
   #endif
   /* generate and sort random-valued array */
   unsigned long idx_size = (start <= end) ? (end - start + 1) : 0;
   array<double> rand_arr(idx_size);
   for (unsigned long n = 0; n < idx_size; n++)
      rand_arr._data[n] = r.generate();
   array<unsigned long> idx = rand_arr.sort_idx();
   /* create array ordered by resulting index */
   array<T,Syn> a(_size);
   for (unsigned long n = 0; n < start; n++)
      a._data[n] = _data[n];
   for (unsigned long n = start; n <= end; n++)
      a._data[n] = _data[idx._data[n - start] + start];
   for (unsigned long n = end + 1; n < _size; n++)
      a._data[n] = _data[n];
   /* take ownership of reordered array */
   array<T,Syn>::swap(*this, a);
   return idx;
}

/*
 * Sort elements in ascending order according to the given comparison functor.
 *
 * Sorting is an O(n + (n/p)*log(n/p)) expected time algorithm, where
 * n is the array size and p is the number of available processors.
 */
template <typename T, typename Syn>
void array<T,Syn>::sort(const comparable_functor<T>& f) {
   auto_write_lock<const Syn> wlock(*this);
   if (_size > 1) {
      sorter s(_data, 0, (_size - 1), f);
      s.run();
   }
}

/*
 * Sort the array.
 * Return an index array mapping sorted position --> original position.
 */
template <typename T, typename Syn>
array<unsigned long> array<T,Syn>::sort_idx(const comparable_functor<T>& f) {
   auto_write_lock<const Syn> wlock(*this);
   /* initialize index array */
   array<unsigned long> idx(_size);
   for (unsigned long n = 0; n < _size; n++)
      idx._data[n] = n;
   /* sort */
   if (_size > 1) {
      sorter_idx s(_data, idx._data, 0, (_size - 1), f);
      s.run();
   }
   return idx;
}

/*
 * Sort the specified subarray within the array.
 */
template <typename T, typename Syn>
void array<T,Syn>::sort_subarray(
   unsigned long start, 
   unsigned long end, 
   const comparable_functor<T>& f)
{
   auto_write_lock<const Syn> wlock(*this);
   #ifdef LANG__ARRAY__CHECK_BOUNDS
      /* check array bounds */
      if (start >= _size)
         throw ex_index_out_of_bounds(
            "array start index out of bounds", 
            start
         );
      else if (end >= _size)
         throw ex_index_out_of_bounds(
            "array end index out of bounds", 
            end
         );
   #endif
   /* sort */
   sorter s(_data, start, end, f);
   s.run();
}

/*
 * Sort the specified subarray within the array.
 * Return an index array mapping sorted position --> original position
 * (within the subarray, so index 0 corresponds to start of subarray).
 */
template <typename T, typename Syn>
array<unsigned long> array<T,Syn>::sort_subarray_idx(
   unsigned long start, 
   unsigned long end, 
   const comparable_functor<T>& f) 
{
   auto_write_lock<const Syn> wlock(*this);
   #ifdef LANG__ARRAY__CHECK_BOUNDS
      /* check array bounds */
      if (start >= _size)
         throw ex_index_out_of_bounds(
            "array start index out of bounds", 
            start
         );
      else if (end >= _size)
         throw ex_index_out_of_bounds(
            "array end index out of bounds", 
            end
         );
   #endif
   /* initialize index array */
   unsigned long idx_size = (start <= end) ? (end - start + 1) : 0;
   array<unsigned long> idx(idx_size);
   for (unsigned long n = 0; n < idx_size; n++)
      idx._data[n] = start + n;
   /* sort */
   sorter_idx s(_data, idx._data - start, start, end, f);
   s.run();
   return idx;
}

/*
 * Remove duplicate elements (as defined by the given comparison functor)
 * and place the remaining unique elements in sorted order.
 */
template <typename T, typename Syn>
void array<T,Syn>::unique(
   const comparable_functor<T>& f)
{
   auto_write_lock<const Syn> wlock(*this);
   this->make_unique(f);
}

template <typename T, typename Syn>
void array<T,Syn>::make_unique(
   const comparable_functor<T>& f, unsigned long terminator_space)
{
   if (_size > 1) {
      /* sort */
      sorter s(_data, 0, (_size - 1), f);
      s.run();
      /* count unique elements */
      unsigned long size_unique = 1;
      for (unsigned long n = 1; n < _size; n++) {
         if (f(_data[n], _data[n-1]) != 0)
            size_unique++;
      }
      /* copy unique elements */
      array<T,Syn> a(size_unique + terminator_space);
      a._data[0] = _data[0];
      for (unsigned long n_unique = 0, n = 1; n < _size; n++) {
         if (f(_data[n], a._data[n_unique]) != 0)
            a._data[++n_unique] = _data[n];
      }
      /* hide terminator */
      a._size -= terminator_space;
      /* take ownership of unique array */
      array<T,Syn>::swap(*this, a);
   }
}

/*
 * Remove duplicate elements (as defined by the given comparison functor)
 * and place the remaining unique elements in sorted order.  In addition,
 * return an index array containing the original positions of the unique
 * elements.
 */
template <typename T, typename Syn>
array<unsigned long> array<T,Syn>::unique_idx(
   const comparable_functor<T>& f)
{
   auto_write_lock<const Syn> wlock(*this);
   return this->make_unique_idx(f);
}

template <typename T, typename Syn>
array<unsigned long> array<T,Syn>::make_unique_idx(
   const comparable_functor<T>& f, unsigned long terminator_space)
{
   if (_size > 1) {
      /* initialize index array */
      array<unsigned long> idx(_size);
      for (unsigned long n = 0; n < _size; n++)
         idx._data[n] = n;
      /* sort */
      sorter_idx s(_data, idx._data, 0, (_size - 1), f);
      s.run();
      /* count unique elements */
      unsigned long size_unique = 1;
      for (unsigned long n = 1; n < _size; n++) {
         if (f(_data[n], _data[n-1]) != 0)
            size_unique++;
      }
      /* copy unique elements */
      array<T,Syn>         a(size_unique + terminator_space);
      array<unsigned long> a_idx(size_unique);
      a._data[0]     = _data[0];
      a_idx._data[0] = idx._data[0];
      for (unsigned long n_unique = 0, n = 1; n < _size; n++) {
         if (f(_data[n], a._data[n_unique]) != 0) {
            ++n_unique;
            a._data[n_unique]     = _data[n];
            a_idx._data[n_unique] = idx._data[n];
         }
      }
      /* hide terminator */
      a._size -= terminator_space;
      /* take ownership of unique array */
      array<T,Syn>::swap(*this, a);
      return a_idx;
   } else {
      return array<unsigned long>(_size);
   }
}

/***************************************************************************
 * Array sorter implementation.
 ***************************************************************************/

/*
 * Constructor.
 */
template <typename T, typename Syn>
array<T,Syn>::sorter::sorter(
   T* data,
   unsigned long low,
   unsigned long high,
   const comparable_functor<T>& f)
 : _data(data),
   _low(low),
   _high(high),
   _f_compare(f)
{ }

/*
 * Destructor.
 */
template <typename T, typename Syn>
array<T,Syn>::sorter::~sorter() {
   /* do nothing */
}

/*
 * Run the sort.
 */
template <typename T, typename Syn>
void array<T,Syn>::sorter::run() {
   /* recursively sort each half of the array */
   if (_low < _high) {
      /* partition array into two halves */
      unsigned long pivot = (_low + _high)/2;
      unsigned long i = _low;
      unsigned long j = _high;
      while (true) {
         while (_f_compare(_data[i], _data[pivot]) < 0)
            i++;
         while (_f_compare(_data[j], _data[pivot]) > 0)
            j--;
         if (i < j) {
            /* swap elements */
            T temp = _data[i];
            _data[i] = _data[j];
            _data[j] = temp;
            /* update pivot index */
            if (pivot == i)
               pivot = j;
            else if (pivot == j)
               pivot = i;
            /* move inward */
            i++;
            j--;
         } else {
            break;
         }
      }
      /* sort each half */
      sorter s0(_data, _low,     j, _f_compare);
      sorter s1(_data,  j+1, _high, _f_compare);
      child_thread::run(s0, s1);
   }
}

/***************************************************************************
 * Array sorter_idx implementation.
 ***************************************************************************/

/*
 * Constructor.
 */
template <typename T, typename Syn>
array<T,Syn>::sorter_idx::sorter_idx(
   T* data,
   unsigned long* idx,
   unsigned long low,
   unsigned long high,
   const comparable_functor<T>& f)
 : sorter(data, low, high, f),
   _idx(idx)
{ }

/*
 * Destructor.
 */
template <typename T, typename Syn>
array<T,Syn>::sorter_idx::~sorter_idx() {
   /* do nothing */
}

/*
 * Run the sort.
 */
template <typename T, typename Syn>
void array<T,Syn>::sorter_idx::run() {
   /* recursively sort each half of the array */
   if (this->_low < this->_high) {
      /* partition array into two halves */
      unsigned long pivot = (this->_low + this->_high)/2;
      unsigned long i = this->_low;
      unsigned long j = this->_high;
      while (true) {
         while (this->_f_compare(this->_data[i], this->_data[pivot]) < 0)
            i++;
         while (this->_f_compare(this->_data[j], this->_data[pivot]) > 0)
            j--;
         if (i < j) {
            /* swap elements */
            T temp = this->_data[i];
            this->_data[i] = this->_data[j];
            this->_data[j] = temp;
            /* swap index */
            unsigned long temp_idx = _idx[i];
            _idx[i] = _idx[j];
            _idx[j] = temp_idx;
            /* update pivot index */
            if (pivot == i)
               pivot = j;
            else if (pivot == j)
               pivot = i;
            /* move inward */
            i++;
            j--;
         } else {
            break;
         }
      }
      /* sort each half */
      sorter_idx s0(this->_data, _idx, this->_low, j, this->_f_compare);
      sorter_idx s1(this->_data, _idx, j+1, this->_high, this->_f_compare);
      child_thread::run(s0, s1);
   }
}

} /* namespace lang */

#endif
