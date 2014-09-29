/*
 * Abstract array.
 */
#ifndef COLLECTIONS__ABSTRACT__ARRAY_HH
#define COLLECTIONS__ABSTRACT__ARRAY_HH

#include "collections/abstract/collection.hh"
#include "functors/comparable_functors.hh"
#include "lang/array.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/random/generators/rand_gen.hh"
#include "math/random/generators/rand_gen_uniform.hh"

namespace collections {
namespace abstract {
/*
 * Imports.
 */
using functors::comparable_functor;
using functors::compare_functors;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;
using math::random::generators::rand_gen;
using math::random::generators::rand_gen_uniform;

/*
 * Abstract base class for arrays.
 */
template <typename T>
class array : virtual public collection<T> {
public:
   /*
    * Destructor.
    */
   virtual ~array() = 0;

   /*
    * Element reference.
    */
   T& operator()(unsigned long) const;
   virtual T& operator[](unsigned long) const = 0;
   
   /*
    * Element replacement.
    * Return a reference to the element that was replaced.
    */
   virtual T& replace(unsigned long, T&) = 0;
   
   /*
    * Subarray operations.
    * Add elements in the specified range to the given collection.
    */
   virtual void subarray(
      unsigned long,                      /* start index */
      unsigned long,                      /* end index   */
      collection<T>&
   ) const = 0;

   virtual void subarray(
      unsigned long,                      /* start index */ 
      unsigned long,                      /* step size   */
      unsigned long,                      /* end index   */
      collection<T>&
   ) const = 0;
   
   virtual void subarray(
      const lang::array<unsigned long>&,  /* indices of desired elements */
      collection<T>&
   ) const = 0;

   /*
    * Collection interface for adding element(s) to array.
    * Return a reference to the array.
    */
   virtual array<T>& add(T&) = 0;
   virtual array<T>& add(const collection<T>&) = 0;
    
   /*
    * Reverse array.
    * Return a reference to the array.
    */
   virtual array<T>& reverse() = 0;

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
    * Randomly permute the elements of the array.
    * Return an index array mapping resulting position --> original position.
    *
    * The generator from which to draw values used in computing the random
    * permutation may optionally be specified.  If unspecified, a generator 
    * is created and used for this purpose.
    */
   lang::array<unsigned long> randperm();
 
   virtual lang::array<unsigned long> randperm(rand_gen<>&) = 0;

   /*
    * Randomly permute the elements of the specified subarray within the array.
    * Return an index array mapping resulting position --> original position
    * (within the subarray, so index 0 corresponds to start of subarray).
    */
   lang::array<unsigned long> randperm_subarray(
      unsigned long,    /* start index */ 
      unsigned long     /* end index   */
   );

   virtual lang::array<unsigned long> randperm_subarray(
      unsigned long,    /* start index */ 
      unsigned long,    /* end index   */
      rand_gen<>&       /* generator   */
   ) = 0;

   /*
    * Sort elements in ascending order according to the given comparison
    * functor.
    */
   virtual void sort(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   ) = 0;

   /*
    * Sort the array.
    * Return an index array mapping sorted position --> original position.
    */
   virtual lang::array<unsigned long> sort_idx(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   ) = 0;
  
   /*
    * Sort the specified subarray within the array.
    */
   virtual void sort_subarray(
      unsigned long,    /* start index */ 
      unsigned long,    /* end index   */
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   ) = 0;

   /*
    * Sort the specified subarray within the array.
    * Return an index array mapping sorted position --> original position
    * (within the subarray, so index 0 corresponds to start of subarray).
    */
   virtual lang::array<unsigned long> sort_subarray_idx(
      unsigned long,    /* start index */ 
      unsigned long,    /* end index   */
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   ) = 0;

   /*
    * Remove duplicate elements (as defined by the given comparison functor)
    * and place the remaining unique elements in sorted order.
    */
   virtual void unique(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   ) = 0;

   /*
    * Remove duplicate elements (as defined by the given comparison functor)
    * and place the remaining unique elements in sorted order.  In addition,
    * return an index array containing the original positions of the unique
    * elements.
    */
   virtual lang::array<unsigned long> unique_idx(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   ) = 0;
};

/*
 * Pure virtual destructor.
 */
template <typename T>
array<T>::~array() { }

/*
 * Element reference (functor form).
 */
template <typename T>
T& array<T>::operator()(unsigned long n) const {
   return this->operator[](n);
}

/*
 * Randomly permute the elements of the array.
 * Return an index array mapping resulting position --> original position.
 */
template <typename T>
lang::array<unsigned long> array<T>::randperm() {
   rand_gen_uniform<> r;
   return this->randperm(r);
}

/*
 * Randomly permute the elements of the specified subarray within the array.
 * Return an index array mapping resulting position --> original position
 * (within the subarray, so index 0 corresponds to start of subarray).
 */
template <typename T>
lang::array<unsigned long> array<T>::randperm_subarray(
   unsigned long start,
   unsigned long end)
{
   rand_gen_uniform<> r;
   return this->randperm_subarray(start, end, r);
}

} /* namespace abstract */
} /* namespace collections */

#endif
