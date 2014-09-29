/*
 * Abstract collection.
 */
#ifndef COLLECTIONS__ABSTRACT__COLLECTION_HH
#define COLLECTIONS__ABSTRACT__COLLECTION_HH

#include "functors/filterable_functors.hh"
#include "functors/foldable_functors.hh"
#include "functors/iterable_functors.hh"
#include "functors/mappable_functors.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"

namespace collections {
namespace abstract {
/*
 * Imports.
 */
using functors::filterable_functor;
using functors::foldable_functor;
using functors::iterable_functor;
using functors::mappable_functor;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

/*
 * Abstract base class for collections.
 */
template <typename T>
class collection {
public:
   /*
    * Destructor.
    */
   virtual ~collection() = 0;

   /*
    * Add element(s) to the collection.
    * Return a reference to the collection.
    */
   virtual collection<T>& add(T&) = 0;
   virtual collection<T>& add(const collection<T>&) = 0;

   /*
    * Size.
    */
   bool is_empty() const;
   virtual unsigned long size() const = 0;

   /*
    * Return iterator over elements.
    */
   virtual auto_ptr< iterator<T> > iter_create() const = 0;

   /*
    * Iterate a functor over the collection.
    */
   void iter(const iterable_functor<T>&) const;

   /*
    * Apply a filter to the collection.
    * Add elements passing the filter to the specified collection.
    */
   void filter(const filterable_functor<T>&, collection<T>&) const;

   /*
    * Check if some/all of the elements pass a filter.
    */
   bool exists(const filterable_functor<T>&) const;
   bool for_all(const filterable_functor<T>&) const;

   /*
    * Apply a fold functor to the collection.
    */
   template <typename U> U& fold(const foldable_functor<T,U>&, U&) const;

   /*
    * Apply a map to elements of the collection.
    * Add each of the mapped elements to the specified collection.
    */
   template <typename U> void map(
      const mappable_functor<T,U>&,
      collection<U>&
   ) const;
};

/*
 * Pure virtual destructor.
 */
template <typename T>
collection<T>::~collection() { }

/*
 * Check if collection is empty (size zero).
 */
template <typename T>
bool collection<T>::is_empty() const {
   return (this->size() == 0);
}

/*
 * Iterate a functor over the collection.
 */
template <typename T>
void collection<T>::iter(const iterable_functor<T>& f) const {
   auto_ptr< iterator<T> > i = this->iter_create();
   while (i->has_next())
      f(i->next());
}

/*
 * Apply a filter to the collection.
 * Add elements passing the filter to the specified collection.
 */
template <typename T>
void collection<T>::filter(
   const filterable_functor<T>& f,
   collection<T>& c) const
{
   auto_ptr< iterator<T> > i = this->iter_create();
   while (i->has_next()) {
      T& t = i->next();
      if (f(t))
         c.add(t);
   }
}

/*
 * Check if at least one of the elements passes a filter.
 */
template <typename T>
bool collection<T>::exists(const filterable_functor<T>& f) const {
   bool found = false;
   auto_ptr< iterator<T> > i = this->iter_create();
   while ((i->has_next()) && (!found))
      found = f(i->next());
   return found;
}

/*
 * Check if all of the elements pass a filter.
 */
template <typename T>
bool collection<T>::for_all(const filterable_functor<T>& f) const {
   bool flag = true;
   auto_ptr< iterator<T> > i = this->iter_create();
   while ((i->has_next()) && (flag))
      flag = f(i->next());
   return flag;
}

/*
 * Apply a fold functor to the collection.
 */
template <typename T>
template <typename U>
U& collection<T>::fold(const foldable_functor<T,U>& f, U& u) const {
   U* u_result = &u;
   auto_ptr< iterator<T> > i = this->iter_create();
   while (i->has_next()) {
      U& u_temp = f(i->next(), *u_result);
      u_result = &u_temp;
   }
   return *u_result;
}
   
/*
 * Apply a map to elements of the collection.
 * Add each of the mapped elements to the specified collection.
 */
template <typename T>
template <typename U>
void collection<T>::map(
   const mappable_functor<T,U>& f,
   collection<U>& c) const
{
   auto_ptr< iterator<T> > i = this->iter_create();
   while (i->has_next())
      c.add(f(i->next()));
}

} /* namespace abstract */
} /* namespace collections */

#endif
