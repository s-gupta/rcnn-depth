/*
 * Filterable functors.
 */
#ifndef FUNCTORS__FILTERABLE_FUNCTORS_HH
#define FUNCTORS__FILTERABLE_FUNCTORS_HH

#include "interfaces/filterable.hh"

namespace functors {
/*
 * Imports.
 */
using interfaces::filterable;

/*
 * Abstract functor that can be used as a filter.
 */
template <typename T>
class filterable_functor {
public:
   virtual ~filterable_functor() { }
   virtual bool operator()(T&) const = 0;
};

/*
 * Abstract functor that can be used as a filter 
 * over const items (and hence, also non-const items).
 */
template <typename T>
class filterable_functor<const T> : public filterable_functor<T> {
public:
   bool operator()(T&) const;
   virtual bool operator()(const T&) const = 0;
};

template <typename T>
bool filterable_functor<const T>::operator()(T& t) const {
   const T& t_const = t;
   return this->operator()(t_const);
}

/*
 * Filter functor defined by a filterable class.
 */
template <typename T>
class filter_functor : public filterable_functor<T> {
public:
   /*
    * Constructor.
    * Initialize functor to use the given filterable object.
    */
   filter_functor(const filterable<T>&);

   /*
    * Filter method.
    * Call the underlying filter_apply method.
    */
   using filterable_functor<T>::operator();
   bool operator()(T&) const;

protected:
   /* object to use for filter_apply method */
   const filterable<T>& _f_obj;
};

/*
 * Constructor.
 * Initialize functor to use the given filterable object.
 */
template <typename T>
filter_functor<T>::filter_functor(const filterable<T>& f_obj)
 : _f_obj(f_obj) 
{ }

/*
 * Filter method.
 * Call the underlying filter_apply method.
 */
template <typename T>
bool filter_functor<T>::operator()(T& t) const {
   return _f_obj.filter_apply(t);
}

/*
 * Filter functor that always returns returns true.
 */
template <typename T>
class filter_functor_true : public filterable_functor<T> {
public:
   using filterable_functor<T>::operator();
   bool operator()(T& t) const {
      return true;
   }
};

/*
 * Filter functor that always returns false.
 */
template <typename T>
class filter_functor_false : public filterable_functor<T> {
public:
   using filterable_functor<T>::operator();
   bool operator()(T& t) const {
      return false;
   }
};

/*
 * Globally accessible set of useful filter functors.
 */
template <typename T>
class filter_functors {
public:
   static const filter_functor_true<T>&  f_true();
   static const filter_functor_false<T>& f_false();
};

template <typename T>
const filter_functor_true<T>& filter_functors<T>::f_true() {
   static const filter_functor_true<T>* f = new filter_functor_true<T>();
   return *f;
}

template <typename T>
const filter_functor_false<T>& filter_functors<T>::f_false() {
   static const filter_functor_false<T>* f = new filter_functor_false<T>();
   return *f;
}

} /* namespace functors */

#endif
