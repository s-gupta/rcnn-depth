/*
 * Foldable functors.
 */
#ifndef FUNCTORS__FOLDABLE_FUNCTORS_HH
#define FUNCTORS__FOLDABLE_FUNCTORS_HH

#include "interfaces/foldable.hh"

namespace functors {
/*
 * Imports.
 */
using interfaces::foldable;

/*
 * Abstract fold functor.
 */
template <typename T, typename U>
class foldable_functor {
public:
   virtual ~foldable_functor() { }
   virtual U& operator()(T&, U&) const = 0;
};

/*
 * Abstract fold functor for use with const items.
 */
template <typename T, typename U>
class foldable_functor<const T, U> : public foldable_functor<T,U> {
public:
   U& operator()(T&, U&) const;
   virtual U& operator()(const T&, U&) const = 0;
};

template <typename T, typename U>
U& foldable_functor<const T, U>::operator()(T& t, U& u) const {
   const T& t_const = t;
   return this->operator()(t_const, u);
}

/*
 * Fold functor defined by a foldable class.
 */
template <typename T, typename U>
class fold_functor : public foldable_functor<T,U> {
public:
   /*
    * Constructor.
    * Initialize functor to use the given foldable object.
    */
   fold_functor(const foldable<T,U>&);

   /*
    * Fold method.
    * Call the underlying fold_apply method.
    */
   using foldable_functor<T,U>::operator();
   U& operator()(T&, U&) const;

protected:
   /* object to use for fold_apply method */
   const foldable<T,U>& _f_obj;
};

/*
 * Constructor.
 * Initialize functor to use the given foldable object.
 */
template <typename T, typename U>
fold_functor<T,U>::fold_functor(const foldable<T,U>& f_obj)
 : _f_obj(f_obj) 
{ }

/*
 * Fold method.
 * Call the underlying fold_apply method.
 */
template <typename T, typename U>
U& fold_functor<T,U>::operator()(T& t, U& u) const {
   return _f_obj.fold_apply(t, u);
}

} /* namespace functors */

#endif
