/*
 * Mappable functors.
 */
#ifndef FUNCTORS__MAPPABLE_FUNCTORS_HH
#define FUNCTORS__MAPPABLE_FUNCTORS_HH

#include "interfaces/mappable.hh"

namespace functors {
/*
 * Imports.
 */
using interfaces::mappable;

/*
 * Abstract functor that can be used as a map.
 */
template <typename T, typename U>
class mappable_functor {
public:
   virtual ~mappable_functor() { }
   virtual U& operator()(T&) const = 0;
};

/*
 * Abstract functor that can be used as a map over
 * const items (and hence, also non-const items).
 */
template <typename T, typename U>
class mappable_functor<const T, U> : public mappable_functor<T,U> {
public:
   U& operator()(T&) const;
   virtual U& operator()(const T&) const = 0;
};

template <typename T, typename U>
U& mappable_functor<const T, U>::operator()(T& t) const {
   const T& t_const = t;
   return this->operator()(t_const);
}

/*
 * Map functor defined by a mappable class.
 */
template <typename T, typename U>
class map_functor : public mappable_functor<T,U> {
public:
   /*
    * Constructor.
    * Initialize functor to use given mappable object.
    */
   map_functor(const mappable<T,U>&);

   /*
    * Map method.
    * Call the underlying map_apply method.
    */
   using mappable_functor<T,U>::operator();
   U& operator()(T&) const;

protected:
   /* object to use for map_apply method */
   const mappable<T,U>& _f_obj;
};

/*
 * Constructor.
 * Initialize functor to use given mappable object.
 */
template <typename T, typename U>
map_functor<T,U>::map_functor(const mappable<T,U>& f_obj) : _f_obj(f_obj) { }

/*
 * Map method.
 * Call the underlying map_apply method.
 */
template <typename T, typename U>
U& map_functor<T,U>::operator()(T& t) const {
   return _f_obj.map_apply(t);
}

/*
 * Identity map functor.
 */
template <typename T, typename U = T>
class map_functor_identity : public mappable_functor<T,U> {
public:
   using mappable_functor<T,U>::operator();
   U& operator()(T& t) const {
      return t;
   }
};

/*
 * Map functor for copying items.
 */
template <typename T, typename U = T>
class map_functor_copy : public mappable_functor<T,U> {
public:
   using mappable_functor<T,U>::operator();
   U& operator()(T& t) const {
      return *(new U(t));
   }
};

/*
 * Globally accessible set of useful map functors.
 */
template <typename T, typename U = T>
class map_functors {
public:
   static const map_functor_identity<T,U>& f_identity();
   static const map_functor_copy<T,U>&     f_copy();
};

template <typename T, typename U>
const map_functor_identity<T,U>& map_functors<T,U>::f_identity() {
   static const map_functor_identity<T,U>* f = new map_functor_identity<T,U>();
   return *f;
}

template <typename T, typename U>
const map_functor_copy<T,U>& map_functors<T,U>::f_copy() {
   static const map_functor_copy<T,U>* f = new map_functor_copy<T,U>();
   return *f;
}

} /* namespace functors */

#endif
