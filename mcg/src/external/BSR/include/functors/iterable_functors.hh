/*
 * Iterable functors.
 */
#ifndef FUNCTORS__ITERABLE_FUNCTORS_HH
#define FUNCTORS__ITERABLE_FUNCTORS_HH

#include "interfaces/iterable.hh"

namespace functors {
/*
 * Imports.
 */
using interfaces::iterable;
   
/*
 * Abstract functor that can be used as an iterator.
 */
template <typename T>
class iterable_functor {
public:
   virtual ~iterable_functor() { }
   virtual void operator()(T&) const = 0;
};

/*
 * Abstract functor that can be used as an iterator 
 * over const items (and hence, also non-const items).
 */
template <typename T>
class iterable_functor<const T> : public iterable_functor<T> {
public:
   void operator()(T&) const;
   virtual void operator()(const T&) const = 0;
};

template <typename T>
void iterable_functor<const T>::operator()(T& t) const {
   const T& t_const = t;
   this->operator()(t_const);
}

/*
 * Iteration functor defined by an iterable class.
 */
template <typename T>
class iter_functor : public iterable_functor<T> {
public:
   /*
    * Constructor.
    * Initialize functor to use given iterable object.
    */
   iter_functor(const iterable<T>&);

   /*
    * Iterate method.
    * Call the underlying iter_apply method.
    */
   using iterable_functor<T>::operator();
   void operator()(T&) const;

protected:
   /* object to use for iter_apply method */
   const iterable<T>& _f_obj;
};

/*
 * Constructor.
 * Initialize functor to use given iterable object.
 */
template <typename T>
iter_functor<T>::iter_functor(const iterable<T>& f_obj)
 : _f_obj(f_obj)
{ }

/*
 * Iterate method.
 * Call the underlying iter_apply method.
 */
template <typename T>
void iter_functor<T>::operator()(T& t) const {
   _f_obj.iter_apply(t);
}

/*
 * Iteration functor for deleting objects.
 */
template <typename T>
class iter_functor_delete : public iterable_functor<T> {
public:
   using iterable_functor<T>::operator();
   void operator()(T& t) const {
      delete &t;
   }
};

/*
 * Iteration functor for deleting const objects.
 * Since const objects can't be deleted, this functor does nothing.
 */
template <typename T>
class iter_functor_delete<const T> : public iterable_functor<const T> {
public:
   using iterable_functor<const T>::operator();
   void operator()(const T& t) const { }
};

/*
 * Iteration functor for performing a no-op.
 */
template <typename T>
class iter_functor_nop : public iterable_functor<T> {
public:
   using iterable_functor<T>::operator();
   void operator()(T& t) const { }
};

/*
 * Globally accessible set of useful iteration functors.
 */
template <typename T>
class iter_functors {
public:
   static const iter_functor_delete<T>& f_delete();
   static const iter_functor_nop<T>&    f_nop();
};

template <typename T>
const iter_functor_delete<T>& iter_functors<T>::f_delete() {
   static const iter_functor_delete<T>* f = new iter_functor_delete<T>();
   return *f;
}

template <typename T>
const iter_functor_nop<T>& iter_functors<T>::f_nop() {
   static const iter_functor_nop<T>* f = new iter_functor_nop<T>();
   return *f;
}

} /* namespace functors */

#endif
