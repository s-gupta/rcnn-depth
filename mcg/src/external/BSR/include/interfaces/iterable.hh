/*
 * Iterable interface.
 */
#ifndef INTERFACES__ITERABLE_HH
#define INTERFACES__ITERABLE_HH

namespace interfaces {

/*
 * Template interface for iterables.
 *
 * If a class extends iterable<T> then it must define a method (iter_apply)
 * that can be used as an iterator over collections of objects of type T.
 */
template <typename T>
class iterable {
public:
   /*
    * Destructor.
    */
   virtual ~iterable();

   /*
    * Iterate method.
    */
   virtual void iter_apply(T&) const = 0;
};

template <typename T>
iterable<T>::~iterable() { }

/*
 * An interable over const T also serves as an interable over T.
 */
template <typename T>
class iterable<const T> : public iterable<T> {
public:
   /*
    * Iterate method.
    * Call the iter_apply version for const items.
    */
   void iter_apply(T&) const;

   /*
    * Iterate method.
    */
   virtual void iter_apply(const T&) const = 0;
};

template <typename T>
void iterable<const T>::iter_apply(T& t) const {
   const T& t_const = t;
   this->iter_apply(t_const);
}

} /* namespace interfaces */

#endif
