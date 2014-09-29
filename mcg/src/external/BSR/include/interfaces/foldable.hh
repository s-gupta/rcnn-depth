/*
 * Foldable interface.
 */
#ifndef INTERFACES__FOLDABLE_HH
#define INTERFACES__FOLDABLE_HH

namespace interfaces {

/*
 * Template interface for foldables.
 *
 * If a class extends foldable<T,U> then it must define a method (fold_apply) 
 * that can be used as a fold function over collections of objects of type T,
 * starting with a given object of type U.
 */
template <typename T, typename U>
class foldable {
public:
   /*
    * Destructor.
    */
   virtual ~foldable();

   /*
    * Fold method.
    */
   virtual U& fold_apply(T&, U&) const = 0;
};

template <typename T, typename U>
foldable<T,U>::~foldable() { }

/*
 * A foldable over const T also serves as a foldable over T.
 */
template <typename T, typename U>
class foldable<const T, U> : public foldable<T,U> {
public:
   /*
    * Fold method.
    * Call the fold_apply version for const items.
    */
   U& fold_apply(T&) const;

   /*
    * Fold method.
    */
   virtual U& fold_apply(const T&, U&) const = 0;
};

template <typename T, typename U>
U& foldable<const T, U>::fold_apply(T& t) const {
   const T& t_const = t;
   return (this->fold_apply(t_const));
}

} /* namespace interfaces */

#endif
