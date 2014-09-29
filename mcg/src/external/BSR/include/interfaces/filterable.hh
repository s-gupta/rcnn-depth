/*
 * Filterable interface.
 */
#ifndef INTERFACES__FILTERABLE_HH
#define INTERFACES__FILTERABLE_HH

namespace interfaces {

/*
 * Template interface for filterables.
 *
 * If a class entends filterable<T> then it must define a method (filter_apply)
 * that can be used as a filter over collections of objects of type T.
 */
template <typename T>
class filterable {
public:
   /*
    * Destructor.
    */
   virtual ~filterable();
   
   /*
    * Filter method.
    */
   virtual bool filter_apply(T&) const = 0;
};

template <typename T>
filterable<T>::~filterable() { }

/*
 * A filterable over const T also serves as a filterable over T.
 */
template <typename T>
class filterable<const T> : public filterable<T> {
public:
   /*
    * Filter method.
    * Call the filter version for const items.
    */
   bool filter_apply(T&) const;

   /*
    * Filter method.
    */
   virtual bool filter_apply(const T&) const = 0;
};

template <typename T>
bool filterable<const T>::filter_apply(T& t) const {
   const T& t_const = t;
   return (this->filter_apply(t_const));
}

} /* namespace interfaces */

#endif
