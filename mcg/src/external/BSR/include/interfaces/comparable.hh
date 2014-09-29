/*
 * Comparable interface.
 */
#ifndef INTERFACES__COMPARABLE_HH
#define INTERFACES__COMPARABLE_HH

#include "interfaces/equalable.hh"

namespace interfaces {

/*
 * Template interface for comparables.
 *
 * If class X extends comparable<T>, then:
 * (1) X must be a T, and
 * (2) X must implement comparisons across type T.
 *
 * Note that the interface is the same regardless of the type const qualifier.
 */
template <typename T>
class comparable : public equalable<const T> {
public:
   /*
    * Destructor.
    */
   virtual ~comparable();
   
   /*
    * Function for comparing two comparable objects.
    */
   static int compare(const comparable<T>&, const comparable<T>&);
   
   /*
    * Method for comparing to other comparable<T> objects.
    * A comparable<T> is cast to a T before calling compare_to.
    */
   int compare_to(const comparable<T>&) const;

   /*
    * Virtual comparison method.
    */
   virtual int compare_to(const T&) const = 0;

   /*
    * Equality test method.
    */
   virtual bool is_equal_to(const T&) const;
};

template <typename T>
class comparable<const T> : public comparable<T> { };

/*
 * Destructor.
 */
template <typename T>
comparable<T>::~comparable() { }

/*
 * Function for comparing two comparable objects.
 */
template <typename T>
int comparable<T>::compare(const comparable<T>& t0, const comparable<T>& t1) {
   return (t0.compare_to(static_cast<const T&>(t1)));
}

/*
 * Method for comparing to other comparable<T> objects.
 * A comparable<T> is cast to a T before calling compare_to.
 */
template <typename T>
int comparable<T>::compare_to(const comparable<T>& that) const {
   return (this->compare_to(static_cast<const T&>(that)));
}

/*
 * Equality test method.
 */
template <typename T>
bool comparable<T>::is_equal_to(const T& that) const {
   return (this->compare_to(that) == 0);
}

} /* namespace interfaces */

#endif
