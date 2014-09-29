/*
 * Equalable interface.
 */
#ifndef INTERFACES__EQUALABLE_HH
#define INTERFACES__EQUALABLE_HH

namespace interfaces {

/*
 * Template interface for equalables.
 *
 * If class X extends equalable<T>, then:
 * (1) X must be a T, and
 * (2) X must implement equality testing across type T.
 *
 * Note that the interface is the same regardless of the type const qualifier.
 */
template <typename T>
class equalable {
public:
   /*
    * Destructor.
    */
   virtual ~equalable();
   
   /*
    * Function for testing equality of two equalable objects.
    */
   static bool is_equal(const equalable<T>&, const equalable<T>&);
   
   /*
    * Method for testing equality with other equalable<T> objects.
    * A equalable<T> is cast to a T before calling is_equal_to.
    */
   bool is_equal_to(const equalable<T>&) const;

   /*
    * Virtual equality test method.
    */
   virtual bool is_equal_to(const T&) const = 0;
};

template <typename T>
class equalable<const T> : public equalable<T> { };

/*
 * Destructor.
 */
template <typename T>
equalable<T>::~equalable() { }

/*
 * Function for testing equality of two equalable objects.
 */
template <typename T>
bool equalable<T>::is_equal(const equalable<T>& t0, const equalable<T>& t1) {
   return (t0.is_equal_to(static_cast<const T&>(t1)));
}

/*
 * Method for testing equality with other equalable<T> objects.
 * A equalable<T> is cast to a T before calling is_equal_to.
 */
template <typename T>
bool equalable<T>::is_equal_to(const equalable<T>& that) const {
   return (this->is_equal_to(static_cast<const T&>(that)));
}

} /* namespace interfaces */

#endif
