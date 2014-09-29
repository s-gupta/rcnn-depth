/*
 * Distanceable interface.
 */
#ifndef INTERFACES__DISTANCEABLE_HH
#define INTERFACES__DISTANCEABLE_HH

namespace interfaces {

/*
 * Template interface for distanceables.
 *
 * If class X extends distanceable<T,V>, then:
 * (1) X must be a T, and
 * (2) X must implement a distance function across type T 
 *     that returns a value of type V.
 *
 * Note that the interface is the same regardless of T's const qualifier.
 */
template <typename T, typename V>
class distanceable {
public:
   /*
    * Destructor.
    */
   virtual ~distanceable();

   /*
    * Function for computing distance between two distanceable objects.
    */
   static V distance(const distanceable<T,V>&, const distanceable<T,V>&);

   /*
    * Method for computing distance to other distanceable<T,V> objects.
    * A distanceable<T,V> is cast to a T before calling distance_to.
    */
   V distance_to(const distanceable<T,V>&) const;

   /*
    * Virtual distance method.
    */
   virtual V distance_to(const T&) const = 0;
};

template <typename T, typename V>
class distanceable<const T,V> : public distanceable<T,V> { };

/*
 * Destructor.
 */
template <typename T, typename V>
distanceable<T,V>::~distanceable() { }

/*
 * Function for computing distance between two distanceable objects.
 */
template <typename T, typename V>
V distanceable<T,V>::distance(
   const distanceable<T,V>& t0,
   const distanceable<T,V>& t1)
{
   return (t0.distance_to(static_cast<const T&>(t1)));
}

/*
 * Method for computing distance to other distanceable<T,V> objects.
 * A distanceable<T,V> is cast to a T before calling distance_to.
 */
template <typename T, typename V>
V distanceable<T,V>::distance_to(const distanceable<T,V>& that) const {
   return (this->distance_to(static_cast<const T&>(that)));
}

} /* namespace interfaces */

#endif
