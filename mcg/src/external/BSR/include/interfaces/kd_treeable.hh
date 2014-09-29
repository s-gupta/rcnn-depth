/*
 * kd-treeable interface.
 */
#ifndef INTERFACES__KD_TREEABLE_HH
#define INTERFACES__KD_TREEABLE_HH

#include "collections/abstract/collection.hh"
#include "interfaces/distanceable.hh"
#include "interfaces/keyable.hh"

namespace interfaces {
/*
 * Imports.
 */
using collections::abstract::collection;

/*
 * kd-tree keys.
 * 
 * A kd-tree key stores a dimension and value along that 
 * dimension at which data is split into two groups.
 *
 * V should be a numeric type for which the standard 
 * comparison and arithmetic operators are available.
 */
template <typename V = double>
class kd_tree_key {
public:
   /*
    * Constructor.
    * Initialize split dimension and value.
    */
   kd_tree_key(unsigned long d, V v)
     : _dimension(d), _value(v) { }

   /*
    * Copy constructor.
    */
   kd_tree_key(const kd_tree_key<V>& k)
     : _dimension(k._dimension), _value(k._value) { }

   /*
    * Return split dimension/value.
    */
   unsigned long split_dimension() const { return _dimension; }
   V split_value() const { return _value; }

protected:
   unsigned long _dimension;  /* split dimension */
   V _value;                  /* split value */
};

/*
 * Template interface for objects that are kd-treeable
 * (in other words, for which kd-trees can be built).
 *
 * kd-trees can be built for objects extending this class 
 * without specifying extra functors at construction time.
 *
 * If a class X extends kd_treeable<T,V> then:
 * (1) X must be a T, and 
 * (2) kd-trees (of dimensionality given by dimensionality(...))
 *     containing X's must be constructable with interior 
 *     key nodes of type kd_tree_key<V>.
 */
template <typename T, typename V = double>
class kd_treeable : public keyable< T, kd_tree_key<V> >,
                    public distanceable<const T,V> {
public:
   /*
    * Destructor.
    */
   virtual ~kd_treeable();

   /*
    * Dimensionality of each object in the collection when placed in a kd-tree.
    * The implementation should be in T::dimensionality(...).
    */
   static unsigned long dimensionality(const collection<T>&);
   
   /*
    * Function for computing distance from an object to a key.
    */
   static V distance(const kd_treeable<T,V>&, const kd_tree_key<V>&);
   
   /*
    * Methods for computing distance between objects.
    */
   using distanceable<const T,V>::distance_to;

   /*
    * Virtual method for computing distance to a key.
    */
   virtual V distance_to(const kd_tree_key<V>&) const = 0;
};

/*
 * Destructor.
 */
template <typename T, typename V>
kd_treeable<T,V>::~kd_treeable() { }

/*
 * Dimensionality of each object in the collection when placed in a kd-tree.
 * The implementation should be in T::dimensionality(...).
 */
template <typename T, typename V>
unsigned long kd_treeable<T,V>::dimensionality(const collection<T>& c) {
   return T::dimensionality(c);
}

/*
 * Function for computing distance from an object to a key.
 */
template <typename T, typename V>
V kd_treeable<T,V>::distance(
   const kd_treeable<T,V>& t,
   const kd_tree_key<V>& k)
{
   return (t.distance_to(k));
}

} /* namespace interfaces */

#endif
