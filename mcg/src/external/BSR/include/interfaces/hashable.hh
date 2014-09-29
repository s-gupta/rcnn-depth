/*
 * Hashable interface.
 */
#ifndef INTERFACES__HASHABLE_HH
#define INTERFACES__HASHABLE_HH

namespace interfaces {

/*
 * Template interface for hashable objects.
 *
 * If a class X extends hashable<T,H> then:
 * (1) X must be a T, and 
 * (2) X must hash to an object of type H via its hash() method.
 *
 * Note that the interface is the same regardless of T's const qualifier.
 */
template <typename T, typename H>
class hashable {
public:
   /*
    * Destructor.
    */
   virtual ~hashable();
   
   /* 
    * Function for hashing a hashable object.
    */
   static H& hash(const hashable<T,H>&);
   
   /*
    * Hash method.
    */
   virtual H& hash() const = 0;
};

template <typename T, typename H>
class hashable<const T,H> : public hashable<T,H> { };

/*
 * Destructor.
 */
template <typename T, typename H>
hashable<T,H>::~hashable() { }

/* 
 * Function for hashing a hashable object.
 */
template <typename T, typename H>
H& hashable<T,H>::hash(const hashable<T,H>& t) {
   return (t.hash());
}

} /* namespace interfaces */

#endif
